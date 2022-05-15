/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 *
 * \brief Defines the MD Graph class
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 *
 * \ingroup module_mdlib
 */

#define CU_INIT_UUID_STATIC
#include "mdgraph_gpu_impl.h"

#include <stdio.h>

#include <cuda_etbl/cuda_graphs.h>
#include <cuda_etbl/graph_update.h>

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/utility/gmxmpi.h"

const bool useGraph = (getenv("GMX_CUDA_GRAPH") != nullptr);

namespace gmx
{

MdGpuGraph::Impl::Impl(const DeviceStreamManager& deviceStreamManager,
                       SimulationWorkload         simulationWork,
                       MPI_Comm                   mpi_comm,
                       gmx_wallcycle*             wcycle) :
    deviceStreamManager_(deviceStreamManager),
    simulationWork_(simulationWork),
    havePPDomainDecomposition_(simulationWork.havePpDomainDecomposition),
    mpi_comm_(mpi_comm),
    wcycle_(wcycle)
{
    event_ = new GpuEventSynchronizer();

    if (havePPDomainDecomposition_)
    {
        MPI_Barrier(mpi_comm_);
        MPI_Comm_size(mpi_comm_, &ppSize_);
        MPI_Comm_rank(mpi_comm_, &ppRank_);
    }
}

MdGpuGraph::Impl::~Impl() = default;


void MdGpuGraph::Impl::enqueueEventFromAllPpRanksToRank0Stream(GpuEventSynchronizer* event,
                                                               const DeviceStream&   stream)
{

    for (int remotePpRank = 1; remotePpRank < ppSize_; remotePpRank++)
    {
        if (ppRank_ == remotePpRank)
        {
            // send event to rank 0
            MPI_Send(&event,
                     sizeof(GpuEventSynchronizer*), //NOLINT(bugprone-sizeof-expression)
                     MPI_BYTE,
                     0,
                     0,
                     mpi_comm_);
        }
        else if (ppRank_ == 0)
        {
            // rank 0 enqueues recieved event
            GpuEventSynchronizer* eventToEnqueue;
            MPI_Recv(&eventToEnqueue,
                     sizeof(GpuEventSynchronizer*), //NOLINT(bugprone-sizeof-expression)
                     MPI_BYTE,
                     remotePpRank,
                     0,
                     mpi_comm_,
                     MPI_STATUS_IGNORE);
            eventToEnqueue->enqueueWaitEvent(stream);
        }
    }

    if (ppRank_ == 0)
    {
        // rank 0 also enqueues its local event
        event->enqueueWaitEvent(stream);
    }
}

void MdGpuGraph::Impl::enqueueRank0EventToAllPpStreams(GpuEventSynchronizer* event, const DeviceStream& stream)
{

    for (int remotePpRank = 1; remotePpRank < ppSize_; remotePpRank++)
    {
        if (ppRank_ == 0)
        {
            // rank 0 sends event to remote rank
            MPI_Send(&event,
                     sizeof(GpuEventSynchronizer*), //NOLINT(bugprone-sizeof-expression)
                     MPI_BYTE,
                     remotePpRank,
                     0,
                     mpi_comm_);
        }
        else if (ppRank_ == remotePpRank)
        {
            // remote rank enqueues recieved event to its stream
            GpuEventSynchronizer* eventToEnqueue;
            MPI_Recv(&eventToEnqueue,
                     sizeof(GpuEventSynchronizer*), //NOLINT(bugprone-sizeof-expression)
                     MPI_BYTE,
                     0,
                     0,
                     mpi_comm_,
                     MPI_STATUS_IGNORE);
            eventToEnqueue->enqueueWaitEvent(stream);
        }
    }

    if (ppRank_ == 0)
    {
        // rank 0 also enqueues event to its local stream
        event->enqueueWaitEvent(stream);
    }
}

void MdGpuGraph::Impl::start(bool                  bNS,
                             bool                  canUseGraphThisStep,
                             bool                  usedGraphLastStep,
                             GpuEventSynchronizer* xReadyOnDeviceEvent)
{
    if (!useGraph)
    {
        return;
    }
    if (bNS)
    {
        graphCreated_ = false;
        return;
    }

    useGraphThisStep_ = canUseGraphThisStep;
    graphIsCapturing_ = (!graphCreated_ && useGraphThisStep_);

    if (graphIsCapturing_)
    {
        wallcycle_start(wcycle_, WallCycleCounter::LaunchGpu);
        wallcycle_sub_start(wcycle_, WallCycleSubCounter::MdGpuGraphCapture);
    }

    if (useGraphThisStep_ && !usedGraphLastStep)
    {
        // Ensure NB local stream on Rank 0 (which will be used for graph capture and/or launch)
        // waits for coordinates to be ready on all ranks
        enqueueEventFromAllPpRanksToRank0Stream(
                xReadyOnDeviceEvent, deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedLocal));
    }

    if (graphIsCapturing_)
    {
        graphCreated_ = true;

        if (ppRank_ == 0)
        {
            stat_ = cudaStreamBeginCapture(
                    deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedLocal).stream(),
                    cudaStreamCaptureModeGlobal);
            CU_RET_ERR(stat_,
                       "cudaStreamBeginCapture in MD graph definition initialization failed.");
        }

        if (havePPDomainDecomposition_)
        {
            // Fork remote NB local streams from Rank 0 NB local stream
            event_->markEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedLocal));
            enqueueRank0EventToAllPpStreams(
                    event_, deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedLocal));

            // Fork NB non-local stream from NB local stream on each rank
            event_->markEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedLocal));
            event_->enqueueWaitEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedNonLocal));
        }

        // Fork update stream from NB local stream on each rank
        event_->markEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedLocal));
        event_->enqueueWaitEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::UpdateAndConstraints));

        if (simulationWork_.useGpuPme)
        {
            // Fork PME stream from NB local stream on each rank
            event_->markEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedLocal));
            event_->enqueueWaitEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::Pme));
        }

        // Re-mark xReadyOnDeviceEvent to allow full isolation within graph capture
        xReadyOnDeviceEvent->markEvent(
                deviceStreamManager_.stream(gmx::DeviceStreamType::UpdateAndConstraints));
    }
};

static int writeFileCount;

void MdGpuGraph::Impl::end(GpuEventSynchronizer* xUpdatedOnDeviceEvent)
{

    if (!useGraphThisStep_)
    {
        return;
    }

    if (graphIsCapturing_)
    {

        if (simulationWork_.useGpuPme)
        {
            // Join PME stream to NB local stream on each rank
            event_->markEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::Pme));
            event_->enqueueWaitEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedLocal));
        }

        // Join update stream to NB local stream on each rank
        event_->markEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::UpdateAndConstraints));
        event_->enqueueWaitEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedLocal));

        if (havePPDomainDecomposition_)
        {
            // Join NB non-local stream to NB local stream on each rank
            event_->markEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedNonLocal));
            event_->enqueueWaitEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedLocal));

            // Join remote NB local streams to Rank 0 NB local stream
            event_->markEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedLocal));
            enqueueEventFromAllPpRanksToRank0Stream(
                    event_, deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedLocal));
        }

        wallcycle_sub_stop(wcycle_, WallCycleSubCounter::MdGpuGraphCapture);

        if (ppRank_ == 0)
        {
            stat_ = cudaStreamEndCapture(
                    deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedLocal).stream(), &graph_);
            CU_RET_ERR(stat_, "cudaStreamEndCapture in MD graph definition finalization failed.");

            const CUetblCudaGraphs* pETBL_CG = nullptr;
            cuGetExportTable((const void**)&pETBL_CG, &CU_ETID_CudaGraphs);
            char filename[80];
            sprintf(filename, "src.%d.dot", writeFileCount);
            pETBL_CG->GraphDebugDotPrint(graph_, filename, CU_GRAPH_DEBUG_DOT_FLAGS_VERBOSE);

            // Instantiate graph, or update a previously instantiated graph (if possible)
            bool updateFailed = false;
            if (updateGraph_)
            {
                wallcycle_sub_start(wcycle_, WallCycleSubCounter::MdGpuGraphUpdate);

                const CUetblGraphExecUpdate* pETBL = nullptr;
                cuGetExportTable((const void**)&pETBL, &CU_ETID_GraphExecUpdate);
                CUgraphExecUpdateResultInfo info;
                stat_ = (cudaError_t)pETBL->GraphExecUpdateByEdge(instance_, graph_, &info);
                // CU_RET_ERR(stat_,
                //         "cudaGraphExecUpdateByEdge in MD graph definition finalization failed.");

                sprintf(filename, "execUpdated.%d.dot", writeFileCount);
                pETBL_CG->GraphExecDebugDotPrint(instance_, filename, CU_GRAPH_DEBUG_DOT_FLAGS_VERBOSE);

                if (stat_ != 0)
                {
                    printf("cudaGraphExecUpdateByEdge failed, will fall back to "
                           "re-instantiation.\n");
                    updateFailed = true;
                }

                wallcycle_sub_stop(wcycle_, WallCycleSubCounter::MdGpuGraphUpdate);
            }
            if (!updateGraph_ || updateFailed || (getenv("GMX_NEVER_USE_UPDATED_GRAPH") != nullptr))
            {
                wallcycle_sub_start(wcycle_, WallCycleSubCounter::MdGpuGraphInstantiate);
                stat_ = cudaGraphInstantiate(&instance_, graph_, nullptr, nullptr, 0);
                wallcycle_sub_stop(wcycle_, WallCycleSubCounter::MdGpuGraphInstantiate);
                CU_RET_ERR(stat_,
                           "cudaGraphInstantiate in MD graph definition finalization failed.");

                sprintf(filename, "execInstantiated.%d.dot", writeFileCount++);
                pETBL_CG->GraphExecDebugDotPrint(instance_, filename, CU_GRAPH_DEBUG_DOT_FLAGS_VERBOSE);
            }

            // With current CUDA, only single-threaded update is possible.
            // Multi-threaded update support will be available in a future CUDA release.
            if (ppSize_ == 1 || (getenv("GMX_CUDA_GRAPH_FORCE_UPDATE") != nullptr))
            {
                updateGraph_ = true;
            }
        }
        if (havePPDomainDecomposition_)
        {
            MPI_Barrier(mpi_comm_);
        }
        wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpu);
    }

    if (ppRank_ == 0)
    {
        wallcycle_start(wcycle_, WallCycleCounter::LaunchGpu);
        wallcycle_sub_start(wcycle_, WallCycleSubCounter::MdGpuGraphLaunch);
        stat_ = cudaGraphLaunch(
                instance_, deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedLocal).stream());
        wallcycle_sub_stop(wcycle_, WallCycleSubCounter::MdGpuGraphLaunch);
        wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpu);
        CU_RET_ERR(stat_, "cudaGraphLaunch in MD graph definition finalization failed.");
        event_->markEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedLocal));
    }
    enqueueRank0EventToAllPpStreams(
            event_, deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedLocal));
    xUpdatedOnDeviceEvent->markEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedLocal));
};

MdGpuGraph::MdGpuGraph(const DeviceStreamManager& deviceStreamManager,
                       SimulationWorkload         simulationWork,
                       MPI_Comm                   mpi_comm,
                       gmx_wallcycle*             wcycle) :
    impl_(new Impl(deviceStreamManager, simulationWork, mpi_comm, wcycle))
{
}

MdGpuGraph::~MdGpuGraph() = default;

void MdGpuGraph::start(bool bNS, bool canUseGraphThisStep, bool usedGraphLastStep, GpuEventSynchronizer* xReadyOnDeviceEvent)
{
    impl_->start(bNS, canUseGraphThisStep, usedGraphLastStep, xReadyOnDeviceEvent);
}

void MdGpuGraph::end(GpuEventSynchronizer* xUpdatedOnDeviceEvent)
{
    impl_->end(xUpdatedOnDeviceEvent);
}

bool MdGpuGraph::useGraphThisStep() const
{
    return impl_->useGraphThisStep();
}

bool MdGpuGraph::graphIsCapturing() const
{
    return impl_->graphIsCapturing();
}

} // namespace gmx
