/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020,2021, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 *
 * \brief Implements backend-agnostic GPU Force Reduction functions
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include "gpuforcereduction_impl.h"

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#if GMX_GPU_CUDA
#    include "gromacs/gpu_utils/gpueventsynchronizer.cuh"
#elif GMX_GPU_SYCL
#    include "gromacs/gpu_utils/gpueventsynchronizer_sycl.h"
#endif
#include "gromacs/mdlib/gpuforcereduction_impl_internal.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

GpuForceReduction::Impl::Impl(const DeviceContext& deviceContext,
                              const DeviceStream&  deviceStream,
                              gmx_wallcycle*       wcycle) :
    baseForce_(),
    deviceContext_(deviceContext),
    deviceStream_(deviceStream),
    nbnxmForceToAdd_(),
    rvecForceToAdd_(),
    wcycle_(wcycle)
{
}

void GpuForceReduction::Impl::reinit(DeviceBuffer<Float3>  baseForcePtr,
                                     const int             numAtoms,
                                     ArrayRef<const int>   cell,
                                     const int             atomStart,
                                     const bool            accumulate,
                                     GpuEventSynchronizer* completionMarker)
{
    GMX_ASSERT((baseForcePtr != nullptr), "Input base force for reduction has no data");
    baseForce_        = baseForcePtr;
    numAtoms_         = numAtoms;
    atomStart_        = atomStart;
    accumulate_       = static_cast<int>(accumulate);
    completionMarker_ = completionMarker;
    cellInfo_.cell    = cell.data();

    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpu);
    reallocateDeviceBuffer(
            &cellInfo_.d_cell, numAtoms_, &cellInfo_.cellSize, &cellInfo_.cellSizeAlloc, deviceContext_);
    copyToDeviceBuffer(&cellInfo_.d_cell,
                       &(cellInfo_.cell[atomStart]),
                       0,
                       numAtoms_,
                       deviceStream_,
                       GpuApiCallBehavior::Async,
                       nullptr);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpu);

    dependencyList_.clear();
};

void GpuForceReduction::Impl::registerNbnxmForce(DeviceBuffer<RVec> forcePtr)
{
    GMX_ASSERT(forcePtr, "Input force for reduction has no data");
    nbnxmForceToAdd_ = forcePtr;
};

void GpuForceReduction::Impl::registerRvecForce(DeviceBuffer<RVec> forcePtr)
{
    GMX_ASSERT(forcePtr, "Input force for reduction has no data");
    rvecForceToAdd_ = forcePtr;
};

void GpuForceReduction::Impl::addDependency(GpuEventSynchronizer* const dependency)
{
    dependencyList_.push_back(dependency);
}

void GpuForceReduction::Impl::execute()
{
    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpu);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuNBFBufOps);

    if (numAtoms_ == 0)
    {
        return;
    }

    GMX_ASSERT(nbnxmForceToAdd_, "Nbnxm force for reduction has no data");

    // Enqueue wait on all dependencies passed
    for (auto* synchronizer : dependencyList_)
    {
        synchronizer->enqueueWaitEvent(deviceStream_);
    }

    const bool addRvecForce = static_cast<bool>(rvecForceToAdd_); // True iff initialized

    launchForceReductionKernel(numAtoms_,
                               atomStart_,
                               addRvecForce,
                               accumulate_,
                               nbnxmForceToAdd_,
                               rvecForceToAdd_,
                               baseForce_,
                               cellInfo_.d_cell,
                               deviceStream_);

    // Mark that kernel has been launched
    if (completionMarker_ != nullptr)
    {
        completionMarker_->markEvent(deviceStream_);
    }

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchGpuNBFBufOps);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpu);
}

GpuForceReduction::Impl::~Impl() = default;

GpuForceReduction::GpuForceReduction(const DeviceContext& deviceContext,
                                     const DeviceStream&  deviceStream,
                                     gmx_wallcycle*       wcycle) :
    impl_(new Impl(deviceContext, deviceStream, wcycle))
{
}

void GpuForceReduction::registerNbnxmForce(DeviceBuffer<RVec> forcePtr)
{
    impl_->registerNbnxmForce(forcePtr);
}

void GpuForceReduction::registerRvecForce(DeviceBuffer<RVec> forcePtr)
{
    impl_->registerRvecForce(forcePtr);
}

void GpuForceReduction::addDependency(GpuEventSynchronizer* const dependency)
{
    impl_->addDependency(dependency);
}

void GpuForceReduction::reinit(DeviceBuffer<RVec>    baseForcePtr,
                               const int             numAtoms,
                               ArrayRef<const int>   cell,
                               const int             atomStart,
                               const bool            accumulate,
                               GpuEventSynchronizer* completionMarker)
{
    impl_->reinit(baseForcePtr, numAtoms, cell, atomStart, accumulate, completionMarker);
}
void GpuForceReduction::execute()
{
    impl_->execute();
}

GpuForceReduction::~GpuForceReduction() = default;

} // namespace gmx
