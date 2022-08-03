/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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

#ifndef GMX_GMXANA_CLUSTER_LINKAGE_H
#define GMX_GMXANA_CLUSTER_LINKAGE_H

#include <stdio.h>

#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "icluster.h"

struct t_clusters;
struct t_mat;

namespace gmx
{

class MDLogger;

class ClusterLinkage : public ICluster
{
public:
    explicit ClusterLinkage(const t_mat* inputMatrix, real rmsdCutOff, int numFrames, const MDLogger& logger) :
        finished_(false), rmsdCutOff_(rmsdCutOff), matrix_(inputMatrix), logger_(logger)
    {
        clusters_.resize(numFrames);
        makeClusters();
    }
    ~ClusterLinkage() override = default;

    ArrayRef<const int> clusterList() const override;

private:
    //! Perform actual clustering.
    void makeClusters();
    //! Did we peform the clustering?
    bool finished_;
    //! Value for RMSD cutoff.
    const real rmsdCutOff_;
    //! Handle to cluster matrix.
    const t_mat* matrix_;
    //! Cluster indices
    std::vector<int> clusters_;
    //! Logger handle
    const MDLogger& logger_;
};

} // namespace gmx

#endif
