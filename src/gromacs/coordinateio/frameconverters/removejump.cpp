/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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
/*!\file
 * \internal
 * \brief
 * Implements gmx::RemoveJump.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */

#include "gmxpre.h"

#include "removejump.h"

#include "gromacs/trajectory/trajectoryframe.h"

namespace gmx
{

void RemoveJump::convertFrame(t_trxframe* input)
{
    int natoms = input->natoms;
    GMX_RELEASE_ASSERT(
            natoms == referenceCoord_.ssize(),
            "Need to have same size of reference coordinates and number of atoms in frame");
    RVec diagonal(0, 0, 0);
    for (int d = 0; d < DIM; d++)
    {
        diagonal[d] = localBox_[d][d];
    }
    for (int i = 0; i < natoms; i++)
    {
        for (int m = DIM - 1; m >= 0; m--)
        {
            if (diagonal[m] > 0)
            {
                while (input->x[i][m] - referenceCoord_[i][m] <= -diagonal[m])
                {
                    for (int d = 0; d <= m; d++)
                    {
                        input->x[i][d] += localBox_[m][d];
                    }
                }
                while (input->x[i][m] - referenceCoord_[i][m] > diagonal[m])
                {
                    for (int d = 0; d <= m; d++)
                    {
                        input->x[i][d] -= localBox_[m][d];
                    }
                }
            }
        }
    }
}

} // namespace gmx
