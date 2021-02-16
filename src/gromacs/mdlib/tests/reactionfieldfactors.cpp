/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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
 * \brief Reaction field parameters tests
 *
 * \author Joe Jordan <ejjordan@kth.se>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "config.h"


#include <gtest/gtest.h>

#include "gromacs/mdlib/reactionfieldfactors.h"
#include "gromacs/mdtypes/md_enums.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

class reactionFieldFactorsTest : public ::testing::Test
{
};

TEST_F(reactionFieldFactorsTest, ReactionFieldFactorPOTSHIFT)
{
    ReactionFieldCoefficients RFCoeffs(1, 1, 1, false, eintmodPOTSHIFT);
    EXPECT_EQ(RFCoeffs.constant_, 0);
    EXPECT_EQ(RFCoeffs.correction_, 1);
}

TEST_F(reactionFieldFactorsTest, ReactionFieldFactorNONE)
{
    ReactionFieldCoefficients RFCoeffs(1, 1, 1, false, eintmodNONE);
    EXPECT_EQ(RFCoeffs.constant_, 0);
    EXPECT_EQ(RFCoeffs.correction_, 0);
}

} // namespace
} // namespace test
} // namespace gmx
