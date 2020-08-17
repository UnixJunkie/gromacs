/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 * \brief
 * Tests to compare that multiple time stepping is (nearly) identical to normal integration.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/mpitest.h"
#include "testutils/setenv.h"
#include "testutils/simulationdatabase.h"

#include "moduletest.h"
#include "simulatorcomparison.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief Test fixture base for two integration schemes
 *
 * This test ensures that integration with(out) different multiple time stepping
 * scheems (called via different mdp options) yield near identical energies,
 * forces and virial at step 0 and similar energies and virial after 4 steps.
 */
using MtsComparisonTestParams = std::tuple<std::string, std::string, int>;
class MtsComparisonTest : public MdrunTestFixture, public ::testing::WithParamInterface<MtsComparisonTestParams>
{
};

TEST_P(MtsComparisonTest, WithinTolerances)
{
    auto params         = GetParam();
    auto simulationName = std::get<0>(params);
    auto mtsScheme      = std::get<1>(params);
    auto numSteps       = std::get<2>(params);

    // Note that there should be no relevant limitation on MPI ranks and OpenMP threads
    SCOPED_TRACE(formatString("Comparing for '%s' no MTS with MTS scheme '%s'",
                              simulationName.c_str(), mtsScheme.c_str()));

    auto sharedMdpOptions = gmx::formatString(
            "integrator   = md\n"
            "dt           = 0.001\n"
            "nsteps       = %d\n"
            "verlet-buffer-tolerance = -1\n"
            "rlist        = 1.0\n"
            "coulomb-type = PME\n"
            "vdw-type     = cut-off\n"
            "rcoulomb     = 0.9\n"
            "rvdw         = 0.9\n"
            "constraints  = h-bonds\n",
            numSteps);

    const int nstfout       = (numSteps == 0 ? 4 : 0);
    auto      refMdpOptions = sharedMdpOptions
                         + gmx::formatString(
                                   "mts       = no\n"
                                   "nstcalcenergy = 4\n"
                                   "nstenergy = 4\n"
                                   "nstxout   = 0\n"
                                   "nstvout   = 0\n"
                                   "nstfout   = %d\n",
                                   nstfout);

    auto mtsMdpOptions = sharedMdpOptions
                         + gmx::formatString(
                                   "mts        = yes\n"
                                   "mts-scheme = %s\n"
                                   "nstcalcenergy = 4\n"
                                   "nstenergy  = 4\n"
                                   "nstxout    = 0\n"
                                   "nstvout    = 0\n"
                                   "nstfout    = %d\n",
                                   mtsScheme.c_str(), nstfout);

    // With numSteps=0 the energy and virial should only differ due to rounding errors
    const real           energyTol = (numSteps == 0 ? 0.001 : (mtsScheme == "pme" ? 0.015 : 0.04));
    const real           virialTol = (numSteps == 0 ? 0.01 : (mtsScheme == "pme" ? 0.1 : 0.2));
    EnergyTermsToCompare energyTermsToCompare{
        { { interaction_function[F_EPOT].longname, relativeToleranceAsFloatingPoint(100.0, energyTol) },
          { "Vir-XX", relativeToleranceAsFloatingPoint(30.0, virialTol) },
          { "Vir-YY", relativeToleranceAsFloatingPoint(30.0, virialTol) },
          { "Vir-ZZ", relativeToleranceAsFloatingPoint(30.0, virialTol) } }
    };

    // Specify how trajectory frame matching must work.
    TrajectoryFrameMatchSettings trajectoryMatchSettings{ true,
                                                          true,
                                                          true,
                                                          ComparisonConditions::NoComparison,
                                                          ComparisonConditions::NoComparison,
                                                          ComparisonConditions::MustCompare };
    TrajectoryTolerances trajectoryTolerances = TrajectoryComparison::s_defaultTrajectoryTolerances;

    // Build the functor that will compare reference and test
    // trajectory frames in the chosen way.
    TrajectoryComparison trajectoryComparison{ trajectoryMatchSettings, trajectoryTolerances };

    // Set file names
    auto simulator1TrajectoryFileName = fileManager_.getTemporaryFilePath("sim1.trr");
    auto simulator1EdrFileName        = fileManager_.getTemporaryFilePath("sim1.edr");
    auto simulator2TrajectoryFileName = fileManager_.getTemporaryFilePath("sim2.trr");
    auto simulator2EdrFileName        = fileManager_.getTemporaryFilePath("sim2.edr");

    // Run grompp
    runner_.tprFileName_ = fileManager_.getTemporaryFilePath("sim.tpr");
    runner_.useTopGroAndNdxFromDatabase(simulationName);
    runner_.useStringAsMdpFile(refMdpOptions);
    runGrompp(&runner_);

    // Do first mdrun
    runner_.fullPrecisionTrajectoryFileName_ = simulator1TrajectoryFileName;
    runner_.edrFileName_                     = simulator1EdrFileName;
    runMdrun(&runner_);

    runner_.useStringAsMdpFile(mtsMdpOptions);
    runGrompp(&runner_);

    // Do second mdrun
    runner_.fullPrecisionTrajectoryFileName_ = simulator2TrajectoryFileName;
    runner_.edrFileName_                     = simulator2EdrFileName;
    runMdrun(&runner_);

    // Compare simulation results
    compareEnergies(simulator1EdrFileName, simulator2EdrFileName, energyTermsToCompare);
    if (numSteps == 0)
    {
        compareTrajectories(simulator1TrajectoryFileName, simulator2TrajectoryFileName, trajectoryComparison);
    }
}

INSTANTIATE_TEST_CASE_P(MultipleTimeSteppingIsNearSingleTimeStepping,
                        MtsComparisonTest,
                        ::testing::Combine(::testing::Values("ala"),
                                           ::testing::Values("pme", "pme-nonbonded-pair-dihedral"),
                                           ::testing::Values(0, 4)));

} // namespace
} // namespace test
} // namespace gmx
