/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2022, by the GROMACS development team, led by
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
 * Tests for the MSM class.
 *
 * \author Cathrine Bergh
 */

#include "gmxpre.h"
#include "gromacs/msm/msm.h"

#include <gtest/gtest.h>

namespace gmx
{
namespace test
{

// class MsmTest : public ::testing::Test
// {
// protected:
//     void SetUp() override
//     {
//       //msm = MarkovModel(4);
//     }
//     //MarkovModel msm;
// };

// // TODO: better reuse data in the tests

// TEST_F(MsmTest, TransitionCountingTest)
// {
//     // TODO: Make test pass or fail
//     //std::vector<int> discretizedTraj = {0, 0, 0, 0, 0, 3, 3, 2};
//     std::vector<int> discretizedTraj = {0, 1, 3, 2, 3, 3, 3, 2};

//     MarkovModel countTestMsm = MarkovModel(4);

//     countTestMsm.countTransitions(discretizedTraj, 1);
//     auto& transitions = countTestMsm.transitionCountsMatrix;

//     const auto& dataView = transitions.asConstView();
//     const int numRows = transitions.extent(0);
//     const int numCols = transitions.extent(1);

//     for (int i = 0; i < numRows; i++)
//     {
//         printf("\n");
//         for (int j=0; j < numCols; j++)
//         {
//             printf("%d ", dataView[i][j]);
//         }
//     }
//     printf("\n");
// }

// TEST_F(MsmTest, TransitionProbabilityTest)
// {
//     //std::vector<int> discretizedTraj = {0, 0, 0, 0, 0, 3, 3, 2};
//     std::vector<int> discretizedTraj = {0, 1, 3, 2, 3, 3, 3, 2};
//     MarkovModel tpmTestMsm = MarkovModel(4);
//     tpmTestMsm.countTransitions(discretizedTraj, 1);

//     tpmTestMsm.computeTransitionProbabilities();
//     auto& probs = tpmTestMsm.transitionProbabilityMatrix;

//     const auto& dataView = probs.asConstView();
//     const int numRows = probs.extent(0);
//     const int numCols = probs.extent(1);

//     for (int i = 0; i < numRows; i++)
//     {
//         printf("\n");
//         for (int j=0; j < numCols; j++)
//         {
//             printf("%f ", dataView[i][j]);
//         }
//     }
//     printf("\n");
// }


TEST(Diagtest, Diagonalize)
{
    /*
    // Should give eigenvalues (3, -1) OK
    MatrixNxM matrix(2,2);
    matrix(0,0)=1;
    matrix(0,1)=2;
    matrix(1,0)=2;
    matrix(1,1)=1;
    */

    // Should give eigenvalues (3, 1) OK
    MatrixNxM matrix(2,2);
    matrix(0,0)=2;
    matrix(0,1)=1;
    matrix(1,0)=1;
    matrix(1,1)=2;

    /*
    // Should give eigenvalues (2+sqrt(2), 2-sqrt(2)) NOT OK
    MatrixNxM matrix(2,2);
    matrix(0,0)=2;
    matrix(0,1)=2;
    matrix(1,0)=1;
    matrix(1,1)=2;
    */

    /*
    // Should give eigenvalues (2+sqrt(5), 2-sqrt(5)) NOT OK
    MatrixNxM matrix(2,2);
    matrix(0,0)=2;
    matrix(0,1)=1;
    matrix(1,0)=5;
    matrix(1,1)=2;
    */

    /*
    // Should give eigenvalues (5, -1) NOT OK
    MatrixNxM matrix(2,2);
    matrix(0,0)=2;
    matrix(0,1)=1;
    matrix(1,0)=9;
    matrix(1,1)=2;
    */

    /*
    // Should get eigenvalues (0, 2) OK
    MatrixNxM matrix(2,2);
    matrix(0,0)=1;
    matrix(0,1)=1;
    matrix(1,0)=1;
    matrix(1,1)=1;
    */

    /*
    // Should give eigenvalues (0, 0, 0, 4) OK
    MatrixNxM matrix(4,4);
    matrix(0,0)=1;
    matrix(0,1)=1;
    matrix(0,2)=1;
    matrix(0,3)=1;
    matrix(1,0)=1;
    matrix(1,1)=1;
    matrix(1,2)=1;
    matrix(1,3)=1;
    matrix(2,0)=1;
    matrix(2,1)=1;
    matrix(2,2)=1;
    matrix(2,3)=1;
    matrix(3,0)=1;
    matrix(3,1)=1;
    matrix(3,2)=1;
    matrix(3,3)=1;
    */

    /*
    // Should give eigenvalues (0, 0, 0.8, 0.5) NOT OK
    MatrixNxM matrix(4,4);
    matrix(0,0)=0.8;
    matrix(0,1)=0.0;
    matrix(0,2)=0.0;
    matrix(0,3)=0.2;
    matrix(1,0)=0.0;
    matrix(1,1)=0.0;
    matrix(1,2)=0.0;
    matrix(1,3)=0.0;
    matrix(2,0)=0.0;
    matrix(2,1)=0.0;
    matrix(2,2)=0.0;
    matrix(2,3)=0.0;
    matrix(3,0)=0.0;
    matrix(3,1)=0.0;
    matrix(3,2)=0.5;
    matrix(3,3)=0.5;
    */

    /*
    // Should give eigenvalues (0, 0, 2, 4) NOT OK
    MatrixNxM matrix(4,4);
    matrix(0,0)=4;
    matrix(0,1)=0;
    matrix(0,2)=0;
    matrix(0,3)=2;
    matrix(1,0)=0;
    matrix(1,1)=0;
    matrix(1,2)=0;
    matrix(1,3)=0;
    matrix(2,0)=0;
    matrix(2,1)=0;
    matrix(2,2)=0;
    matrix(2,3)=0;
    matrix(3,0)=0;
    matrix(3,1)=0;
    matrix(3,2)=3;
    matrix(3,3)=2;
    */

    const auto result = diagonalize(matrix);

    fprintf(stderr,"\n eigenvalue 1 = %f ", result.eigenvalues_[0]);
    fprintf(stderr,"\n eigenvalue 2 = %f ", result.eigenvalues_[1]);
    fprintf(stderr,"\n eigenvalue 3 = %f ", result.eigenvalues_[2]);
    fprintf(stderr,"\n eigenvalue 4 = %f ", result.eigenvalues_[3]);

}

// TEST_F(MsmTest, DiagonalizationTest)
// {
//     std::vector<int> discretizedTraj = {0, 0, 0, 0, 0, 3, 3, 2}; //runs, should get eigenvalues (0.8, 0.5, 0, 0)
//     //std::vector<int> discretizedTraj = {2, 2, 0, 0, 0, 3, 3, 2}; //runs
//     //std::vector<int> discretizedTraj = {0, 1, 1, 2, 3, 3, 3, 2}; //floating point exception
//     MarkovModel diagTestMsm = MarkovModel(4);
//     diagTestMsm.countTransitions(discretizedTraj, 1);
//     diagTestMsm.computeTransitionProbabilities();

//     auto tpm = diagTestMsm.transitionProbabilityMatrix;

//     // TODO: if diagonalization is not called before
//     // attributes are collected, they will be zero.
//     // Better to initialize them by running the method?
//     diagTestMsm.diagonalizeMatrix(tpm);

//     auto eigvals = diagTestMsm.eigenvalues;
//     auto eigvecs = diagTestMsm.eigenvectors;

//     for (int i = 0; i < eigvals.size(); ++i)
//     {
//         printf("Val elm %d: %f\n", i, eigvals[i]);
//     }

//     for (int i = 0; i < eigvecs.size(); ++i)
//     {
//         printf("Vec elm %d: %f\n", i, eigvecs[i]);
//     }

// }

} //namespace test
} //namespace gmx
