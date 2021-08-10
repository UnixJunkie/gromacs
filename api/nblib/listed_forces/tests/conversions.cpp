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
 * \brief
 * This implements basic nblib utility tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include <numeric>

#include <gtest/gtest.h>

#include "nblib/box.h"
#include "nblib/listed_forces/conversions.hpp"
#include "nblib/listed_forces/tests/listedtesthelpers.h"

#include "testutils/testasserts.h"

namespace nblib
{
namespace test
{
namespace
{

ListedInteractionData someBondsAndAngles()
{
    ListedInteractionData         interactions;
    HarmonicBondType              bond1{ 10, 0.1 };
    HarmonicBondType              bond2{ 20, 0.2 };
    std::vector<HarmonicBondType> bonds{ bond1, bond2 };
    pickType<HarmonicBondType>(interactions).parameters = bonds;

    HarmonicAngle              angle1(100, Degrees(100));
    HarmonicAngle              angle2(200, Degrees(101));
    std::vector<HarmonicAngle> angles{ angle1, angle2 };
    pickType<HarmonicAngle>(interactions).parameters = angles;

    std::vector<InteractionIndex<HarmonicBondType>> bondIndices{ { 0, 1, 0 }, { 1, 2, 0 }, { 2, 3, 1 } };
    pickType<HarmonicBondType>(interactions).indices = std::move(bondIndices);

    std::vector<InteractionIndex<HarmonicAngle>> angleIndices{ { 0, 1, 2, 0 }, { 1, 2, 3, 1 } };
    pickType<HarmonicAngle>(interactions).indices = std::move(angleIndices);

    return interactions;
}

TEST(ListedShims, ParameterConversion)
{
    ListedInteractionData interactions = someBondsAndAngles();

    auto [idef, gmx_params] = createFFparams(interactions);

    EXPECT_EQ(gmx_params->iparams.size(), 4);
    EXPECT_EQ(gmx_params->iparams[0].harmonic.rA,
              pickType<HarmonicBondType>(interactions).parameters[0].equilConstant());
    EXPECT_REAL_EQ_TOL(gmx_params->iparams[2].harmonic.rA,
                       pickType<HarmonicAngle>(interactions).parameters[0].equilConstant() / DEG2RAD,
                       gmx::test::defaultRealTolerance());

    EXPECT_EQ(idef->il[F_BONDS].iatoms.size(), 9);
    std::vector<int> bondIatoms{ 0, 0, 1, 0, 1, 2, 1, 2, 3 };
    EXPECT_EQ(idef->il[F_BONDS].iatoms, bondIatoms);
    std::vector<int> angleIatoms{ 2, 0, 1, 2, 3, 1, 2, 3 };
    EXPECT_EQ(idef->il[F_ANGLES].iatoms, angleIatoms);
    idef->clear();
}

std::vector<gmx::RVec> twoCenterCoordinates = { { 1.382, 1.573, 1.482 }, { 1.281, 1.559, 1.596 } };

std::vector<gmx::RVec> threeCenterCoordinates = { { 1.382, 1.573, 1.482 },
                                                  { 1.281, 1.559, 1.596 },
                                                  { 1.292, 1.422, 1.663 } };

std::vector<gmx::RVec> fourCentercoordinates = { { 1.382, 1.573, 1.482 },
                                                 { 1.281, 1.559, 1.596 },
                                                 { 1.292, 1.422, 1.663 },
                                                 { 1.189, 1.407, 1.775 } };

std::array<std::vector<gmx::RVec>, 3> NCenterCoordinates{ twoCenterCoordinates,
                                                          threeCenterCoordinates,
                                                          fourCentercoordinates };

template<class Interaction>
struct TypeInput
{
    typedef Interaction         type; // needed for pickType
    ListedTypeData<Interaction> interactionData;
    // assign coordinates depending on the number of centers in the interaction type from the array above
    std::vector<gmx::RVec> coordinates = std::get<NCenter<Interaction>{} - 2>(NCenterCoordinates);
};

std::tuple TestInput{
    // Two Center Types
    TypeInput<HarmonicBondType>{
            { { { HarmonicBondType(500.0, 0.15) } }, { indexVector<HarmonicBondType>() } } },
    TypeInput<G96BondType>{ { { { G96BondType(50.0, 0.15) } }, { indexVector<G96BondType>() } } },
    TypeInput<CubicBondType>{ { { { CubicBondType(50.0, 2.0, 0.16) } }, { indexVector<CubicBondType>() } } },
    TypeInput<MorseBondType>{ { { { MorseBondType(30.0, 2.7, 0.15) } }, { indexVector<MorseBondType>() } } },
    TypeInput<FENEBondType>{ { { { FENEBondType(5.0, 0.4) } }, { indexVector<FENEBondType>() } } },
    // Three Center Types
    TypeInput<HarmonicAngle>{
            { { { HarmonicAngle(2.2, Degrees(91.0)) } }, { indexVector<HarmonicAngle>() } } },
    TypeInput<G96Angle>{ { { { G96Angle(50.0, Degrees(100)) } }, { indexVector<G96Angle>() } } },
    TypeInput<RestrictedAngle>{
            { { { RestrictedAngle(50.0, Degrees(100)) } }, { indexVector<RestrictedAngle>() } } },
    TypeInput<LinearAngle>{ { { { LinearAngle(50.0, 0.4) } }, { indexVector<LinearAngle>() } } },
    TypeInput<QuarticAngle>{ { { { QuarticAngle(1.1, 2.3, 4.6, 7.8, 9.2, Degrees(87)) } },
                               { indexVector<QuarticAngle>() } } },
    TypeInput<CrossBondBond>{ { { { CrossBondBond(45.0, 0.8, 0.7) } }, { indexVector<CrossBondBond>() } } },
    TypeInput<CrossBondAngle>{
            { { { CrossBondAngle(45.0, 0.8, 0.7, 0.3) } }, { indexVector<CrossBondAngle>() } } },
    // Four Center Types
    TypeInput<ProperDihedral>{
            { { { ProperDihedral(Degrees(45), 2.3, 1) } }, { indexVector<ProperDihedral>() } } }
};

template<class... Ts>
ListedInteractionData combineTestInput(std::tuple<Ts...> testInput)
{
    ListedInteractionData interactionData;
    // transfer all elements of testInput into the returned ListedInteractionData
    // use a lambda + for_each_tuple
    auto copyParamsOneType = [&interactionData](const auto& typeInput) {
        for (size_t i = 0; i < typeInput.interactionData.parameters.size(); i++)
        {
            auto interactionParams = typeInput.interactionData.parameters[i];
            using InteractionType  = decltype(interactionParams);
            pickType<InteractionType>(interactionData).parameters.push_back(interactionParams);

            auto indices = typeInput.interactionData.indices[i];
            pickType<InteractionType>(interactionData).indices.push_back(indices);
        }
    };
    for_each_tuple(copyParamsOneType, testInput);

    return interactionData;
}

TEST(NBlibTest, EndToEndListedComparison)
{
    int                   numParticles    = 4;
    ListedInteractionData interactionData = combineTestInput(TestInput);

    Box                    box(1.0);
    std::vector<gmx::RVec> coordinates = {
        { 1.382, 1.573, 1.482 }, { 1.281, 1.559, 1.596 }, { 1.292, 1.422, 1.663 }, { 1.189, 1.407, 1.775 }
    };

    compareNblibAndGmxListedImplementations(interactionData, coordinates, numParticles, 1, box, 1e-2);
}

template<typename Interaction>
class NblibGmxListed : public testing::Test
{
public:
    void compareNblibAndGmx()
    {
        ListedInteractionData  interactionData;
        TypeInput<Interaction> typeInput       = pickType<Interaction>(TestInput);
        pickType<Interaction>(interactionData) = typeInput.interactionData;
        compareNblibAndGmxListedImplementations(
                interactionData, typeInput.coordinates, typeInput.coordinates.size(), 1, Box(1.0), 1e-4);
    }
};

// Extract a typelist interaction types from the TestInput tuple
template<class TestInputofInteractionType>
using ExtractType = typename TestInputofInteractionType::type;
using ListedTypes = Map<ExtractType, decltype(TestInput)>;
using TestTypes   = Reduce<::testing::Types, ListedTypes>;

TYPED_TEST_CASE_P(NblibGmxListed);

TYPED_TEST_P(NblibGmxListed, SameForcesOnBoth)
{
    this->compareNblibAndGmx();
}

REGISTER_TYPED_TEST_CASE_P(NblibGmxListed, SameForcesOnBoth);

/* The following macro enclosure is needed to suppress a compiler warning from gtest
 * when instantiating a parametrized test
 * This is a known issue which is resolved in gtest v1.10 and won't be needed after update
 * https://github.com/google/googletest/pull/2316
 */
#if defined(__GNUC__) || defined(__clang__)
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wpragmas"
#    pragma GCC diagnostic ignored "-Wused-but-marked-unused"

INSTANTIATE_TYPED_TEST_CASE_P(CompareEachTypeInNblibAndGmx, NblibGmxListed, TestTypes);

#    pragma GCC diagnostic pop

#elif
INSTANTIATE_TYPED_TEST_CASE_P(CompareEachTypeInNblibAndGmx, NblibGmxListed, TestTypes);
#endif
} // namespace
} // namespace test
} // namespace nblib
