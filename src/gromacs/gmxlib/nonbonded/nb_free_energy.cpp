/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020,2021, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "nb_free_energy.h"

#include "config.h"

#include <cmath>

#include <algorithm>

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/arrayref.h"


//! Scalar (non-SIMD) data types.
struct ScalarDataTypes
{
    using RealType                     = real; //!< The data type to use as real.
    using IntType                      = int;  //!< The data type to use as int.
    static constexpr int simdRealWidth = 1;    //!< The width of the RealType.
    static constexpr int simdIntWidth  = 1;    //!< The width of the IntType.
};

#if GMX_SIMD_HAVE_REAL && GMX_SIMD_HAVE_INT32_ARITHMETICS
//! SIMD data types.
struct SimdDataTypes
{
    using RealType                     = gmx::SimdReal;         //!< The data type to use as real.
    using IntType                      = gmx::SimdInt32;        //!< The data type to use as int.
    static constexpr int simdRealWidth = GMX_SIMD_REAL_WIDTH;   //!< The width of the RealType.
    static constexpr int simdIntWidth  = GMX_SIMD_FINT32_WIDTH; //!< The width of the IntType.
};
#endif

//! Computes r^(1/p) and 1/r^(1/p) for the standard p=6
template<class RealType>
static inline void pthRoot(const RealType r, RealType* pthRoot, RealType* invPthRoot)
{
    *invPthRoot = gmx::invsqrt(std::cbrt(r));
    *pthRoot    = 1 / (*invPthRoot);
}

template<class RealType>
static inline RealType calculateRinv6(const RealType rInvV)
{
    RealType rInv6 = rInvV * rInvV;
    return (rInv6 * rInv6 * rInv6);
}

template<class RealType>
static inline RealType calculateVdw6(const RealType c6, const RealType rInv6)
{
    return (c6 * rInv6);
}

template<class RealType>
static inline RealType calculateVdw12(const RealType c12, const RealType rInv6)
{
    return (c12 * rInv6 * rInv6);
}

/* reaction-field electrostatics */
template<class RealType>
static inline RealType reactionFieldScalarForce(const RealType qq,
                                                const RealType rInv,
                                                const RealType r,
                                                const real     krf,
                                                const real     two)
{
    return (qq * (rInv - two * krf * r * r));
}
template<class RealType>
static inline RealType reactionFieldPotential(const RealType qq,
                                              const RealType rInv,
                                              const RealType r,
                                              const real     krf,
                                              const real     potentialShift)
{
    return (qq * (rInv + krf * r * r - potentialShift));
}

/* Ewald electrostatics */
template<class RealType>
static inline RealType ewaldScalarForce(const RealType coulomb, const RealType rInv)
{
    return (coulomb * rInv);
}
template<class RealType>
static inline RealType ewaldPotential(const RealType coulomb, const RealType rInv, const real potentialShift)
{
    return (coulomb * (rInv - potentialShift));
}

/* cutoff LJ */
template<class RealType>
static inline RealType lennardJonesScalarForce(const RealType v6, const RealType v12)
{
    return (v12 - v6);
}
template<class RealType>
static inline RealType lennardJonesPotential(const RealType v6,
                                             const RealType v12,
                                             const RealType c6,
                                             const RealType c12,
                                             const real     repulsionShift,
                                             const real     dispersionShift,
                                             const real     oneSixth,
                                             const real     oneTwelfth)
{
    return ((v12 + c12 * repulsionShift) * oneTwelfth - (v6 + c6 * dispersionShift) * oneSixth);
}

/* Ewald LJ */
static inline real ewaldLennardJonesGridSubtract(const real c6grid, const real potentialShift, const real oneSixth)
{
    return (c6grid * potentialShift * oneSixth);
}

/* LJ Potential switch */
template<class RealType>
static inline RealType potSwitchScalarForceMod(const RealType fScalarInp,
                                               const RealType potential,
                                               const RealType sw,
                                               const RealType r,
                                               const RealType rVdw,
                                               const RealType dsw,
                                               const real     zero)
{
    if (r < rVdw)
    {
        real fScalar = fScalarInp * sw - r * potential * dsw;
        return (fScalar);
    }
    return (zero);
}
template<class RealType>
static inline RealType potSwitchPotentialMod(const RealType potentialInp,
                                             const RealType sw,
                                             const RealType r,
                                             const RealType rVdw,
                                             const real     zero)
{
    if (r < rVdw)
    {
        real potential = potentialInp * sw;
        return (potential);
    }
    return (zero);
}


//! Templated free-energy non-bonded kernel
template<typename DataTypes, bool useSoftCore, bool scLambdasOrAlphasDiffer, bool vdwInteractionTypeIsEwald, bool elecInteractionTypeIsEwald, bool vdwModifierIsPotSwitch>
static void nb_free_energy_kernel(const t_nblist&                nlist,
                                  gmx::ArrayRef<const gmx::RVec> coords,
                                  gmx::ForceWithShiftForces*     forceWithShiftForces,
                                  const t_forcerec&              fr,
                                  gmx::ArrayRef<const real>      chargeA,
                                  gmx::ArrayRef<const real>      chargeB,
                                  gmx::ArrayRef<const int>       typeA,
                                  gmx::ArrayRef<const int>       typeB,
                                  int                            flags,
                                  gmx::ArrayRef<const real>      lambda,
                                  gmx::ArrayRef<real>            dvdl,
                                  gmx::ArrayRef<real>            energygrp_elec,
                                  gmx::ArrayRef<real>            energygrp_vdw,
                                  t_nrnb* gmx_restrict nrnb)
{
    constexpr int STATE_A = 0;
    constexpr int STATE_B = 1;
    constexpr int NSTATES = 2;

    using RealType = typename DataTypes::RealType;
    using IntType  = typename DataTypes::IntType;

    /* FIXME: How should these be handled with SIMD? */
    constexpr real oneTwelfth = 1.0 / 12.0;
    constexpr real oneSixth   = 1.0 / 6.0;
    constexpr real zero       = 0.0;
    constexpr real half       = 0.5;
    constexpr real one        = 1.0;
    constexpr real two        = 2.0;
    constexpr real six        = 6.0;

    /* Extract pointer to non-bonded interaction constants */
    const interaction_const_t* ic = fr.ic.get();

    const auto& scParams      = *ic->softCoreParameters;
    const bool  doForces      = ((flags & GMX_NONBONDED_DO_FORCE) != 0);
    const bool  doShiftForces = ((flags & GMX_NONBONDED_DO_SHIFTFORCE) != 0);
    const bool  doPotential   = ((flags & GMX_NONBONDED_DO_POTENTIAL) != 0);

    // Note that the nbnxm kernels do not support Coulomb potential switching at all
    GMX_ASSERT(ic->coulomb_modifier != InteractionModifiers::PotSwitch,
               "Potential switching is not supported for Coulomb with FEP");

    real vdw_swV3, vdw_swV4, vdw_swV5, vdw_swF2, vdw_swF3, vdw_swF4;
    if (vdwModifierIsPotSwitch)
    {
        const real d = ic->rvdw - ic->rvdw_switch;
        vdw_swV3     = -10.0 / (d * d * d);
        vdw_swV4     = 15.0 / (d * d * d * d);
        vdw_swV5     = -6.0 / (d * d * d * d * d);
        vdw_swF2     = -30.0 / (d * d * d);
        vdw_swF3     = 60.0 / (d * d * d * d);
        vdw_swF4     = -30.0 / (d * d * d * d * d);
    }
    else
    {
        /* Avoid warnings from stupid compilers (looking at you, Clang!) */
        vdw_swV3 = vdw_swV4 = vdw_swV5 = vdw_swF2 = vdw_swF3 = vdw_swF4 = 0.0;
    }

    NbkernelElecType icoul;
    if (ic->eeltype == CoulombInteractionType::Cut || EEL_RF(ic->eeltype))
    {
        icoul = NbkernelElecType::ReactionField;
    }
    else
    {
        icoul = NbkernelElecType::None;
    }

    real rcutoff_max2 = std::max(ic->rcoulomb, ic->rvdw);
    rcutoff_max2      = rcutoff_max2 * rcutoff_max2;

    const real* tab_ewald_F_lj           = nullptr;
    const real* tab_ewald_V_lj           = nullptr;
    const real* ewtab                    = nullptr;
    real        coulombTableScale        = 0;
    real        coulombTableScaleInvHalf = 0;
    real        vdwTableScale            = 0;
    real        vdwTableScaleInvHalf     = 0;
    real        sh_ewald                 = 0;
    if (elecInteractionTypeIsEwald || vdwInteractionTypeIsEwald)
    {
        sh_ewald = ic->sh_ewald;
    }
    if (elecInteractionTypeIsEwald)
    {
        const auto& coulombTables = *ic->coulombEwaldTables;
        ewtab                     = coulombTables.tableFDV0.data();
        coulombTableScale         = coulombTables.scale;
        coulombTableScaleInvHalf  = half / coulombTableScale;
    }
    if (vdwInteractionTypeIsEwald)
    {
        const auto& vdwTables = *ic->vdwEwaldTables;
        tab_ewald_F_lj        = vdwTables.tableF.data();
        tab_ewald_V_lj        = vdwTables.tableV.data();
        vdwTableScale         = vdwTables.scale;
        vdwTableScaleInvHalf  = half / vdwTableScale;
    }

    /* For Ewald/PME interactions we cannot easily apply the soft-core component to
     * reciprocal space. When we use non-switched Ewald interactions, we
     * assume the soft-coring does not significantly affect the grid contribution
     * and apply the soft-core only to the full 1/r (- shift) pair contribution.
     *
     * However, we cannot use this approach for switch-modified since we would then
     * effectively end up evaluating a significantly different interaction here compared to the
     * normal (non-free-energy) kernels, either by applying a cutoff at a different
     * position than what the user requested, or by switching different
     * things (1/r rather than short-range Ewald). For these settings, we just
     * use the traditional short-range Ewald interaction in that case.
     */
    GMX_RELEASE_ASSERT(!(vdwInteractionTypeIsEwald && vdwModifierIsPotSwitch),
                       "Can not apply soft-core to switched Ewald potentials");

    real dvdlCoul = 0;
    real dvdlVdw  = 0;

    /* Lambda factor for state A, 1-lambda*/
    const real lambda_coul = lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)];
    const real lambda_vdw  = lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)];
    real       LFC[NSTATES], LFV[NSTATES];
    LFC[STATE_A] = one - lambda_coul;
    LFV[STATE_A] = one - lambda_vdw;

    /* Lambda factor for state B, lambda*/
    LFC[STATE_B] = lambda_coul;
    LFV[STATE_B] = lambda_vdw;

    /*derivative of the lambda factor for state A and B */
    real DLF[NSTATES];
    DLF[STATE_A] = -1;
    DLF[STATE_B] = 1;

    real           lFacCoul[NSTATES], dlFacCoul[NSTATES], lFacVdw[NSTATES], dlFacVdw[NSTATES];
    constexpr real sc_r_power = 6.0_real;
    for (int i = 0; i < NSTATES; i++)
    {
        lFacCoul[i]  = (scParams.lambdaPower == 2 ? (1 - LFC[i]) * (1 - LFC[i]) : (1 - LFC[i]));
        dlFacCoul[i] = DLF[i] * scParams.lambdaPower / sc_r_power
                       * (scParams.lambdaPower == 2 ? (1 - LFC[i]) : 1);
        lFacVdw[i]  = (scParams.lambdaPower == 2 ? (1 - LFV[i]) * (1 - LFV[i]) : (1 - LFV[i]));
        dlFacVdw[i] = DLF[i] * scParams.lambdaPower / sc_r_power
                      * (scParams.lambdaPower == 2 ? (1 - LFV[i]) : 1);
    }

    // TODO: We should get rid of using pointers to real
    const real* x             = coords[0];
    real* gmx_restrict f      = &(forceWithShiftForces->force()[0][0]);
    real* gmx_restrict fshift = &(forceWithShiftForces->shiftForces()[0][0]);

    const real rlistSquared = gmx::square(fr.rlist);

    int numExcludedPairsBeyondRlist = 0;

    for (int n = 0; n < nlist.nri; n++)
    {
        int npair_within_cutoff = 0;

        const int  is    = nlist.shift[n];
        const int  is3   = DIM * is;
        const real shX   = fr.shift_vec[is][XX];
        const real shY   = fr.shift_vec[is][YY];
        const real shZ   = fr.shift_vec[is][ZZ];
        const int  nj0   = nlist.jindex[n];
        const int  nj1   = nlist.jindex[n + 1];
        const int  ii    = nlist.iinr[n];
        const int  ii3   = 3 * ii;
        const real ix    = shX + x[ii3 + 0];
        const real iy    = shY + x[ii3 + 1];
        const real iz    = shZ + x[ii3 + 2];
        const real iqA   = ic->epsfac * chargeA[ii];
        const real iqB   = ic->epsfac * chargeB[ii];
        const int  ntiA  = 2 * fr.ntype * typeA[ii];
        const int  ntiB  = 2 * fr.ntype * typeB[ii];
        real       vCTot = 0;
        real       vVTot = 0;
        real       fIX   = 0;
        real       fIY   = 0;
        real       fIZ   = 0;

        for (int k = nj0; k < nj1; k++)
        {
            int            tj[NSTATES];
            const int      jnr = nlist.jjnr[k];
            const int      j3  = 3 * jnr;
            RealType       c6[NSTATES], c12[NSTATES], qq[NSTATES], vCoul[NSTATES], vVdw[NSTATES];
            RealType       r, rInv, rp, rpm2;
            RealType       alphaVdwEff, alphaCoulEff, sigma6[NSTATES];
            const RealType dX  = ix - x[j3];
            const RealType dY  = iy - x[j3 + 1];
            const RealType dZ  = iz - x[j3 + 2];
            const RealType rSq = dX * dX + dY * dY + dZ * dZ;
            RealType       fScalC[NSTATES], fScalV[NSTATES];
            /* Check if this pair on the exlusions list.*/
            const bool bPairIncluded = nlist.excl_fep.empty() || nlist.excl_fep[k];

            if (rSq >= rcutoff_max2 && bPairIncluded)
            {
                /* We save significant time by skipping all code below.
                 * Note that with soft-core interactions, the actual cut-off
                 * check might be different. But since the soft-core distance
                 * is always larger than r, checking on r here is safe.
                 * Exclusions outside the cutoff can not be skipped as
                 * when using Ewald: the reciprocal-space
                 * Ewald component still needs to be subtracted.
                 */

                continue;
            }
            npair_within_cutoff++;

            if (rSq > rlistSquared)
            {
                numExcludedPairsBeyondRlist++;
            }

            if (rSq > 0)
            {
                /* Note that unlike in the nbnxn kernels, we do not need
                 * to clamp the value of rSq before taking the invsqrt
                 * to avoid NaN in the LJ calculation, since here we do
                 * not calculate LJ interactions when C6 and C12 are zero.
                 */

                rInv = gmx::invsqrt(rSq);
                r    = rSq * rInv;
            }
            else
            {
                /* The force at r=0 is zero, because of symmetry.
                 * But note that the potential is in general non-zero,
                 * since the soft-cored r will be non-zero.
                 */
                rInv = 0;
                r    = 0;
            }

            if (useSoftCore)
            {
                rpm2 = rSq * rSq;  /* r4 */
                rp   = rpm2 * rSq; /* r6 */
            }
            else
            {
                /* The soft-core power p will not affect the results
                 * with not using soft-core, so we use power of 0 which gives
                 * the simplest math and cheapest code.
                 */
                rpm2 = rInv * rInv;
                rp   = 1;
            }

            RealType fScal = 0;

            qq[STATE_A] = iqA * chargeA[jnr];
            qq[STATE_B] = iqB * chargeB[jnr];

            tj[STATE_A] = ntiA + 2 * typeA[jnr];
            tj[STATE_B] = ntiB + 2 * typeB[jnr];

            if (bPairIncluded)
            {
                c6[STATE_A] = fr.nbfp[tj[STATE_A]];
                c6[STATE_B] = fr.nbfp[tj[STATE_B]];

                for (int i = 0; i < NSTATES; i++)
                {
                    c12[i] = fr.nbfp[tj[i] + 1];
                    if (useSoftCore)
                    {
                        if ((c6[i] > 0) && (c12[i] > 0))
                        {
                            /* c12 is stored scaled with 12.0 and c6 is scaled with 6.0 - correct for this */
                            sigma6[i] = half * c12[i] / c6[i];
                            if (sigma6[i] < scParams.sigma6Minimum) /* for disappearing coul and vdw with soft core at the same time */
                            {
                                sigma6[i] = scParams.sigma6Minimum;
                            }
                        }
                        else
                        {
                            sigma6[i] = scParams.sigma6WithInvalidSigma;
                        }
                    }
                }

                if (useSoftCore)
                {
                    /* only use softcore if one of the states has a zero endstate - softcore is for avoiding infinities!*/
                    if ((c12[STATE_A] > 0) && (c12[STATE_B] > 0))
                    {
                        alphaVdwEff  = 0;
                        alphaCoulEff = 0;
                    }
                    else
                    {
                        alphaVdwEff  = scParams.alphaVdw;
                        alphaCoulEff = scParams.alphaCoulomb;
                    }
                }

                for (int i = 0; i < NSTATES; i++)
                {
                    fScalC[i] = 0;
                    fScalV[i] = 0;
                    vCoul[i]  = 0;
                    vVdw[i]   = 0;

                    RealType rInvC, rInvV, rC, rV, rPInvC, rPInvV;

                    /* Only spend time on A or B state if it is non-zero */
                    if ((qq[i] != 0) || (c6[i] != 0) || (c12[i] != 0))
                    {
                        /* this section has to be inside the loop because of the dependence on sigma6 */
                        if (useSoftCore)
                        {
                            rPInvC = one / (alphaCoulEff * lFacCoul[i] * sigma6[i] + rp);
                            pthRoot(rPInvC, &rInvC, &rC);
                            if (scLambdasOrAlphasDiffer)
                            {
                                rPInvV = one / (alphaVdwEff * lFacVdw[i] * sigma6[i] + rp);
                                pthRoot(rPInvV, &rInvV, &rV);
                            }
                            else
                            {
                                /* We can avoid one expensive pow and one / operation */
                                rPInvV = rPInvC;
                                rInvV  = rInvC;
                                rV     = rC;
                            }
                        }
                        else
                        {
                            rPInvC = 1;
                            rInvC  = rInv;
                            rC     = r;

                            rPInvV = 1;
                            rInvV  = rInv;
                            rV     = r;
                        }

                        /* Only process the coulomb interactions if we have charges,
                         * and if we either include all entries in the list (no cutoff
                         * used in the kernel), or if we are within the cutoff.
                         */
                        bool computeElecInteraction =
                                (elecInteractionTypeIsEwald && r < ic->rcoulomb)
                                || (!elecInteractionTypeIsEwald && rC < ic->rcoulomb);

                        if ((qq[i] != 0) && computeElecInteraction)
                        {
                            if (elecInteractionTypeIsEwald)
                            {
                                vCoul[i]  = ewaldPotential(qq[i], rInvC, sh_ewald);
                                fScalC[i] = ewaldScalarForce(qq[i], rInvC);
                            }
                            else
                            {
                                vCoul[i] = reactionFieldPotential(
                                        qq[i], rInvC, rC, ic->reactionFieldCoefficient, ic->reactionFieldShift);
                                fScalC[i] = reactionFieldScalarForce(
                                        qq[i], rInvC, rC, ic->reactionFieldCoefficient, two);
                            }
                        }

                        /* Only process the VDW interactions if we have
                         * some non-zero parameters, and if we either
                         * include all entries in the list (no cutoff used
                         * in the kernel), or if we are within the cutoff.
                         */
                        bool computeVdwInteraction = (vdwInteractionTypeIsEwald && r < ic->rvdw)
                                                     || (!vdwInteractionTypeIsEwald && rV < ic->rvdw);
                        if ((c6[i] != 0 || c12[i] != 0) && computeVdwInteraction)
                        {
                            RealType rInv6;
                            if (useSoftCore)
                            {
                                rInv6 = rPInvV;
                            }
                            else
                            {
                                rInv6 = calculateRinv6(rInvV);
                            }
                            RealType vVdw6  = calculateVdw6(c6[i], rInv6);
                            RealType vVdw12 = calculateVdw12(c12[i], rInv6);

                            vVdw[i]   = lennardJonesPotential(vVdw6,
                                                            vVdw12,
                                                            c6[i],
                                                            c12[i],
                                                            ic->repulsion_shift.cpot,
                                                            ic->dispersion_shift.cpot,
                                                            oneSixth,
                                                            oneTwelfth);
                            fScalV[i] = lennardJonesScalarForce(vVdw6, vVdw12);

                            if (vdwInteractionTypeIsEwald)
                            {
                                /* Subtract the grid potential at the cut-off */
                                vVdw[i] += ewaldLennardJonesGridSubtract(
                                        fr.ljpme_c6grid[tj[i]], ic->sh_lj_ewald, oneSixth);
                            }

                            if (vdwModifierIsPotSwitch)
                            {
                                RealType d        = rV - ic->rvdw_switch;
                                d                 = (d > zero) ? d : zero;
                                const RealType d2 = d * d;
                                const RealType sw =
                                        one + d2 * d * (vdw_swV3 + d * (vdw_swV4 + d * vdw_swV5));
                                const RealType dsw = d2 * (vdw_swF2 + d * (vdw_swF3 + d * vdw_swF4));

                                fScalV[i] = potSwitchScalarForceMod(
                                        fScalV[i], vVdw[i], sw, rV, ic->rvdw, dsw, zero);
                                vVdw[i] = potSwitchPotentialMod(vVdw[i], sw, rV, ic->rvdw, zero);
                            }
                        }

                        /* fScalC (and fScalV) now contain: dV/drC * rC
                         * Now we multiply by rC^-p, so it will be: dV/drC * rC^1-p
                         * Further down we first multiply by r^p-2 and then by
                         * the vector r, which in total gives: dV/drC * (r/rC)^1-p
                         */
                        fScalC[i] *= rPInvC;
                        fScalV[i] *= rPInvV;
                    }
                } // end for (int i = 0; i < NSTATES; i++)

                /* Assemble A and B states */
                for (int i = 0; i < NSTATES; i++)
                {
                    vCTot += LFC[i] * vCoul[i];
                    vVTot += LFV[i] * vVdw[i];

                    fScal += LFC[i] * fScalC[i] * rpm2;
                    fScal += LFV[i] * fScalV[i] * rpm2;

                    if (useSoftCore)
                    {
                        dvdlCoul += vCoul[i] * DLF[i]
                                    + LFC[i] * alphaCoulEff * dlFacCoul[i] * fScalC[i] * sigma6[i];
                        dvdlVdw += vVdw[i] * DLF[i]
                                   + LFV[i] * alphaVdwEff * dlFacVdw[i] * fScalV[i] * sigma6[i];
                    }
                    else
                    {
                        dvdlCoul += vCoul[i] * DLF[i];
                        dvdlVdw += vVdw[i] * DLF[i];
                    }
                }
            } // end if (bPairIncluded)
            else if (icoul == NbkernelElecType::ReactionField)
            {
                /* For excluded pairs, which are only in this pair list when
                 * using the Verlet scheme, we don't use soft-core.
                 * As there is no singularity, there is no need for soft-core.
                 */
                const real FF = -two * ic->reactionFieldCoefficient;
                RealType   VV = ic->reactionFieldCoefficient * rSq - ic->reactionFieldShift;

                if (ii == jnr)
                {
                    VV *= half;
                }

                for (int i = 0; i < NSTATES; i++)
                {
                    vCTot += LFC[i] * qq[i] * VV;
                    fScal += LFC[i] * qq[i] * FF;
                    dvdlCoul += DLF[i] * qq[i] * VV;
                }
            }

            if (elecInteractionTypeIsEwald && (r < ic->rcoulomb || !bPairIncluded))
            {
                /* See comment in the preamble. When using Ewald interactions
                 * (unless we use a switch modifier) we subtract the reciprocal-space
                 * Ewald component here which made it possible to apply the free
                 * energy interaction to 1/r (vanilla coulomb short-range part)
                 * above. This gets us closer to the ideal case of applying
                 * the softcore to the entire electrostatic interaction,
                 * including the reciprocal-space component.
                 */
                real v_lr, f_lr;

                const RealType ewrt   = r * coulombTableScale;
                IntType        ewitab = static_cast<IntType>(ewrt);
                const RealType eweps  = ewrt - ewitab;
                ewitab                = 4 * ewitab;
                f_lr                  = ewtab[ewitab] + eweps * ewtab[ewitab + 1];
                v_lr = (ewtab[ewitab + 2] - coulombTableScaleInvHalf * eweps * (ewtab[ewitab] + f_lr));
                f_lr *= rInv;

                /* Note that any possible Ewald shift has already been applied in
                 * the normal interaction part above.
                 */

                if (ii == jnr)
                {
                    /* If we get here, the i particle (ii) has itself (jnr)
                     * in its neighborlist. This can only happen with the Verlet
                     * scheme, and corresponds to a self-interaction that will
                     * occur twice. Scale it down by 50% to only include it once.
                     */
                    v_lr *= half;
                }

                for (int i = 0; i < NSTATES; i++)
                {
                    vCTot -= LFC[i] * qq[i] * v_lr;
                    fScal -= LFC[i] * qq[i] * f_lr;
                    dvdlCoul -= (DLF[i] * qq[i]) * v_lr;
                }
            }

            if (vdwInteractionTypeIsEwald && (r < ic->rvdw || !bPairIncluded))
            {
                /* See comment in the preamble. When using LJ-Ewald interactions
                 * (unless we use a switch modifier) we subtract the reciprocal-space
                 * Ewald component here which made it possible to apply the free
                 * energy interaction to r^-6 (vanilla LJ6 short-range part)
                 * above. This gets us closer to the ideal case of applying
                 * the softcore to the entire VdW interaction,
                 * including the reciprocal-space component.
                 */
                /* We could also use the analytical form here
                 * iso a table, but that can cause issues for
                 * r close to 0 for non-interacting pairs.
                 */

                const RealType rs   = rSq * rInv * vdwTableScale;
                const IntType  ri   = static_cast<IntType>(rs);
                const RealType frac = rs - ri;
                const RealType f_lr = (1 - frac) * tab_ewald_F_lj[ri] + frac * tab_ewald_F_lj[ri + 1];
                /* TODO: Currently the Ewald LJ table does not contain
                 * the factor 1/6, we should add this.
                 */
                const RealType FF = f_lr * rInv / six;
                RealType       VV =
                        (tab_ewald_V_lj[ri] - vdwTableScaleInvHalf * frac * (tab_ewald_F_lj[ri] + f_lr))
                        / six;

                if (ii == jnr)
                {
                    /* If we get here, the i particle (ii) has itself (jnr)
                     * in its neighborlist. This can only happen with the Verlet
                     * scheme, and corresponds to a self-interaction that will
                     * occur twice. Scale it down by 50% to only include it once.
                     */
                    VV *= half;
                }

                for (int i = 0; i < NSTATES; i++)
                {
                    const real c6grid = fr.ljpme_c6grid[tj[i]];
                    vVTot += LFV[i] * c6grid * VV;
                    fScal += LFV[i] * c6grid * FF;
                    dvdlVdw += (DLF[i] * c6grid) * VV;
                }
            }

            if (doForces)
            {
                const real tX = fScal * dX;
                const real tY = fScal * dY;
                const real tZ = fScal * dZ;
                fIX           = fIX + tX;
                fIY           = fIY + tY;
                fIZ           = fIZ + tZ;
                /* OpenMP atomics are expensive, but this kernels is also
                 * expensive, so we can take this hit, instead of using
                 * thread-local output buffers and extra reduction.
                 *
                 * All the OpenMP regions in this file are trivial and should
                 * not throw, so no need for try/catch.
                 */
#pragma omp atomic
                f[j3] -= tX;
#pragma omp atomic
                f[j3 + 1] -= tY;
#pragma omp atomic
                f[j3 + 2] -= tZ;
            }
        } // end for (int k = nj0; k < nj1; k++)

        /* The atomics below are expensive with many OpenMP threads.
         * Here unperturbed i-particles will usually only have a few
         * (perturbed) j-particles in the list. Thus with a buffered list
         * we can skip a significant number of i-reductions with a check.
         */
        if (npair_within_cutoff > 0)
        {
            if (doForces)
            {
#pragma omp atomic
                f[ii3] += fIX;
#pragma omp atomic
                f[ii3 + 1] += fIY;
#pragma omp atomic
                f[ii3 + 2] += fIZ;
            }
            if (doShiftForces)
            {
#pragma omp atomic
                fshift[is3] += fIX;
#pragma omp atomic
                fshift[is3 + 1] += fIY;
#pragma omp atomic
                fshift[is3 + 2] += fIZ;
            }
            if (doPotential)
            {
                int ggid = nlist.gid[n];
#pragma omp atomic
                energygrp_elec[ggid] += vCTot;
#pragma omp atomic
                energygrp_vdw[ggid] += vVTot;
            }
        }
    } // end for (int n = 0; n < nlist.nri; n++)

#pragma omp atomic
    dvdl[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)] += dvdlCoul;
#pragma omp atomic
    dvdl[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)] += dvdlVdw;

    /* Estimate flops, average for free energy stuff:
     * 12  flops per outer iteration
     * 150 flops per inner iteration
     */
    atomicNrnbIncrement(nrnb, eNR_NBKERNEL_FREE_ENERGY, nlist.nri * 12 + nlist.jindex[nlist.nri] * 150);

    if (numExcludedPairsBeyondRlist > 0)
    {
        gmx_fatal(FARGS,
                  "There are %d perturbed non-bonded pair interactions beyond the pair-list cutoff "
                  "of %g nm, which is not supported. This can happen because the system is "
                  "unstable or because intra-molecular interactions at long distances are "
                  "excluded. If the "
                  "latter is the case, you can try to increase nstlist or rlist to avoid this."
                  "The error is likely triggered by the use of couple-intramol=no "
                  "and the maximal distance in the decoupled molecule exceeding rlist.",
                  numExcludedPairsBeyondRlist,
                  fr.rlist);
    }
}

typedef void (*KernelFunction)(const t_nblist&                nlist,
                               gmx::ArrayRef<const gmx::RVec> coords,
                               gmx::ForceWithShiftForces*     forceWithShiftForces,
                               const t_forcerec&              fr,
                               gmx::ArrayRef<const real>      chargeA,
                               gmx::ArrayRef<const real>      chargeB,
                               gmx::ArrayRef<const int>       typeA,
                               gmx::ArrayRef<const int>       typeB,
                               int                            flags,
                               gmx::ArrayRef<const real>      lambda,
                               gmx::ArrayRef<real>            dvdl,
                               gmx::ArrayRef<real>            energygrp_elec,
                               gmx::ArrayRef<real>            energygrp_vdw,
                               t_nrnb* gmx_restrict nrnb);

template<bool useSoftCore, bool scLambdasOrAlphasDiffer, bool vdwInteractionTypeIsEwald, bool elecInteractionTypeIsEwald, bool vdwModifierIsPotSwitch>
static KernelFunction dispatchKernelOnUseSimd(const bool useSimd)
{
    if (useSimd)
    {
#if GMX_SIMD_HAVE_REAL && GMX_SIMD_HAVE_INT32_ARITHMETICS && GMX_USE_SIMD_KERNELS
        /* FIXME: Here SimdDataTypes should be used to enable SIMD. So far, the code in nb_free_energy_kernel is not adapted to SIMD */
        return (nb_free_energy_kernel<ScalarDataTypes, useSoftCore, scLambdasOrAlphasDiffer, vdwInteractionTypeIsEwald, elecInteractionTypeIsEwald, vdwModifierIsPotSwitch>);
#else
        return (nb_free_energy_kernel<ScalarDataTypes, useSoftCore, scLambdasOrAlphasDiffer, vdwInteractionTypeIsEwald, elecInteractionTypeIsEwald, vdwModifierIsPotSwitch>);
#endif
    }
    else
    {
        return (nb_free_energy_kernel<ScalarDataTypes, useSoftCore, scLambdasOrAlphasDiffer, vdwInteractionTypeIsEwald, elecInteractionTypeIsEwald, vdwModifierIsPotSwitch>);
    }
}

template<bool useSoftCore, bool scLambdasOrAlphasDiffer, bool vdwInteractionTypeIsEwald, bool elecInteractionTypeIsEwald>
static KernelFunction dispatchKernelOnVdwModifier(const bool vdwModifierIsPotSwitch, const bool useSimd)
{
    if (vdwModifierIsPotSwitch)
    {
        return (dispatchKernelOnUseSimd<useSoftCore, scLambdasOrAlphasDiffer, vdwInteractionTypeIsEwald, elecInteractionTypeIsEwald, true>(
                useSimd));
    }
    else
    {
        return (dispatchKernelOnUseSimd<useSoftCore, scLambdasOrAlphasDiffer, vdwInteractionTypeIsEwald, elecInteractionTypeIsEwald, false>(
                useSimd));
    }
}

template<bool useSoftCore, bool scLambdasOrAlphasDiffer, bool vdwInteractionTypeIsEwald>
static KernelFunction dispatchKernelOnElecInteractionType(const bool elecInteractionTypeIsEwald,
                                                          const bool vdwModifierIsPotSwitch,
                                                          const bool useSimd)
{
    if (elecInteractionTypeIsEwald)
    {
        return (dispatchKernelOnVdwModifier<useSoftCore, scLambdasOrAlphasDiffer, vdwInteractionTypeIsEwald, true>(
                vdwModifierIsPotSwitch, useSimd));
    }
    else
    {
        return (dispatchKernelOnVdwModifier<useSoftCore, scLambdasOrAlphasDiffer, vdwInteractionTypeIsEwald, false>(
                vdwModifierIsPotSwitch, useSimd));
    }
}

template<bool useSoftCore, bool scLambdasOrAlphasDiffer>
static KernelFunction dispatchKernelOnVdwInteractionType(const bool vdwInteractionTypeIsEwald,
                                                         const bool elecInteractionTypeIsEwald,
                                                         const bool vdwModifierIsPotSwitch,
                                                         const bool useSimd)
{
    if (vdwInteractionTypeIsEwald)
    {
        return (dispatchKernelOnElecInteractionType<useSoftCore, scLambdasOrAlphasDiffer, true>(
                elecInteractionTypeIsEwald, vdwModifierIsPotSwitch, useSimd));
    }
    else
    {
        return (dispatchKernelOnElecInteractionType<useSoftCore, scLambdasOrAlphasDiffer, false>(
                elecInteractionTypeIsEwald, vdwModifierIsPotSwitch, useSimd));
    }
}

template<bool useSoftCore>
static KernelFunction dispatchKernelOnScLambdasOrAlphasDifference(const bool scLambdasOrAlphasDiffer,
                                                                  const bool vdwInteractionTypeIsEwald,
                                                                  const bool elecInteractionTypeIsEwald,
                                                                  const bool vdwModifierIsPotSwitch,
                                                                  const bool useSimd)
{
    if (scLambdasOrAlphasDiffer)
    {
        return (dispatchKernelOnVdwInteractionType<useSoftCore, true>(
                vdwInteractionTypeIsEwald, elecInteractionTypeIsEwald, vdwModifierIsPotSwitch, useSimd));
    }
    else
    {
        return (dispatchKernelOnVdwInteractionType<useSoftCore, false>(
                vdwInteractionTypeIsEwald, elecInteractionTypeIsEwald, vdwModifierIsPotSwitch, useSimd));
    }
}

static KernelFunction dispatchKernel(const bool                 scLambdasOrAlphasDiffer,
                                     const bool                 vdwInteractionTypeIsEwald,
                                     const bool                 elecInteractionTypeIsEwald,
                                     const bool                 vdwModifierIsPotSwitch,
                                     const bool                 useSimd,
                                     const interaction_const_t& ic)
{
    if (ic.softCoreParameters->alphaCoulomb == 0 && ic.softCoreParameters->alphaVdw == 0)
    {
        return (dispatchKernelOnScLambdasOrAlphasDifference<false>(scLambdasOrAlphasDiffer,
                                                                   vdwInteractionTypeIsEwald,
                                                                   elecInteractionTypeIsEwald,
                                                                   vdwModifierIsPotSwitch,
                                                                   useSimd));
    }
    else
    {
        return (dispatchKernelOnScLambdasOrAlphasDifference<true>(scLambdasOrAlphasDiffer,
                                                                  vdwInteractionTypeIsEwald,
                                                                  elecInteractionTypeIsEwald,
                                                                  vdwModifierIsPotSwitch,
                                                                  useSimd));
    }
}


void gmx_nb_free_energy_kernel(const t_nblist&                nlist,
                               gmx::ArrayRef<const gmx::RVec> coords,
                               gmx::ForceWithShiftForces*     ff,
                               const t_forcerec&              fr,
                               gmx::ArrayRef<const real>      chargeA,
                               gmx::ArrayRef<const real>      chargeB,
                               gmx::ArrayRef<const int>       typeA,
                               gmx::ArrayRef<const int>       typeB,
                               int                            flags,
                               gmx::ArrayRef<const real>      lambda,
                               gmx::ArrayRef<real>            dvdl,
                               gmx::ArrayRef<real>            energygrp_elec,
                               gmx::ArrayRef<real>            energygrp_vdw,
                               t_nrnb*                        nrnb)
{
    const interaction_const_t& ic = *fr.ic;
    GMX_ASSERT(EEL_PME_EWALD(ic.eeltype) || ic.eeltype == CoulombInteractionType::Cut || EEL_RF(ic.eeltype),
               "Unsupported eeltype with free energy");
    GMX_ASSERT(ic.softCoreParameters, "We need soft-core parameters");

    const auto& scParams                   = *ic.softCoreParameters;
    const bool  vdwInteractionTypeIsEwald  = (EVDW_PME(ic.vdwtype));
    const bool  elecInteractionTypeIsEwald = (EEL_PME_EWALD(ic.eeltype));
    const bool  vdwModifierIsPotSwitch     = (ic.vdw_modifier == InteractionModifiers::PotSwitch);
    bool        scLambdasOrAlphasDiffer    = true;
    const bool  useSimd                    = fr.use_simd_kernels;

    if (scParams.alphaCoulomb == 0 && scParams.alphaVdw == 0)
    {
        scLambdasOrAlphasDiffer = false;
    }
    else
    {
        if (lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)]
                    == lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)]
            && scParams.alphaCoulomb == scParams.alphaVdw)
        {
            scLambdasOrAlphasDiffer = false;
        }
    }

    KernelFunction kernelFunc;
    kernelFunc = dispatchKernel(scLambdasOrAlphasDiffer,
                                vdwInteractionTypeIsEwald,
                                elecInteractionTypeIsEwald,
                                vdwModifierIsPotSwitch,
                                useSimd,
                                ic);
    kernelFunc(nlist, coords, ff, fr, chargeA, chargeB, typeA, typeB, flags, lambda, dvdl, energygrp_elec, energygrp_vdw, nrnb);
}
