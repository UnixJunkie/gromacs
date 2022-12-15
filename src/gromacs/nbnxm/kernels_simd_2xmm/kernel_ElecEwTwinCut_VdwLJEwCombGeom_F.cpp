/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
/*
 * Note: this file was generated by the Verlet kernel generator for
 * kernel type 2xmm.
 */

/* Some target architectures compile kernels for only some NBNxN
 * kernel flavours, but the code is generated before the target
 * architecture is known. So compilation is conditional upon
 * KernelLayout::r2xMM, so that this file reduces to a stub
 * function definition when the kernel will never be called.
 */
#include "gmxpre.h"

#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/nbnxm/nbnxm_simd.h"
#if GMX_HAVE_NBNXM_SIMD_2XMM
#    include "gromacs/nbnxm/simd_kernel.h"

namespace gmx
{

template void
nbnxmKernelSimd<KernelLayout::r2xMM, KernelCoulombType::EwaldAnalytical, true, LJCombinationRule::None, InteractionModifiers::PotShift, true, EnergyOutput::None>(
        const NbnxnPairlistCpu*    nbl,
        const nbnxn_atomdata_t*    nbat,
        const interaction_const_t* ic,
        const rvec*                shift_vec,
        nbnxn_atomdata_output_t*   out);

} // namespace gmx

#endif
