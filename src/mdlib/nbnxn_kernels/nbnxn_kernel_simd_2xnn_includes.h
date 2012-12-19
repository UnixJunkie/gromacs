/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
 * Copyright (c) 2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

/* This files includes all x86 SIMD kernel flavors.
 * Only the Electrostatics type and optionally the VdW cut-off check
 * need to be set before including this file.
 */

/* Include the force+energy kernels */
#define CALC_ENERGIES
#define LJ_COMB_GEOM
#include "nbnxn_kernel_simd_2xnn_outer.h"
#undef LJ_COMB_GEOM
#define LJ_COMB_LB
#include "nbnxn_kernel_simd_2xnn_outer.h"
#undef LJ_COMB_LB
#include "nbnxn_kernel_simd_2xnn_outer.h"
#undef CALC_ENERGIES

/* Include the force+energygroups kernels */
#define CALC_ENERGIES
#define ENERGY_GROUPS
#define LJ_COMB_GEOM
#include "nbnxn_kernel_simd_2xnn_outer.h"
#undef LJ_COMB_GEOM
#define LJ_COMB_LB
#include "nbnxn_kernel_simd_2xnn_outer.h"
#undef LJ_COMB_LB
#include "nbnxn_kernel_simd_2xnn_outer.h"
#undef ENERGY_GROUPS
#undef CALC_ENERGIES

/* Include the force only kernels */
#define LJ_COMB_GEOM
#include "nbnxn_kernel_simd_2xnn_outer.h"
#undef LJ_COMB_GEOM
#define LJ_COMB_LB
#include "nbnxn_kernel_simd_2xnn_outer.h"
#undef LJ_COMB_LB
#include "nbnxn_kernel_simd_2xnn_outer.h"
