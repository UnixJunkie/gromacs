/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2018,2019,2021, by the GROMACS development team, led by
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

#include "gromacs/utility/enumerationhelpers.h"
#include "timecontrol.h"

#include <mutex>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"

/* The source code in this file should be thread-safe.
         Please keep it that way. */

/* Globals for trajectory input */
struct t_timecontrol
{
    t_timecontrol(real inputT, bool inputSet) : t(inputT), bSet(inputSet) {}
    real t;
    bool bSet;
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static gmx::EnumerationArray<TimeControl, t_timecontrol> timecontrol = { t_timecontrol(0, false),
                                                                         t_timecontrol(0, false),
                                                                         t_timecontrol(0, false) };

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static std::mutex g_timeControlMutex;

gmx_bool bTimeSet(TimeControl tcontrol)
{
    gmx_bool ret;

    const std::lock_guard<std::mutex> lock(g_timeControlMutex);
    ret = timecontrol[tcontrol].bSet;

    return ret;
}

real rTimeValue(TimeControl tcontrol)
{
    real ret;

    const std::lock_guard<std::mutex> lock(g_timeControlMutex);
    ret = timecontrol[tcontrol].t;
    return ret;
}

void setTimeValue(TimeControl tcontrol, real value)
{
    const std::lock_guard<std::mutex> lock(g_timeControlMutex);
    timecontrol[tcontrol].t    = value;
    timecontrol[tcontrol].bSet = TRUE;
}
