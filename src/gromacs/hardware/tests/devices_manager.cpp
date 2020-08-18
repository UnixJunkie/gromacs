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
 * Tests for DevicesManager
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_hardware
 */
#include "gmxpre.h"

#include "gromacs/hardware/devices_manager.h"

#include "config.h"

#include <algorithm>

#include <gtest/gtest.h>

#include "gromacs/hardware/device_information.h"
#include "gromacs/utility/inmemoryserializer.h"
#include "gromacs/utility/stringutil.h"

namespace
{

TEST(DevicesManagerTest, Serialization)
{
    std::vector<std::unique_ptr<DeviceInformation>> deviceInfosIn = DevicesManager::findDevices();
    gmx::InMemorySerializer                         writer;
    DevicesManager::serializeDeviceInformations(deviceInfosIn, &writer);
    auto buffer = writer.finishAndGetBuffer();

    gmx::InMemoryDeserializer                       reader(buffer, false);
    std::vector<std::unique_ptr<DeviceInformation>> deviceInfosOut =
            DevicesManager::deserializeDeviceInformations(&reader);

    EXPECT_EQ(deviceInfosOut.size(), deviceInfosIn.size())
            << "Number of accessible devices changed after serialization/deserialization.";

    for (int deviceId = 0; deviceId < static_cast<int>(deviceInfosIn.size()); deviceId++)
    {
        EXPECT_FALSE(deviceInfosIn[deviceId] == nullptr) << gmx::formatString(
                "Device #%d information is nullptr before serialization.", deviceId);
        EXPECT_FALSE(deviceInfosOut[deviceId] == nullptr) << gmx::formatString(
                "Device #%d information is nullptr after serialization.", deviceId);

        const DeviceInformation& deviceInfoIn  = *deviceInfosIn[deviceId];
        const DeviceInformation& deviceInfoOut = *deviceInfosOut[deviceId];
        EXPECT_EQ(deviceInfoIn.status, deviceInfoOut.status) << gmx::formatString(
                "Device status changed after serialization/deserialization for device #%d.", deviceId);

        EXPECT_EQ(deviceInfoIn.id, deviceInfoOut.id) << gmx::formatString(
                "Device id changed after serialization/deserialization for device #%d.", deviceId);

#if GMX_GPU_OPENCL
        EXPECT_EQ(deviceInfoIn.oclPlatformId, deviceInfoOut.oclPlatformId) << gmx::formatString(
                "Device OpenCL platform ID changed after serialization/deserialization for device "
                "#%d.",
                deviceId);

#endif // GMX_GPU_OPENCL
    }
}

} // namespace
