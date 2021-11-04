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
/*! \libinternal \file
 * \brief
 * Declares utilities for template metaprogramming
 *
 * \author Roland Schulz <roland.schulz@intel.com>
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_TEMPLATE_MP_H
#define GMX_UTILITY_TEMPLATE_MP_H

#include <cassert>
#include <cstddef>

#include <optional>
#include <tuple>
#include <utility>

#include "gromacs/compat/mp11.h"

namespace gmx
{

/*! \internal \brief
 * Helper function to select appropriate template based on runtime values.
 *
 * Can use enums or booleans for template parameters.
 * These enums must have a member \c Count indicating the total number of valid values.
 *
 * Example usage:
 * \code
    enum class Options {
        Op1 = 0,
        Op2 = 1,
        Count = 2
    };

    template<bool p0, Options p1, Options p2>
    bool foo(int i);

    bool bar(bool p0, Options p1, Options p2, int i) {
        return dispatchTemplatedFunction(
            [=](auto p0, auto p1, auto p2) {
                return foo<p0, p1, p2>(i);
            },
            p0, p1, p2);
    }
 * \endcode
 *
 * \tparam Function Type of \p f.
 * \param f Function to call.
 * \return The result of calling \c f().
*/
template<class Function>
auto dispatchTemplatedFunction(Function&& f)
{
    return std::forward<Function>(f)();
}

// Recursive templates confuse Doxygen
//! \cond
template<class Function, class Enum, class... Enums>
auto dispatchTemplatedFunction(Function&& f, Enum e, Enums... es)
{
    return dispatchTemplatedFunction(
            [&](auto... es_) {
                return compat::mp_with_index<size_t(Enum::Count)>(size_t(e), [&](auto e_) {
                    return std::forward<Function>(f)(
                            std::integral_constant<Enum, static_cast<Enum>(size_t(e_))>(), es_...);
                });
            },
            es...);
}

template<class Function, class... Enums>
auto dispatchTemplatedFunction(Function&& f, bool e, Enums... es)
{
    return dispatchTemplatedFunction(
            [&](auto... es_) {
                return compat::mp_with_index<2>(size_t(e), [&](auto e_) {
                    return std::forward<Function>(f)(std::bool_constant<static_cast<bool>(e_)>(), es_...);
                });
            },
            es...);
}
//! \endcond

template<int CurrentArgIndex, class Object, typename... AllArgs, typename... ResolvedArgs>
auto constructObjectWithVariadicOptionalInitializerHelper(std::tuple<std::optional<AllArgs>...> allArgs,
                                                          ResolvedArgs... resolvedArgs)
{
    if constexpr (CurrentArgIndex < 0)
    {
        return Object{ resolvedArgs... };
    }
    else
    {
        static_assert(std::tuple_size_v<decltype(allArgs)> >= CurrentArgIndex);
        auto currentArg = std::get<CurrentArgIndex>(allArgs);
        if (currentArg.has_value())
        {
            return constructObjectWithVariadicOptionalInitializerHelper<CurrentArgIndex - 1, Object>(
                    allArgs, currentArg.value(), resolvedArgs...);
        }
        else
        {
            return constructObjectWithVariadicOptionalInitializerHelper<CurrentArgIndex - 1, Object>(
                    allArgs, resolvedArgs...);
        }
    }
}

/*! \internal \brief
 * Helper function to create a class with initializer list when some of the arguments are optional.
 *
 * Example usage:
 * \code
class Foo {
    template <typename... Args>
    Foo(Args... args) {}
}

Foo bar(bool haveOne, bool haveTwo) {
    std::optional<int> one = haveOne ? 1 : std::nullopt;
    std::optional<int> two = haveTwo ? 2 : std::nullopt;
    return constructObjectWithVariadicOptionalInitializer<Foo>(one, two);
}
 * \endcode
 */
template<typename Object, typename... AllArgs>
auto constructObjectWithVariadicOptionalInitializer(std::optional<AllArgs>... allArgs)
{
    auto          allArgsTuple = std::make_tuple(allArgs...);
    constexpr int numArgs      = std::tuple_size_v<decltype(allArgsTuple)>;
    static_assert(numArgs > 0);
    return constructObjectWithVariadicOptionalInitializerHelper<numArgs - 1, Object>(std::move(allArgsTuple));
}


} // namespace gmx

#endif // GMX_UTILITY_TEMPLATE_MP_H
