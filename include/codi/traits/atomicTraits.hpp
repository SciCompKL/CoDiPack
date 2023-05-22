/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * CoDiPack is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU
 * General Public License along with CoDiPack.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <type_traits>

/** \copydoc codi::Namespace */
namespace codi {

  /// Traits for atomic types.
  namespace AtomicTraits {

    /**
     * @brief Indicate whether a type is atomic.
     *
     * @tparam Type  The type to be checked.
     */
    template<typename Type>
    struct IsAtomic : std::false_type {};

    /// Enable if abbreviation for IsAtomic.
    template<typename Type>
    using EnableIfAtomic = typename std::enable_if<IsAtomic<Type>::value>::type;

    /**
     * @brief Convert an atomic type into the underlying type. Implementation.
     *
     * Default for non-atomic values.
     *
     * @tparam T_NotAtomic  Non-atomic type.
     */
    template<typename T_NotAtomic, typename = void>
    struct RemoveAtomicImpl {
      public:
        using Type = T_NotAtomic;  ///< See RemoveAtomicImpl.
    };

    /**
     * @brief Specialization for atomic types.
     *
     * @tparam Atomic  Atomic type.
     */
    template<typename Atomic>
    struct RemoveAtomicImpl<Atomic, EnableIfAtomic<Atomic>> {
      public:
        using Type = typename Atomic::Type;  ///< See RemoveAtomicImpl.
    };

    /// Wrapper for removing atomic from a type.
    template<typename Type>
    using RemoveAtomic = typename RemoveAtomicImpl<Type>::Type;
  }
}
