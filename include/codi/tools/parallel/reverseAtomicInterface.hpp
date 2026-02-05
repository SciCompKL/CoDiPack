/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), RPTU University Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, RPTU University Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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
 *  - SciComp, RPTU University Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Provides a data type on which += update operations are performed atomically.
   *
   * In a multithreaded environment, data races on adjoint variables are fixed by performing updates atomically, whereas
   * other read and write operations do not need to be atomic. This template takes an ordinary adjoint variable type,
   * like a floating point type or a CoDiPack forward type, and ensures that the corresponding += update operation is
   * performed atomically.
   *
   * Implementations likely require template specializations with respect to the underlying type, especially if it is an
   * active CoDiPack type.
   *
   * An implementation should preserve the memory footprint of the underlying type, e.g., by inheriting from the under-
   * lying type or by having a variable of the underlying type as the only member variable.
   *
   * @tparam T_Type  The underlying data type.
   * @tparam T_Impl  Implementing class.
   */
  template<typename T_Type, typename T_Impl>
  struct ReverseAtomicInterface {
    public:
      using Type = T_Type;  ///< See ReverseAtomicInterface.
      using Impl = T_Impl;  ///< See ReverseAtomicInterface.

      CODI_INLINE ReverseAtomicInterface() {}                               ///< Constructor
      CODI_INLINE ReverseAtomicInterface(ReverseAtomicInterface const&) {}  ///< Constructor
      CODI_INLINE ReverseAtomicInterface(Type const&) {}                    ///< Constructor
      ~ReverseAtomicInterface() {}                                          ///< Destructor

      /// Assignment operator with implementing type as rhs. Not atomic.
      CODI_INLINE Impl& operator=(Impl const& other);
      CODI_INLINE Impl& operator=(Type const& other);  ///< Assignment operator with underlying type as rhs. Not atomic.

      CODI_INLINE void operator+=(Impl const& other);  ///< Atomic incremental update with implementing type as rhs.
      CODI_INLINE void operator+=(Type const& other);  ///< Atomic incremental update with underlying type as rhs.

      CODI_INLINE operator Type() const;  ///< Implicit cast to underlying type for rhs access. Not atomic.
  };

#if CODI_IDE
  /// Helper for IDE code completion.
  template<typename Type>
  using CODI_DEFAULT_REVERSE_ATOMIC = ReverseAtomicInterface<Type, CODI_IMPLEMENTATION>;
#endif
}
