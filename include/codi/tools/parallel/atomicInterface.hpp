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

#include "../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Provides a data type on which operations are performed atomically.
   *
   * If used with an underlying floating point type or an active CoDiPack type, this data type is suitable as an adjoint
   * variable type. Provides also increment and decrement operators for the use case of an underlying integer type. The
   * increment and decrement operators don't have to be implemented for underlying non-integer types.
   *
   * Implementations likely require template specializations with respect to the underlying type, especially if it is an
   * active CoDiPack type.
   *
   * @tparam T_Type  The underlying data type.
   * @tparam T_Impl  Implementing class.
   */
  template<typename T_Type, typename T_Impl>
  struct AtomicInterface {
    public:
      using Type = T_Type;  ///< See AtomicInterface.
      using Impl = T_Impl;  ///< See AtomicInterface.

      CODI_INLINE AtomicInterface() {}                        ///< Constructor
      CODI_INLINE AtomicInterface(AtomicInterface const&) {}  ///< Constructor
      CODI_INLINE AtomicInterface(Type const&) {}             ///< Constructor
      ~AtomicInterface() {}                                   ///< Destructor

      CODI_INLINE Impl& operator=(Impl const& other);  ///< Assignment operator with implementing type as rhs.
      CODI_INLINE Impl& operator=(Type const& other);  ///< Assignment operator with underlying type as rhs.

      CODI_INLINE Type operator+=(Impl const& other);  ///< Incremental update with implementing type as rhs.
      CODI_INLINE Type operator+=(Type const& other);  ///< Incremental update with underlying type as rhs.

      CODI_INLINE Type operator++();     ///< Pre-increment operator.
      CODI_INLINE Type operator++(int);  ///< Post-increment operator.
      CODI_INLINE Type operator--();     ///< Pre-decrement operator.
      CODI_INLINE Type operator--(int);  ///< Post-decrement operator.

      CODI_INLINE operator Type() const;  ///< Implicit cast to underlying type for rhs access.
  };

#if CODI_IDE
  /// Helper for IDE code completion.
  template<typename Type>
  using CODI_DEFAULT_ATOMIC = AtomicInterface<Type, CODI_IMPLEMENTATION>;
#endif
}
