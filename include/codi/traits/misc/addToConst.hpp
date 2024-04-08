/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2024 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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

  /**
   * @brief Turn pointers into pointers to const, and references into references to const.
   * @tparam T_Type Type to modify, if needed.
   *
   * Default implementation does not modify the type.
   */
  template<typename T_Type>
  struct AddToConstImplementation {
    public:
      using Type = T_Type;          ///< See AddToConstImplementation.
      using ModifiedType = T_Type;  ///< Modified type. Pointers/references become pointers/references to const.
  };

#ifndef DOXYGEN_DISABLE
  /// Specialization for pointer types.
  template<typename T_Type>
  struct AddToConstImplementation<T_Type*> {
    public:
      using Type = T_Type*;
      using ModifiedType = typename std::add_const<typename std::remove_const<T_Type>::type>::type*;
  };

  /// Specialization for reference types.
  template<typename T_Type>
  struct AddToConstImplementation<T_Type&> {
    public:
      using Type = T_Type&;
      using ModifiedType = typename std::add_const<typename std::remove_const<T_Type>::type>::type&;
  };
#endif

  /**
   * @brief Pointers/references become pointers/references to const.
   */
  template<typename Type>
  using AddToConst = typename AddToConstImplementation<Type>::ModifiedType;
}
