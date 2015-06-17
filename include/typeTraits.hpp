/**
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing, TU Kaiserslautern
 *
 * This file is part of CoDiPack.
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
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
 * Authors: Max Sagebaum, Tim Albring, TU Kaiserslautern.
 */

#pragma once

namespace codi {

  /**
   * @brief Provides information about the types which are used in the active types.
   *
   * This is the general implementation for all types which are used in the active type
   * template parameter. It is used in CoDiPack to gather information about the the specific
   * type.
   *
   * @tparam T The type used in the active types.
   */
  template<typename T>
  class TypeTraits {
  public:
    /**
     * @brief The passive value of the type.
     *
     * The default implementation defines the type itself as the passive type.
     */
    typedef T PassiveReal;

    /**
     * @brief Get the primal base value of the type.
     *
     * The default implementation returns the identity.
     * @param t   The value from which the base value is extracted.
     * @return The base value of the type.
     */
    static const T getBaseValue(const T& t) { return t;}
  };
}
