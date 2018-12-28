/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#pragma once

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
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
     * @brief The real value of the type.
     *
     * The default implementation defines the type itself as the type.
     */
    typedef T Real;

    /**
     * @brief The maximum derivative order of the type
     *
     * The default implementation assumes that the type contains only primal values.
     */
    static const size_t MaxDerivativeOrder = 0;

    /**
     * @brief Get the primal base value of the type.
     *
     * The default implementation returns the identity.
     * @param t   The value from which the base value is extracted.
     * @return The base value of the type.
     */
    static const T& getBaseValue(const T& t) { return t;}
  };
}
