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

#include <vector>

#include "../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * See \ref sec_namingConventions for the naming conventions.
   *
   * In the access function:
   *  - i: Output value of the function. Range: [0, m).
   *  - j: Input value of the function. First derivative direction. Range: [0,n).
   *  - k: Input value of the function. Second derivative direction. Range: [0,n).
   *
   * @tparam T_T  Data type.
   */
  /**
   * @brief General interface for Hessian access in CoDiPack.
   *
   * Helper methods which store or read data from a Hessian expect it to implement this interface.
   *
   * See \ref sec_namingConventions for the mathematical nomenclature of the arguments and components.
   *
   * @tparam T_T  The data type in the Hessian.
   */
  template<typename T_T>
  struct HessianInterface {
    public:

      using T = CODI_DD(T_T, double);  ///< See HessianInterface.

      virtual ~HessianInterface() {}  ///< Destructor

      size_t getM() const;  ///< Get size of rows (output variables).
      size_t getN() const;  ///< Get size of columns (input variables).

      /**
       * @brief Value access.
       *
       * @param[in] i  Output value of the function. Range: [0, m).
       * @param[in] j  Input value of the function. First derivative direction. Range: [0,n).
       * @param[in] k  Input value of the function. Second derivative direction. Range: [0,n).
       */
      T operator()(size_t const i, size_t const j, size_t const k) const;

      /**
       * @brief Reference access.
       *
       * \copydetails operator()(size_t const i, size_t const j, size_t const k) const
       */
      T& operator()(size_t const i, size_t const j, size_t const k);

      void resize(size_t const m, size_t const n);  ///< Resize to the new dimensions.
      size_t size() const;                          ///< Get total size of the Hessian.
  };
}
