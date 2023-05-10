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

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../../traits/misc/enableIfHelpers.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief General interface for Jacobian access in CoDiPack.
   *
   * Helper methods which store or read data from a Jacobian expect it to implement this interface.
   *
   * See \ref sec_namingConventions for the mathematical nomenclature of the arguments and components.
   *
   * @tparam T_T  The data type in the Jacobian.
   */
  template<typename T_T>
  struct JacobianInterface {
    public:

      using T = CODI_DD(T_T, double);  ///< See JacobianInterface.

      size_t getM() const;  ///< Get size of rows (output variables).
      size_t getN() const;  ///< Get size of columns (input variables).

      /// Value access, i in [0, ..., m), j in [0, ..., n).
      T operator()(size_t const i, size_t const j) const;

      /// Reference access, i in [0, ..., m), j in [0, ..., n).
      T& operator()(size_t const i, size_t const j);

      void resize(size_t const m, size_t const n);  ///< Resize the Jacobian.
      size_t size() const;                          ///< Get total size of the Jacobian.

      /// Interface for the JacobianDelayAccessor
      void setLogic(size_t const i, size_t const j, T const& v);
  };

  /**
   * @brief Output a Jacobian on the data stream.
   *
   * The format is (Matlab):
   * [1, 2, 3;
   *  4, 5, 6;
   *  7, 8, 8]
   */
  template<typename Stream, typename Jac, typename = enable_if_base_of<Jac, JacobianInterface<typename Jac::T>>>
  Stream& operator<<(Stream& out, CODI_DD(Jac, CODI_T(JacobianInterface<double>)) const& jacobian) {
    out << "[";
    for (size_t i = 0; i < jacobian.getM(); ++i) {
      if (i != 0) {
        out << " ";  // Padding for the '['.
      }

      for (size_t j = 0; j < jacobian.getN(); ++j) {
        if (j != 0) {
          out << ", ";
        }
        out << jacobian(i, j);
      }

      if (i + 1 < jacobian.getM()) {
        out << ";\n";
      } else {
        out << "]";
      }
    }

    return out;
  }
}
