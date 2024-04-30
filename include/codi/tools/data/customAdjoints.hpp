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

#include <map>

#include "../../misc/macros.hpp"
#include "../../traits/adjointVectorTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Implementation of adjoints via a map.
   *
   * Useful for the evaluation of tape parts with non-contiguous, far-apart identifiers. Can be used as custom adjoints
   * for tape evaluations, see CustomAdjointVectorEvaluationTapeInterface. Can also be used for Jacobian computations
   * with custom adjoints, see Algorithms::computeJacobianCustomAdjoints().
   *
   * The implementation ensures that any identifier can be used to access the map, both in constant and non-constant
   * calling contexts. Entries that do not exist yet are created on the fly and default-initialized.
   *
   * @tparam T_Identifier  Identifier type, usually chosen as Tape::Identifier.
   * @tparam T_Gradient    Gradient type, corresponds to the Adjoint template parameter in the referred-to functions.
   */
  template<typename T_Identifier, typename T_Gradient>
  struct MappedAdjoints {
    public:
      using Identifier = CODI_DD(T_Identifier, int);  ///< See MappedAdjoints.
      using Gradient = CODI_DD(T_Gradient, double);   ///< See MappedAdjoints.

      /// Internal map implementation, kept public for flexibility.
      /// Mutable so that access in constant call contexts does not prevent the insertion of new elements.
      mutable std::map<Identifier, Gradient> adjoints;

      /// Access operator in constant call contexts.
      Gradient const& operator[](Identifier const& i) const {
        return adjoints[i];
      }

      /// Access operator in non-constant call contexts.
      Gradient& operator[](Identifier const& i) {
        return adjoints[i];
      }
  };

#ifndef DOXYGEN_DISABLE
  // Specialize adjoint vector traits.
  namespace AdjointVectorTraits {
    template<typename T_Identifier, typename T_Gradient>
    struct GradientImplementation<MappedAdjoints<T_Identifier, T_Gradient>> {
      public:
        using Gradient = T_Gradient;
    };
  }
#endif
}
