/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 *  - SciComp, TU Kaiserslautern:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../../traits/realTraits.hpp"
#include "../unaryExpression.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Definitions for the ArrayAccessExpression.
  ///
  /// @tparam T_Real Return type of the expression.
  template<typename T_Real>
  struct ArrayAccessExpressionImpl {
    public:
      using Real = CODI_DD(T_Real, CODI_ANY);  ///< See ArrayAccessExpressionImpl.

      using Traits = RealTraits::AggregatedTypeTraits<Real>;  ///< Traits of the aggregated type.

      using InnerReal = typename Traits::InnerType;  ///< Inner type of the aggregate.

      /// Implementation of the array access operator for a specific element.
      /// @tparam T_element Element that is accessed.
      template<size_t T_element>
      struct ArrayAccessOperationImpl {
        public:
          /// Operation for array access.
          /// @tparam T_OpReal Real value of the operator.
          template<typename T_OpReal>
          struct type : public UnaryOperation<T_OpReal> {
            public:

              using OpReal = CODI_DD(T_OpReal, double);     ///< See type.
              static size_t constexpr element = T_element;  ///< See ArrayAccessOperationImpl.

              using Jacobian = Real;  ///< Jacobian is the aggregated type.

              /// \copydoc codi::UnaryOperation::primal()
              template<typename Arg>
              static CODI_INLINE OpReal primal(Arg const& arg) {
                return Traits::template arrayAccess<element>(arg);
              }

              /// \copydoc codi::UnaryOperation::gradient()
              template<typename Arg>
              static CODI_INLINE Jacobian gradient(Arg const& arg, OpReal const& result) {
                CODI_UNUSED(result);
                return Traits::template adjointOfArrayAccess<element>(arg, 1.0);
              }
          };
      };

      /// Definition of the array access expression.
      template<size_t element, typename Arg>
      using type = UnaryExpression<InnerReal, Arg, ArrayAccessOperationImpl<element>::template type>;
  };

  /// Expression that performs a[element] in a compile time context.
  template<typename Real, size_t element, typename Arg>
  using ArrayAccessExpression = typename ArrayAccessExpressionImpl<Real>::template type<element, Arg>;

}
