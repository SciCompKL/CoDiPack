/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../../traits/realTraits.hpp"
#include "../computeExpression.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Expression that performs a[element] in a compile time context.
  ///
  /// Based on the array access operators defined in AggregatedTypeTraits.
  ///
  /// @tparam T_Real Aggregated type of a. (E.g. std::complex<double>)
  /// @tparam T_element Element that is accessed.
  template<typename T_Aggregated, size_t T_element>
  struct ArrayAccessExpressionImpl {
    public:
      using Aggregated = CODI_DD(T_Aggregated, CODI_ANY);  ///< See ArrayAccessExpressionImpl.
      static size_t constexpr element = T_element;         ///< See ArrayAccessExpressionImpl.

      using Traits = RealTraits::AggregatedTypeTraits<Aggregated>;  ///< Traits of the aggregated type.

      using InnerReal = typename Traits::InnerType;  ///< Inner type of the aggregate.

      /// Operation for array access.
      /// @tparam T_OpReal Real value of the operator.
      template<typename T_Real>
      struct ArrayAccessOperation : public UnaryOperation<T_Real, ArrayAccessOperation<T_Real>> {
        public:

          using Real = CODI_DD(T_Real, double);  ///< See ArrayAccessOperation.

          using Jacobian = Aggregated;  ///< Jacobian is the aggregated type.

          /// \copydoc codi::UnaryOperation::primal()
          template<typename Arg>
          static CODI_INLINE Real primal(Arg const& arg) {
            return Traits::template arrayAccess<element>(arg);
          }

          /// \copydoc UnaryOperation::applyTangentArg
          template<typename Tangent, typename Arg>
          static CODI_INLINE auto applyTangentArg(Tangent const& tangent, Real const& result, Arg const& arg);
          // TODO: Implement

          /// \copydoc UnaryOperation::applyAdjointArg
          template<typename Adjoint, typename Arg>
          static CODI_INLINE auto applyAdjointArg(Adjoint const& adjoint, Real const& result, Arg const& arg) {
            CODI_UNUSED(arg);

            return Traits::template adjointOfArrayAccess<element>(result, adjoint);
          }
      };

      /// Definition of the array access expression.
      template<typename Arg>
      using Expression = ComputeExpression<InnerReal, ArrayAccessOperation, Arg>;
  };

  /// Expression that performs a[element] in a compile time context.
  template<typename Aggregated, size_t element, typename Arg>
  using ArrayAccessExpression = typename ArrayAccessExpressionImpl<Aggregated, element>::template Expression<Arg>;

}
