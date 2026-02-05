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

#include "../../config.h"
#include "../../misc/compileTimeLoop.hpp"
#include "../../misc/eventSystem.hpp"
#include "../../misc/macros.hpp"
#include "../../traits/realTraits.hpp"
#include "../lhsExpressionInterface.hpp"
#include "../logic/constructStaticContext.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Defines an aggregated type via an array and implements the ExpressionInterface.
   *
   * See AggregatedActiveType for details.
   *
   * @tparam T_Real Real value of the aggregated type.
   * @tparam T_InnerActiveType CoDiPack type that composes the aggregated type.
   * @tparam T_Impl The final implementation of the aggregated type.
   * @tparam T_isStatic If the aggregated type is created in a static context.
   */
  template<typename T_Real, typename T_InnerActiveType, typename T_Impl, bool T_isStatic>
  struct AggregatedActiveTypeBase : public ExpressionInterface<T_Real, T_Impl> {
    public:

      using Real = T_Real;  ///< See AggregatedActiveTypeBase.
      using InnerActiveType = CODI_DD(T_InnerActiveType,
                                      CODI_DEFAULT_LHS_EXPRESSION);  ///< See AggregatedActiveTypeBase.
      using Impl = CODI_DD(T_Impl, AggregatedActiveTypeBase);        ///< See AggregatedActiveTypeBase.
      static bool constexpr isStatic = T_isStatic;                   ///< See AggregatedActiveTypeBase.

      using Tape = typename InnerActiveType::Tape;            ///< The tape of the inner active type.
      using Traits = RealTraits::AggregatedTypeTraits<Real>;  ///< The traits for the aggregated type.
      static int constexpr Elements = Traits::Elements;       ///< The number of elements in the aggregated type.

      using InnerReal = typename Traits::InnerType;       ///< Inner real type of the active type.
      using InnerIdentifier = typename Tape::Identifier;  ///< Identifier of the underlying tape.

      InnerActiveType values[Elements];  ///< Array representation.

      CODI_INLINE AggregatedActiveTypeBase() = default;                                 ///< Constructor.
      CODI_INLINE AggregatedActiveTypeBase(AggregatedActiveTypeBase const&) = default;  ///< Constructor.
      CODI_INLINE ~AggregatedActiveTypeBase() = default;                                ///< Destructor.

      /*******************************************************************************/
      /// Implementation of ExpressionInterface

      /// \copydoc codi::ExpressionInterface::StoreAs
      using StoreAs = typename std::conditional<isStatic, Impl, Impl const&>::type;
      using ADLogic = Tape;  ///< \copydoc codi::ExpressionInterface::ADLogic

      /// \copydoc codi::ExpressionInterface::getValue()
      CODI_INLINE Real const getValue() const {
        Real value{};
        static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
          Traits::template arrayAccess<i.value>(value) = values[i.value].getValue();
        });
        return value;
      }

      /// \copydoc codi::ExpressionInterface::applyTangent
      template<size_t argNumber>
      CODI_INLINE Real applyTangent(InnerReal const& tangent) const {
        CODI_UNUSED(tangent);
        return Real();  // TODO: return Traits::template tangentOfConstructor<argNumber>(getValue(), adjoint);
      }

      /// \copydoc codi::ExpressionInterface::applyAdjoint
      ///
      /// Calls Traits::adjointOfConstructor.
      template<size_t argNumber>
      CODI_INLINE InnerReal applyAdjoint(Real const& adjoint) const {
        return Traits::template adjointOfConstructor<argNumber>(getValue(), adjoint);
      }

      /// \copydoc codi::ComputeOperation::getMathRep
      static CODI_INLINE std::string getMathRep() {
        return Traits::getMathRep();
      }

      /*******************************************************************************/
      /// Implementation of NodeInterface

      static size_t constexpr LinkCount = Elements;  ///< \copydoc codi::NodeInterface::LinkCount

      /// \copydoc codi::NodeInterface::getLink
      template<size_t argNumber>
      CODI_INLINE InnerActiveType const& getLink() const {
        return values[argNumber];
      }

    protected:

      /// Cast to implementation.
      CODI_INLINE Impl const& cast() const {
        return static_cast<Impl const&>(*this);
      }

      /// Cast to implementation.
      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }
  };

  /**
   * @brief Represents a concrete aggregated lvalue int the CoDiPack expression tree.
   *
   * An aggregated type, in the sense of CoDiPack, is a structure that can be expressed by a set of double values. E.g.
   * std::complex<double> can be represented by two double values. The use case of this class is to be able to add such
   * types to the CoDiPack expression tree. If you use std::complex<codi::RealReverse>, then the expression tree for
   * operations on the complex number is usually quite short. When a and b are complex numbers with a CoDiPack type,
   * then <code>c = sin(a) + cos(b) </code> would record at least 6 statements. If complex numbers are
   * integrated into the expression tree, then only one or two statements are recorded.
   *
   * In order to add an aggregated type to the CoDiPack expression tree, the traits
   * codi::RealTraits::AggregatedTypeTraits need to be specialized for the aggregated type. The helper
   * codi::RealTraits::ArrayAggregatedTypeTraitsBase can be used if the aggregated type can be interpreted as an array
   * of values. In addition this class needs to be extended and special constructors and assignment operators need
   * to be implemented.
   *
   * An example is described in \ref Tutorial_07_Aggregated_type_implementation.
   *
   * @tparam T_Real Real value of the aggregated type. (E.g. std::complex<double>)
   * @tparam T_InnerActiveType CoDiPack type that composes the aggregated type. (E.g. codi::RealReverse)
   * @tparam T_Impl The final implementation of the aggregated type.
   */
  template<typename T_Real, typename T_InnerActiveType, typename T_Impl>
  struct AggregatedActiveType : public AggregatedActiveTypeBase<T_Real, T_InnerActiveType, T_Impl, false> {
    public:
      using Real = T_Real;                                                              ///< See AggregatedActiveType.
      using InnerActiveType = CODI_DD(T_InnerActiveType, CODI_DEFAULT_LHS_EXPRESSION);  ///< See AggregatedActiveType.
      using Impl = CODI_DD(T_Impl, AggregatedActiveType);                               ///< See AggregatedActiveType.

      using Base =
          AggregatedActiveTypeBase<T_Real, T_InnerActiveType, T_Impl, false>;  ///< Abbreviation for base class.
      using Traits = RealTraits::AggregatedTypeTraits<Real>;                   ///< The traits for the aggregated type.
      using PassiveReal = RealTraits::PassiveReal<Real>;                       ///< Passive value type of the real.

      using Base::Base;                               ///< Use base constructors.
      CODI_INLINE ~AggregatedActiveType() = default;  ///< Destructor.

      /// Constructor.
      template<typename Expr>
      CODI_INLINE AggregatedActiveType(ExpressionInterface<Real, Expr> const& expr) : Base() {
        store(expr, EventHints::Statement::Expression);
      }

      /// Constructor.
      CODI_INLINE AggregatedActiveType(AggregatedActiveType const& expr) : Base() {
        store(expr, EventHints::Statement::Copy);
      }

      /// Constructor.
      CODI_INLINE AggregatedActiveType(PassiveReal const& expr) : Base() {
        static_for<Base::Elements>([&](auto i) CODI_LAMBDA_INLINE {
          Base::values[i.value] = Traits::template arrayAccess<i.value>(expr);
        });
      }

      /// Assign operation.
      template<typename Expr>
      CODI_INLINE Impl& operator=(ExpressionInterface<Real, Expr> const& expr) {
        store(expr, EventHints::Statement::Expression);

        return Base::cast();
      }

      /// Assign operation.
      CODI_INLINE Impl& operator=(AggregatedActiveType const& expr) {
        store(expr, EventHints::Statement::Copy);

        return Base::cast();
      }

      /// Assign operation.
      CODI_INLINE Impl& operator=(PassiveReal const& expr) {
        static_for<Base::Elements>([&](auto i) CODI_LAMBDA_INLINE {
          Base::values[i.value] = Traits::template arrayAccess<i.value>(expr);
        });

        return Base::cast();
      }

    protected:

      /// \copydoc codi::InternalStatementRecordingTapeInterface::store()
      template<typename Rhs>
      CODI_INLINE void store(ExpressionInterface<Real, Rhs> const& rhs, EventHints::Statement const& eventType) {
        static_for<Traits::Elements>([&](auto i) {
          InnerActiveType& lhs = Base::values[i.value];
          EventSystem<typename InnerActiveType::Tape>::notifyStatementPrimalListeners(
              InnerActiveType::getTape(), lhs.getValue(), lhs.getIdentifier(),
              Traits::template arrayAccess<i.value>(rhs.cast().getValue()), eventType);
        });

        InnerActiveType::getTape().store(*this, rhs.cast());
      }
  };
}
