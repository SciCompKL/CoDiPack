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

#include "../config.h"
#include "../misc/macros.hpp"
#include "../tapes/interfaces/fullTapeInterface.hpp"
#include "../traits/realTraits.hpp"
#include "assignmentOperators.hpp"
#include "incrementOperators.hpp"
#include "lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Represents the base implementation concrete lvalue in the CoDiPack expression tree.
   *
   * The class uses members for storing the value and the identifier.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * The storage of the underlying tape and the access to it is left to the implementing class.
   *
   * @tparam T_Tape  The tape that manages all expressions created with this type.
   * @tparam T_Impl  Implementing class.
   */
  template<typename T_Tape, typename T_Impl>
  struct ActiveTypeBase
      : public LhsExpressionInterface<typename T_Tape::Real, typename T_Tape::Gradient, T_Tape, T_Impl>,
        public AssignmentOperators<T_Tape, T_Impl>,
        public IncrementOperators<T_Tape, T_Impl> {
    public:

      /// See ActiveTypeBase.
      /// For reverse AD, the tape must implement ReverseTapeInterface.
      /// For forward AD, the 'tape' (that is not a tape, technically) must implement
      /// InternalStatementRecordingTapeInterface and GradientAccessTapeInterface.
      using Tape = CODI_DD(T_Tape, CODI_DEFAULT_TAPE);

      /// Abbreviation for the implementing class.
      using Impl = CODI_DD(T_Impl, CODI_DEFAULT_LHS_EXPRESSION);

      using Real = typename Tape::Real;                   ///< See LhsExpressionInterface.
      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.
      using Identifier = typename Tape::Identifier;       ///< See LhsExpressionInterface.
      using Gradient = typename Tape::Gradient;           ///< See LhsExpressionInterface.

      using Base = LhsExpressionInterface<Real, Gradient, T_Tape, T_Impl>;  ///< Base class abbreviation.

    private:

      Real primalValue;
      Identifier identifier;

    public:

      /// Cast to the implementation.
      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

      /// Constructor
      CODI_INLINE ActiveTypeBase() : primalValue(), identifier() {
        Base::init(Real(), EventHints::Statement::Passive);
      }

      /// Constructor
      CODI_INLINE ActiveTypeBase(ActiveTypeBase const& v) : primalValue(), identifier() {
        Base::init(v.getValue(), EventHints::Statement::Copy);
        cast().getTape().store(*this, v);
      }

      /// Constructor
      CODI_INLINE ActiveTypeBase(Real const& value) : primalValue(value), identifier() {
        Base::init(value, EventHints::Statement::Passive);
      }

      /// Constructor
      template<typename U = Real, typename = RealTraits::EnableIfNotPassiveReal<U>>
      CODI_INLINE ActiveTypeBase(PassiveReal const& value) : primalValue(value), identifier() {
        Base::init(value, EventHints::Statement::Passive);
      }

      /// Constructor
      template<typename Rhs>
      CODI_INLINE ActiveTypeBase(ExpressionInterface<Real, Rhs> const& rhs) : primalValue(), identifier() {
        Base::init(rhs.cast().getValue(), EventHints::Statement::Expression);
        cast().getTape().store(*this, rhs.cast());
      }

      /// Constructor
      template<typename Rhs, typename U = Real, typename = RealTraits::EnableIfNotPassiveReal<U>>
      CODI_INLINE ActiveTypeBase(ExpressionInterface<typename U::Real, Rhs> const& rhs)
          : primalValue(rhs.cast()), identifier() {
        Base::init(rhs.cast().getValue(), EventHints::Statement::Passive);
      }

      /// Destructor
      CODI_INLINE ~ActiveTypeBase() {
        Base::destroy();
      }

      /// See LhsExpressionInterface::operator=(LhsExpressionInterface const&).
      CODI_INLINE Impl& operator=(ActiveTypeBase const& v) {
        static_cast<LhsExpressionInterface<Real, Gradient, Tape, Impl>&>(*this) = v;
        return cast();
      }
      using Base::operator=;

      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      using StoreAs = Impl const&;  ///< \copydoc codi::ExpressionInterface::StoreAs
      using ActiveResult = Impl;    ///< \copydoc codi::ExpressionInterface::ActiveResult

      /// @}
      /*******************************************************************************/
      /// @name Implementation of LhsExpressionInterface
      /// @{

      /// \copydoc codi::LhsExpressionInterface::getIdentifier()
      CODI_INLINE Identifier& getIdentifier() {
        return identifier;
      }

      /// \copydoc codi::LhsExpressionInterface::getIdentifier() const
      CODI_INLINE Identifier const& getIdentifier() const {
        return identifier;
      }

      /// \copydoc codi::LhsExpressionInterface::value()
      CODI_INLINE Real& value() {
        return primalValue;
      }

      /// \copydoc codi::LhsExpressionInterface::value() const
      CODI_INLINE Real const& value() const {
        return primalValue;
      }

      // getTape must be implemented in the derived type

      /// @}
  };
}
