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

#include "../misc/macros.hpp"
#include "../config.h"
#include "../tapes/interfaces/fullTapeInterface.hpp"
#include "../traits/realTraits.hpp"
#include "assignmentOperators.hpp"
#include "incrementOperators.hpp"
#include "lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Represents a concrete lvalue in the CoDiPack expression tree.
   *
   * The class uses members for storing the value and the identifier.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam T_Tape  The tape that manages all expressions created with this type.
   */
  template<typename T_Tape>
  struct ActiveType
      : public LhsExpressionInterface<typename T_Tape::Real, typename T_Tape::Gradient, T_Tape, ActiveType<T_Tape>>,
        public AssignmentOperators<T_Tape, ActiveType<T_Tape>>,
        public IncrementOperators<T_Tape, ActiveType<T_Tape>> {
    public:

      /// See ActiveType.
      /// For reverse AD, the tape must implement ReverseTapeInterface.
      /// For forward AD, the 'tape' (that is not a tape, technically) must implement
      /// InternalStatementRecordingTapeInterface and GradientAccessTapeInterface.
      using Tape = CODI_DD(T_Tape, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));

      using Real = typename Tape::Real;                   ///< See LhsExpressionInterface.
      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.
      using Identifier = typename Tape::Identifier;       ///< See LhsExpressionInterface.
      using Gradient = typename Tape::Gradient;           ///< See LhsExpressionInterface.

      using Base = LhsExpressionInterface<Real, Gradient, Tape, ActiveType>;  ///< Base class abbreviation.

    private:

      Real primalValue;
      Identifier identifier;

      static Tape tape;

    public:

      /// Constructor
      CODI_INLINE ActiveType() : primalValue(), identifier() {
        Base::init();
      }

      /// Constructor
      CODI_INLINE ActiveType(ActiveType<Tape> const& v) : primalValue(), identifier() {
        Base::init();
        this->getTape().store(*this, v);
      }

      /// Constructor
      CODI_INLINE ActiveType(Real const& value) : primalValue(value), identifier() {
        Base::init();
      }

      /// Constructor
      template<typename U = Real, typename = RealTraits::EnableIfNotPassiveReal<U>>
      CODI_INLINE ActiveType(PassiveReal const& value) :
        primalValue(value), identifier() {
        Base::init();
      }

      /// Constructor
      template<typename Rhs>
      CODI_INLINE ActiveType(ExpressionInterface<Real, Rhs> const& rhs) : primalValue(), identifier() {
        Base::init();
        this->getTape().store(*this, rhs.cast());
      }

      /// Constructor
      template<typename Rhs, typename U = Real, typename = RealTraits::EnableIfNotPassiveReal<U>>
      CODI_INLINE ActiveType(ExpressionInterface<typename U::Real, Rhs> const& rhs) : primalValue(rhs.cast()), identifier() {
        Base::init();
      }

      /// Destructor
      CODI_INLINE ~ActiveType() {
        Base::destroy();
      }

      /// See LhsExpressionInterface::operator =(ExpressionInterface const&).
      CODI_INLINE ActiveType<Tape>& operator=(ActiveType<Tape> const& v) {
        static_cast<LhsExpressionInterface<Real, Gradient, Tape, ActiveType>&>(*this) = v;
        return *this;
      }
      using LhsExpressionInterface<Real, Gradient, Tape, ActiveType>::operator=;

      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      using StoreAs = ActiveType const&;  ///< \copydoc codi::ExpressionInterface::StoreAs
      using ActiveResult = ActiveType;    ///< \copydoc codi::ExpressionInterface::ActiveResult

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

      /// \copydoc codi::LhsExpressionInterface::getTape()
      static CODI_INLINE Tape& getTape() {
        return tape;
      }

      /// @}
  };

  template<typename Tape>
  Tape ActiveType<Tape>::tape{};
}
