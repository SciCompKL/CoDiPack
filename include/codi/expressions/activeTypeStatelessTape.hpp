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
   * @brief Represents a concrete lvalue in the CoDiPack expression tree.
   *
   * See also LhsExpressionInterface.
   *
   * This active type does not work with a fixed tape. Instead, getTape() constructs a new temporary tape on every call.
   * In particular, tapes for this active type can not have a persistent state.
   *
   * @tparam T_Tape  The tape that manages all expressions created with this type.
   */
  template<typename T_Tape>
  struct ActiveTypeStatelessTape : public LhsExpressionInterface<typename T_Tape::Real, typename T_Tape::Gradient,
                                                                 T_Tape, ActiveTypeStatelessTape<T_Tape>>,
                                   public AssignmentOperators<typename T_Tape::Real, T_Tape::AllowJacobianOptimization,
                                                              ActiveTypeStatelessTape<T_Tape>>,
                                   public IncrementOperators<T_Tape, ActiveTypeStatelessTape<T_Tape>> {
    public:

      using Tape = CODI_DD(T_Tape, CODI_DEFAULT_TAPE);  ///< See ActiveTypeStatelessTape.

      using Real = typename Tape::Real;                    ///< See LhsExpressionInterface.
      using PassiveReal = RealTraits::PassiveReal<Real>;   ///< Basic computation type.
      using Identifier = typename Tape::Identifier;        ///< See LhsExpressionInterface.
      using Gradient = typename Tape::Gradient;            ///< See LhsExpressionInterface.
      using TapeData = typename Tape::ActiveTypeTapeData;  ///< See IdentifierInformationTapeInterface.

      using Base = LhsExpressionInterface<Real, Gradient, T_Tape, ActiveTypeStatelessTape<T_Tape>>;  ///< Base class
                                                                                                     ///< abbreviation.

    private:

      Real primalValue;
      TapeData tapeData;

    public:

      /// @brief Constructor
      /// @details CUDA compiler has problems when this function is annotated with \c __device__.
      constexpr CODI_INLINE_NO_FA ActiveTypeStatelessTape() = default;

      /// Constructor
      constexpr CODI_INLINE ActiveTypeStatelessTape(PassiveReal const& value) : primalValue(value), tapeData() {}

      /// Constructor
      CODI_INLINE ActiveTypeStatelessTape(ActiveTypeStatelessTape const& v) : primalValue(), tapeData() {
        Base::init(v.getValue(), EventHints::Statement::Copy);
        getTape().store(*this, v);
      }

      /// Constructor
      template<typename Rhs>
      CODI_INLINE ActiveTypeStatelessTape(ExpressionInterface<Real, Rhs> const& rhs) : primalValue(), tapeData() {
        Base::init(rhs.cast().getValue(), EventHints::Statement::Expression);
        getTape().store(*this, rhs.cast());
      }

      /*******************************************************************************/
      /// @name Assignment operators (all forwarding to the base class)
      /// @{

      /// See ActiveTypeStatelessTape::operator=(ActiveTypeStatelessTape const&).
      CODI_INLINE ActiveTypeStatelessTape& operator=(ActiveTypeStatelessTape const& v) {
        static_cast<Base&>(*this) = static_cast<Base const&>(v);
        return *this;
      }

      using Base::operator=;

      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      using StoreAs = ActiveTypeStatelessTape const&;  ///< \copydoc codi::ExpressionInterface::StoreAs
      using ADLogic = Tape;                            ///< \copydoc codi::ExpressionInterface::ADLogic

      /// @}
      /*******************************************************************************/
      /// @name Implementation of LhsExpressionInterface
      /// @{

      /// \copydoc codi::LhsExpressionInterface::getIdentifier()
      CODI_INLINE Identifier& getIdentifier() {
        return getTape().getIdentifier(tapeData);
      }

      /// \copydoc codi::LhsExpressionInterface::getIdentifier() const
      CODI_INLINE Identifier const& getIdentifier() const {
        return getTape().getIdentifier(tapeData);
      }

      /// \copydoc codi::LhsExpressionInterface::getTapeData()
      CODI_INLINE TapeData& getTapeData() {
        return tapeData;
      }

      /// \copydoc codi::LhsExpressionInterface::getTapeData() const
      CODI_INLINE TapeData const& getTapeData() const {
        return tapeData;
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
      static CODI_INLINE Tape getTape() {
        return Tape();
      }

      /// @}
  };
}
