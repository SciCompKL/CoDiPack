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
#include "../../misc/macros.hpp"
#include "../../tapes/interfaces/gradientAccessTapeInterface.hpp"
#include "../../tapes/interfaces/internalStatementRecordingTapeInterface.hpp"
#include "../../traits/realTraits.hpp"
#include "../expressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Replacement type of LhsExpressionInterface types in ConstructStaticContext.
   *
   * See ConstructStaticContext for detailed information.
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam T_Tape  The tape that create the original expression.
   */
  template<typename T_Tape>
  struct StaticContextActiveType : public ExpressionInterface<typename T_Tape::Real, StaticContextActiveType<T_Tape>> {
    public:

      using Tape = CODI_DD(
          T_Tape, CODI_T(CODI_UNION<InternalStatementRecordingTapeInterface<int>,
                                    GradientAccessTapeInterface<double, int>>));  ///< See StaticContextActiveType.

      using Real = CODI_DD(typename Tape::Real, double);  ///< See TapeTypesInterface.
      using Identifier = typename Tape::Identifier;       ///< See TapeTypesInterface.

      using Base = ExpressionInterface<Real, StaticContextActiveType>;  ///< Base class abbreviation.

    private:

      Real const primal;
      Identifier const identifier;

    public:

      /// Constructor
      CODI_INLINE StaticContextActiveType(Real const& primal, Identifier const& identifier)
          : primal(primal), identifier(identifier) {}

      /// Copy constructor
      CODI_INLINE StaticContextActiveType(StaticContextActiveType const& other)
          : Base(static_cast<Base const&>(other)), primal(other.primal), identifier(other.identifier) {}

      /// Empty constructor for delayed construction that is overwritten with placement new operators
      CODI_INLINE StaticContextActiveType() : primal(), identifier() {}

      /*******************************************************************************/
      /// @name Partial implementation of LhsExpressionInterface
      /// @{

      /// \copydoc codi::LhsExpressionInterface::getIdentifier() const
      CODI_INLINE Identifier const& getIdentifier() const {
        return identifier;
      }

      /// @}
      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      using StoreAs = StaticContextActiveType;  ///< \copydoc codi::ExpressionInterface::StoreAs
      using ADLogic = Tape;                     ///< \copydoc codi::ExpressionInterface::ADLogic

      /// \copydoc codi::ExpressionInterface::getValue()
      CODI_INLINE Real const getValue() const {
        return primal;
      }

      /// @}
      /*******************************************************************************/
      /// @name Implementation of NodeInterface
      /// @{

      static size_t constexpr LinkCount = 0;  ///< \copydoc codi::NodeInterface::LinkCount

      /// @}

    private:
      StaticContextActiveType& operator=(StaticContextActiveType const&) = delete;
  };

  /// Specialization of ActiveResultImpl for active types in a static context.
  template<typename T_Real, typename T_Tape>
  struct ExpressionTraits::ActiveResultImpl<T_Real, T_Tape, true> {
    public:

      using Real = CODI_DD(T_Real, CODI_ANY);  ///< See ExpressionTraits::ActiveResultImpl.
      using Tape = CODI_DD(T_Tape, CODI_ANY);  ///< See ExpressionTraits::ActiveResultImpl.

      /// The resulting active type of an expression.
      using ActiveResult = StaticContextActiveType<Tape>;
  };
}
