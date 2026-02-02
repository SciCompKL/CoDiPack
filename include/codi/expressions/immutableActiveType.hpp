/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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
#include "activeType.hpp"
#include "assignmentOperators.hpp"
#include "incrementOperators.hpp"
#include "lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Creates a pseudo active type from a data value. Can be used to overlay existing data with immutable active
   * types.
   *
   * The class stores copies to the value and identifier. The identifier is taken as it is and not initialized or
   * destroyed. The class only wraps the data in a CoDiPack expression.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam T_ActiveType  The type of the active type which is wrapped.
   */
  template<typename T_ActiveType>
  struct ImmutableActiveType
      : public LhsExpressionInterface<typename T_ActiveType::Real, typename T_ActiveType::Gradient,
                                      typename T_ActiveType::Tape, ImmutableActiveType<T_ActiveType>>,
        public AssignmentOperators<typename T_ActiveType::Tape::Real, T_ActiveType::Tape::AllowJacobianOptimization,
                                   ImmutableActiveType<T_ActiveType>>,
        public IncrementOperators<typename T_ActiveType::Tape, ImmutableActiveType<T_ActiveType>> {
    public:

      using ActiveType = CODI_DD(T_ActiveType, CODI_T(ActiveType<CODI_DEFAULT_TAPE>));  ///< See ImmutableActiveType.
      using Tape = typename ActiveType::Tape;                                           ///< See ActiveType.

      using Real = typename Tape::Real;                    ///< See LhsExpressionInterface.
      using PassiveReal = RealTraits::PassiveReal<Real>;   ///< Basic computation type.
      using Identifier = typename Tape::Identifier;        ///< See LhsExpressionInterface.
      using Gradient = typename Tape::Gradient;            ///< See LhsExpressionInterface.
      using TapeData = typename Tape::ActiveTypeTapeData;  ///< See IdentifierInformationTapeInterface.

      using Base = LhsExpressionInterface<Real, Gradient, Tape, ImmutableActiveType>;  ///< Base class abbreviation.

    private:

      Real const primalValue;
      TapeData const tapeData;

    public:

      /// The tape data is not initialized. It is assumed to be a valid tape data (either default or assigned by an
      /// expression) and has to be valid throughout the lifespan of this object.
      CODI_INLINE ImmutableActiveType(Real const& value, TapeData const& tapeData)
          : primalValue(value), tapeData(tapeData) {
        // deliberately left empty
      }

      /// Create an immutable copy of an active type. It is assumed that the tape data is valid throughout the lifespan
      /// of this object.
      CODI_INLINE ImmutableActiveType(ActiveType const& value)
          : primalValue(value.getValue()), tapeData(value.getTapeData()) {
        // deliberately left empty
      }

      /// The tape data is not destroyed. It is assumed to be still valid, since this is only an immutable copy of
      /// the actual value.
      CODI_INLINE ~ImmutableActiveType() {
        // deliberately left empty
      }

      /// This class is immutable, delete all assignment operators.
      CODI_INLINE ImmutableActiveType<ActiveType>& operator=(ImmutableActiveType<ActiveType> const& v) = delete;

      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      using StoreAs = ImmutableActiveType const&;  ///< \copydoc codi::ExpressionInterface::StoreAs
      using ADLogic = Tape;                        ///< \copydoc codi::ExpressionInterface::ADLogic

      /// @}
      /*******************************************************************************/
      /// @name Implementation of LhsExpressionInterface
      ///
      /// Only the const access functions are implemented.
      ///
      /// @{

      /// \copydoc codi::LhsExpressionInterface::getIdentifier() const
      CODI_INLINE Identifier const& getIdentifier() const {
        return getTape().getIdentifier(tapeData);
      }

      /// \copydoc codi::LhsExpressionInterface::getTapeData() const
      CODI_INLINE TapeData const& getTapeData() const {
        return tapeData;
      }

      /// \copydoc codi::LhsExpressionInterface::value() const
      CODI_INLINE Real const& value() const {
        return primalValue;
      }

      /// \copydoc codi::LhsExpressionInterface::getTape()
      static CODI_INLINE Tape& getTape() {
        return ActiveType::getTape();
      }

      /// @}
  };
}  // end of namespace codi
