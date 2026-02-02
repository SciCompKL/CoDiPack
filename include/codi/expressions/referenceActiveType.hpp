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
#include "../tapes/interfaces/gradientAccessTapeInterface.hpp"
#include "../traits/realTraits.hpp"
#include "assignmentOperators.hpp"
#include "incrementOperators.hpp"
#include "lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Holds a reference to an ActiveType for manual optimization of common arguments.
   *
   * See the \ref Example_14_ReferenceActiveType for an example use.
   *
   * @tparam T_Type  The type of the reference which is captured.
   */
  template<typename T_Type>
  struct ReferenceActiveType
      : public LhsExpressionInterface<typename T_Type::Real, typename T_Type::Gradient, typename T_Type::Tape,
                                      ReferenceActiveType<T_Type>>,
        public AssignmentOperators<typename T_Type::Tape::Real, T_Type::Tape::AllowJacobianOptimization,
                                   ReferenceActiveType<T_Type>>,
        public IncrementOperators<typename T_Type::Tape, ReferenceActiveType<T_Type>> {
    public:

      /// See ReferenceActiveType.
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);
      using Tape = typename Type::Tape;  ///< See LhsExpressionInterface.

      using Real = typename Tape::Real;                    ///< See LhsExpressionInterface.
      using PassiveReal = RealTraits::PassiveReal<Real>;   ///< Basic computation type.
      using Identifier = typename Tape::Identifier;        ///< See LhsExpressionInterface.
      using Gradient = typename Tape::Gradient;            ///< See LhsExpressionInterface.
      using TapeData = typename Tape::ActiveTypeTapeData;  ///< See IdentifierInformationTapeInterface.

    protected:

      Type& reference;  ///< Reference to the underlying active type.

    public:

      // TODO: Implement const variant

      /// Used by Jacobian tapes to optimize for reoccurring arguments.
      mutable Real jacobian;

      /// Constructor
      CODI_INLINE ReferenceActiveType(Type& v) : reference(v), jacobian() {}

      /// Copy constructor
      CODI_INLINE ReferenceActiveType(ReferenceActiveType const& o) : reference(o.reference), jacobian() {}

      /// See LhsExpressionInterface::operator=(ExpressionInterface const&).
      CODI_INLINE ReferenceActiveType& operator=(ReferenceActiveType const& v) {
        static_cast<LhsExpressionInterface<Real, Gradient, Tape, ReferenceActiveType>&>(*this) = v;
        return *this;
      }
      using LhsExpressionInterface<Real, Gradient, Tape, ReferenceActiveType>::operator=;

      /*******************************************************************************/
      /// @name Implementation of LhsExpressionInterface
      /// @{

      using StoreAs = ReferenceActiveType const&;  ///< \copydoc codi::ExpressionInterface::StoreAs
      using ADLogic = Tape;                        ///< \copydoc codi::ExpressionInterface::ADLogic

      /// \copydoc codi::LhsExpressionInterface::getIdentifier()
      CODI_INLINE Identifier& getIdentifier() {
        return reference.getIdentifier();
      }

      /// \copydoc codi::LhsExpressionInterface::getIdentifier() const
      CODI_INLINE Identifier const& getIdentifier() const {
        return reference.getIdentifier();
      }

      /// \copydoc codi::LhsExpressionInterface::getTapeData()
      CODI_INLINE TapeData& getTapeData() {
        return reference.getTapeData();
      }

      /// \copydoc codi::LhsExpressionInterface::getTapeData() const
      CODI_INLINE TapeData const& getTapeData() const {
        return reference.getTapeData();
      }

      /// \copydoc codi::LhsExpressionInterface::value()
      CODI_INLINE Real& value() {
        return reference.value();
      }

      /// \copydoc codi::LhsExpressionInterface::value() const
      CODI_INLINE Real const& value() const {
        return reference.value();
      }

      /// \copydoc codi::LhsExpressionInterface::getTape()
      static CODI_INLINE decltype(Type::getTape()) getTape() {
        return Type::getTape();
      }
  };
}
