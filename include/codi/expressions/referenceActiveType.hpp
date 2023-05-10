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
  struct ReferenceActiveType : public LhsExpressionInterface<typename T_Type::Real, typename T_Type::Gradient,
                                                             typename T_Type::Tape, ReferenceActiveType<T_Type>>,
                               public AssignmentOperators<T_Type, ReferenceActiveType<T_Type>>,
                               public IncrementOperators<T_Type, ReferenceActiveType<T_Type>> {
    public:

      /// See ReferenceActiveType.
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);
      using Tape = typename Type::Tape;  ///< See LhsExpressionInterface.

      using Real = typename Tape::Real;                   ///< See LhsExpressionInterface.
      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.
      using Identifier = typename Tape::Identifier;       ///< See LhsExpressionInterface.
      using Gradient = typename Tape::Gradient;           ///< See LhsExpressionInterface.

    private:

      Type& reference;

    public:

      // TODO: Implement const variant

      /// Used by Jacobian tapes to optimize for reoccurring arguments.
      mutable Real jacobian;

      /// Constructor
      CODI_INLINE ReferenceActiveType(Type& v) : reference(v), jacobian() {}

      /// See LhsExpressionInterface::operator=(ExpressionInterface const&).
      CODI_INLINE ReferenceActiveType<Tape>& operator=(ReferenceActiveType<Tape> const& v) {
        static_cast<LhsExpressionInterface<Real, Gradient, Tape, ReferenceActiveType>&>(*this) = v;
        return *this;
      }
      using LhsExpressionInterface<Real, Gradient, Tape, ReferenceActiveType>::operator=;

      /*******************************************************************************/
      /// @name Implementation of LhsExpressionInterface
      /// @{

      using StoreAs = ReferenceActiveType const&;        ///< \copydoc codi::ExpressionInterface::StoreAs
      using ActiveResult = typename Type::ActiveResult;  ///< \copydoc codi::ExpressionInterface::ActiveResult

      /// \copydoc codi::LhsExpressionInterface::getIdentifier()
      CODI_INLINE Identifier& getIdentifier() {
        return reference.getIdentifier();
      }

      /// \copydoc codi::LhsExpressionInterface::getIdentifier() const
      CODI_INLINE Identifier const& getIdentifier() const {
        return reference.getIdentifier();
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
      static CODI_INLINE Tape& getTape() {
        return Type::getTape();
      }
  };
}
