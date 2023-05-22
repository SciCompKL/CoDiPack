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
#include "expressionInterface.hpp"
#include "logic/nodeInterface.hpp"
#include "logic/traversalLogic.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Helper class for the constant data conversion in primal value tapes.
  ///
  /// @tparam T_StoreData  Type of the constant data store in the tape.
  template<typename T_StoreData>
  struct ConstantDataConversion {
    public:

      using StoreData = CODI_DD(T_StoreData, double);  ///< See ConstantDataConversion.
      using ArgumentData = StoreData;                  ///< Defined by specializations.

      /// Convert the data from the store to the argument of the constant expression.
      static ArgumentData const& fromDataStore(StoreData const& v) {
        return v;
      }

      /// Convert the data from the constant expression to the store.
      static StoreData const& toDataStore(StoreData const& v) {
        return v;
      }
  };

  /**
   * @brief Represents constant values in the expression tree.
   *
   * All values that are not a CoDiPack type are considered constant, for example values like 4.0 or double a.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam T_Real  Original primal value of the statement/expression.
   * @tparam T_ConversionOperator  Functions for the conversion of the constant data for primal value tape stores.
   */
  template<typename T_Real, template<typename> class T_ConversionOperator = ConstantDataConversion>
  struct ConstantExpression : public ExpressionInterface<T_Real, ConstantExpression<T_Real, T_ConversionOperator>> {
    public:

      using Real = T_Real;  ///< See ConstantExpression.

      /// See ConstantExpression.
      template<typename T>
      using ConversionOperator = CODI_DD(CODI_T(T_ConversionOperator<T>), CODI_T(ConstantDataConversion<T>));

    private:
      Real primalValue;

    public:

      /// Constructor
      CODI_INLINE ConstantExpression(Real const& v) : primalValue(v) {}

      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      using StoreAs = ConstantExpression;  ///< \copydoc codi::ExpressionInterface::StoreAs
      using ActiveResult = void;           ///< \copydoc codi::ExpressionInterface::ActiveResult

      /// \copydoc codi::ExpressionInterface::getValue
      CODI_INLINE Real const& getValue() const {
        return primalValue;
      }

      /// \copydoc codi::ExpressionInterface::getJacobian()
      template<size_t argNumber>
      CODI_INLINE Real getJacobian() const {
        return Real();
      }

      /// @}
      /*******************************************************************************/
      /// @name Implementation of NodeInterface
      /// @{

      static bool constexpr EndPoint = true;  ///< \copydoc codi::NodeInterface::EndPoint

      /// \copydoc codi::NodeInterface::forEachLink()
      template<typename Logic, typename... Args>
      CODI_INLINE void forEachLink(TraversalLogic<Logic>& logic, Args&&... args) const {
        CODI_UNUSED(logic, args...);
      }

      /// \copydoc codi::NodeInterface::forEachLinkConstExpr()
      template<typename Logic, typename... Args>
      CODI_INLINE static typename Logic::ResultType constexpr forEachLinkConstExpr(Args&&... CODI_UNUSED_ARG(args)) {
        return Logic::NeutralElement;
      }

      /// @}
  };
}
