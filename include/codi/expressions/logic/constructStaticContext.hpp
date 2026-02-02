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

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../../tapes/interfaces/reverseTapeInterface.hpp"
#include "../../traits/expressionTraits.hpp"
#include "../computeExpression.hpp"
#include "../constantExpression.hpp"
#include "../emptyExpression.hpp"
#include "../expressionInterface.hpp"
#include "../static/staticContextActiveType.hpp"
#include "nodeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Helper class for the construction of an expression in a different context.
   *
   * Converts the leaf nodes of the expression into the static context replacements. The initialization is
   * performed via three arrays.
   *
   * Conversion and initialization is done for:
   *  - LhsExpressionInterface -> StaticContextActiveType: id = identifiers[primalValueOffset]
   *                                                       primal = primalVector[id]
   *  - ConstantExpression -> ConstantExpression: value = constantData[constantValueOffset]
   *
   * The offsets are computed from the corresponding expression traits NumberOfActiveTypeArguments and
   * NumberOfConstantTypeArguments. They are evaluated on each sub graph.
   *
   * @tparam T_Rhs  The expression type. Needs to implement ExpressionInterface.
   * @tparam T_Tape  The tape which stored the expression.
   *
   */
  template<typename T_Rhs, typename T_Tape, size_t T_primalValueOffset, size_t T_constantValueOffset, typename = void>
  struct ConstructStaticContextLogic {
    public:

      using Rhs = CODI_DD(T_Rhs, CODI_T(ExpressionInterface<double, CODI_ANY>));  ///< See ConstructStaticContextLogic.
      using Tape = CODI_DD(
          T_Tape, CODI_T(ReverseTapeInterface<double, double, CODI_ANY>));  ///< See ConstructStaticContextLogic.
      static constexpr size_t primalValueOffset =
          CODI_DD(T_primalValueOffset, 0);  ///< See ConstructStaticContextLogic.
      static constexpr size_t constantValueOffset =
          CODI_DD(T_constantValueOffset, 0);  ///< See ConstructStaticContextLogic.

      using Real = typename Tape::Real;                ///< See TapeTypesInterface.
      using Identifier = typename Tape::Identifier;    ///< See TapeTypesInterface.
      using PassiveReal = typename Tape::PassiveReal;  ///< Basic computation type.

      /// The resulting expression type after all nodes are replaced.
      using ResultType = CODI_DD(T_Rhs, CODI_T(ExpressionInterface<double, CODI_ANY>));

      /**
       * @brief Perform the construction
       *
       * See ConstructStaticContextLogic on how the arguments are used and which conversion are performed.
       */
      static ResultType construct(Real* primalVector, Identifier const* const identifiers,
                                  PassiveReal const* const constantData);
  };

#ifndef DOXYGEN_DISABLE

  template<typename T_Rhs, typename T_Tape, size_t T_primalValueOffset, size_t T_constantValueOffset>
  struct ConstructStaticContextLogic<T_Rhs, T_Tape, T_primalValueOffset, T_constantValueOffset,
                                     ExpressionTraits::EnableIfLhsExpression<T_Rhs>> {
    public:

      using Rhs = T_Rhs;
      using Tape = T_Tape;
      static constexpr size_t primalValueOffset = T_primalValueOffset;
      static constexpr size_t constantValueOffset = T_constantValueOffset;

      using Real = typename Tape::Real;
      using Identifier = typename Tape::Identifier;
      using PassiveReal = typename Tape::PassiveReal;

      /// Conversion from LhsExpressionInterface to StaticContextActiveType.
      using ResultType = StaticContextActiveType<Tape>;

      /// Uses primalVector[identifiers[primalValueOffset]] and identifiers[primalValueOffset].
      CODI_INLINE static ResultType construct(Real* primalVector, Identifier const* const identifiers,
                                              PassiveReal const* const constantData) {
        CODI_UNUSED(constantData);

        Identifier const identifier = identifiers[primalValueOffset];
        Real const primal = primalVector[identifier];

        return ResultType(primal, identifier);
      }
  };

  template<typename T_Rhs, typename T_Tape, size_t T_primalValueOffset, size_t T_constantValueOffset>
  struct ConstructStaticContextLogic<T_Rhs, T_Tape, T_primalValueOffset, T_constantValueOffset,
                                     ExpressionTraits::EnableIfConstantExpression<T_Rhs>> {
    public:

      using Rhs = T_Rhs;
      using Tape = T_Tape;
      static constexpr size_t primalValueOffset = T_primalValueOffset;
      static constexpr size_t constantValueOffset = T_constantValueOffset;

      using Real = typename Rhs::Real;
      using InnerReal = typename Tape::Real;
      using Identifier = typename Tape::Identifier;
      using PassiveInnerReal = typename Tape::PassiveReal;

      /// Conversion from ConstantExpression to ConstantExpression.
      using ResultType = ConstantExpression<Real>;

      /// Uses constantData[constantValueOffset].
      CODI_INLINE static ResultType construct(InnerReal* primalVector, Identifier const* const identifiers,
                                              PassiveInnerReal const* const constantData) {
        CODI_UNUSED(primalVector, identifiers);

        using ConversionOperator = typename Rhs::template ConversionOperator<PassiveInnerReal>;
        using AggregateTraits = RealTraits::AggregatedTypeTraits<Real>;

        Real value{};
        static_for<AggregateTraits::Elements>([&](auto i) CODI_LAMBDA_INLINE {
          AggregateTraits::template arrayAccess<i.value>(value) =
              ConversionOperator::fromDataStore(constantData[constantValueOffset + i.value]);
        });

        return ResultType(value);
      }
  };

  template<typename T_Real, template<typename> class T_Operation, typename T_Tape, size_t T_primalValueOffset,
           size_t T_constantValueOffset, typename... T_Args>
  struct ConstructStaticContextLogic<ComputeExpression<T_Real, T_Operation, T_Args...>, T_Tape, T_primalValueOffset,
                                     T_constantValueOffset> {
    public:

      using OpReal = T_Real;
      template<typename T>
      using Operation = T_Operation<T>;
      using Tape = T_Tape;
      static constexpr size_t primalValueOffset = T_primalValueOffset;
      static constexpr size_t constantValueOffset = T_constantValueOffset;
      using Args = std::tuple<T_Args...>;

      using Real = typename Tape::Real;
      using Identifier = typename Tape::Identifier;
      using PassiveReal = typename Tape::PassiveReal;

      // Recursively compute the primal value offset. The template parameter T prevents a full template
      // specialization.
      template<typename T, size_t ArgNum>
      struct PrimalValueOffsetForArg
          : public std::integral_constant<size_t, PrimalValueOffsetForArg<T, ArgNum - 1>::value +
                                                      ExpressionTraits::NumberOfActiveTypeArguments<
                                                          typename std::tuple_element<ArgNum - 1, Args>::type>::value> {
      };

      // Primal value offset for the first argument. Here we return the global offset of the construction logic.
      template<typename T>
      struct PrimalValueOffsetForArg<T, 0> : public std::integral_constant<size_t, primalValueOffset> {};

      // Recursively compute the constant value offset. The template parameter T prevents a full template
      // specialization.
      template<typename T, size_t ArgNum>
      struct ConstantValueOffsetForArg
          : public std::integral_constant<size_t, ConstantValueOffsetForArg<T, ArgNum - 1>::value +
                                                      ExpressionTraits::NumberOfConstantTypeArguments<
                                                          typename std::tuple_element<ArgNum - 1, Args>::type>::value> {
      };

      // Constant value offset for the first argument. Here we return the global offset of the construction logic.
      template<typename T>
      struct ConstantValueOffsetForArg<T, 0> : public std::integral_constant<size_t, constantValueOffset> {};

      // Static context constructor for the given argument.
      template<size_t ArgNum>
      using ArgConstructor = ConstructStaticContextLogic<std::tuple_element_t<ArgNum, Args>, Tape,
                                                         PrimalValueOffsetForArg<double, ArgNum>::value,
                                                         ConstantValueOffsetForArg<double, ArgNum>::value>;

      // Modified argument type for the result expression.
      template<size_t ArgNum>
      using ArgMod = typename ArgConstructor<ArgNum>::ResultType;

      // Helper for ResultType definition. This allows us to expand on the index sequence Is.
      template<size_t... Is>
      static ComputeExpression<OpReal, T_Operation /* Make clang happy */, ArgMod<Is>...> ResultTypeHelper(
          std::index_sequence<Is...>);

      // Definition of the return type as the return value of ResultTypeHelper.
      using ResultType = remove_all<decltype(ResultTypeHelper(std::index_sequence_for<T_Args...>()))>;

      // Helper for construct definition. This allows us to expand on the index sequence Is.
      template<size_t... Is>
      CODI_INLINE static ResultType constructHelper(Real* primalVector, Identifier const* const identifiers,
                                                    PassiveReal const* const constantData, std::index_sequence<Is...>) {
        CODI_UNUSED(primalVector, identifiers, constantData);  // Required for empty sequence Is
        return ResultType(ArgConstructor<Is>::construct(primalVector, identifiers, constantData)...);
      }

      CODI_INLINE static ResultType construct(Real* primalVector, Identifier const* const identifiers,
                                              PassiveReal const* const constantData) {
        return constructHelper(primalVector, identifiers, constantData, std::index_sequence_for<T_Args...>());
      }
  };
#endif
}
