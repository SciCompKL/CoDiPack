#pragma once

#include "../../config.h"
#include "../../aux/macros.hpp"
#include "../../tapes/interfaces/reverseTapeInterface.hpp"
#include "../../traits/expressionTraits.hpp"
#include "../binaryExpression.hpp"
#include "../expressionInterface.hpp"
#include "../static/staticContextActiveType.hpp"
#include "../unaryExpression.hpp"
#include "nodeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Helper class for the construction of an expression in a different context.
   *
   * Converts the termination nodes of the expression into the static context replacements. The initialization is
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
   * @tparam _Rhs  The expression type. Needs to implement the expression ExpressionInterface.
   * @tparam _Tape  The tape which stored the expression (Used for the definition of the primal tapes
   *
   */
  template<typename _Rhs, typename _Tape, size_t _primalValueOffset, size_t _constantValueOffset, typename = void>
  struct ConstructStaticContextLogic {
    public:

      using Rhs = CODI_DD(_Rhs, CODI_T(ExpressionInterface<double, CODI_ANY>)); ///< See ConstructStaticContextLogic
      using Tape = CODI_DD(_Tape, CODI_T(ReverseTapeInterface<double, double, CODI_ANY>)); ///< See ConstructStaticContextLogic
      static constexpr size_t primalValueOffset = CODI_DD(_primalValueOffset, 0); ///< See ConstructStaticContextLogic
      static constexpr size_t constantValueOffset = CODI_DD(_constantValueOffset, 0); ///< See ConstructStaticContextLogic

      using Real = typename Tape::Real;               ///< See TapeTypesInterface.
      using Identifier = typename Tape::Identifier;   ///< See TapeTypesInterface.
      using PassiveReal = typename Tape::PassiveReal; ///< Basic computation type

      /// The resulting expression type after all nodes are replaced.
      using ResultType = CODI_DD(_Rhs, CODI_T(ExpressionInterface<double, CODI_ANY>));

      /**
       * @brief Perform the construction
       *
       * See ConstructStaticContextLogic on how the arguments are used and which conversion are performed.
       */
      static ResultType construct(Real* primalVector, Identifier const* const identifiers, PassiveReal const* const constantData);
  };

#ifndef DOXYGEN_DISABLE

  template<typename _Rhs, typename _Tape, size_t _primalValueOffset, size_t _constantValueOffset>
  struct ConstructStaticContextLogic<_Rhs, _Tape, _primalValueOffset, _constantValueOffset, ExpressionTraits::EnableIfLhsExpression<_Rhs>> {
    public:

      using Rhs = _Rhs;
      using Tape = _Tape;
      static constexpr size_t primalValueOffset = _primalValueOffset;
      static constexpr size_t constantValueOffset = _constantValueOffset;

      using Real = typename Tape::Real;
      using Identifier = typename Tape::Identifier;
      using PassiveReal = typename Tape::PassiveReal;

      /// Conversion from LhsExpressionInterface to StaticContextActiveType
      using ResultType = StaticContextActiveType<Tape>;

      /// Uses primalVector[identifiers[primalValueOffset]] and identifiers[primalValueOffset]
      static ResultType construct(Real* primalVector, Identifier const* const identifiers, PassiveReal const* const constantData) {
        CODI_UNUSED(constantData);

        Identifier const identifier = identifiers[primalValueOffset];
        Real const primal = primalVector[identifier];

        return ResultType(primal, identifier);
      }
  };

  template<typename _Rhs, typename _Tape, size_t _primalValueOffset, size_t _constantValueOffset>
  struct ConstructStaticContextLogic<_Rhs, _Tape, _primalValueOffset, _constantValueOffset, ExpressionTraits::EnableIfConstantExpression<_Rhs>> {
    public:

      using Rhs = _Rhs;
      using Tape = _Tape;
      static constexpr size_t primalValueOffset = _primalValueOffset;
      static constexpr size_t constantValueOffset = _constantValueOffset;

      using Real = typename Tape::Real;
      using Identifier = typename Tape::Identifier;
      using PassiveReal = typename Tape::PassiveReal;

      /// Conversion from ConstantExpression to ConstantExpression
      using ResultType = ConstantExpression<PassiveReal>;

      /// Uses constantData[constantValueOffset]
      static ResultType construct(Real* primalVector, Identifier const* const identifiers, PassiveReal const* const constantData) {
        CODI_UNUSED(primalVector, identifiers);

        return ResultType(constantData[constantValueOffset]);
      }
  };

  template<typename _Real, typename _ArgA, typename _ArgB, template<typename> class _Operation, typename _Tape, size_t _primalValueOffset, size_t _constantValueOffset>
  struct ConstructStaticContextLogic< BinaryExpression<_Real, _ArgA, _ArgB, _Operation>, _Tape, _primalValueOffset, _constantValueOffset> {
    public:

      using OpReal = _Real;
      using ArgA = _ArgA;
      using ArgB = _ArgB;
      template<typename T> using Operation = _Operation<T>;
      using Tape = _Tape;
      static constexpr size_t primalValueOffset = _primalValueOffset;
      static constexpr size_t constantValueOffset = _constantValueOffset;

      using Real = typename Tape::Real;
      using Identifier = typename Tape::Identifier;
      using PassiveReal = typename Tape::PassiveReal;

      /// Unmodified offsets for first argument
      using ArgAConstructor = ConstructStaticContextLogic<ArgA, Tape, primalValueOffset, constantValueOffset>;
      using ArgAMod = typename ArgAConstructor::ResultType;

      /// Shift offsets by the number of occurrences in the first sub tree.
      static size_t constexpr primalValueOffsetArgB = primalValueOffset + ExpressionTraits::NumberOfActiveTypeArguments<ArgA>::value;
      static size_t constexpr constantValueOffsetArgB = constantValueOffset + ExpressionTraits::NumberOfConstantTypeArguments<ArgA>::value;
      using ArgBConstructor = ConstructStaticContextLogic<ArgB, Tape, primalValueOffsetArgB, constantValueOffsetArgB>;
      using ArgBMod = typename ArgBConstructor::ResultType;

      using ResultType = BinaryExpression<OpReal, ArgAMod, ArgBMod, Operation>;

      static ResultType construct(Real* primalVector, Identifier const* const identifiers, PassiveReal const* const constantData) {
        return ResultType(
              ArgAConstructor::construct(primalVector, identifiers, constantData),
              ArgBConstructor::construct(primalVector, identifiers, constantData));
      }
  };

  template<typename _Real, typename _Arg, template<typename> class _Operation, typename _Tape, size_t _primalValueOffset, size_t _constantValueOffset>
  struct ConstructStaticContextLogic< UnaryExpression<_Real, _Arg, _Operation>, _Tape, _primalValueOffset, _constantValueOffset> {
    public:

      using OpReal = _Real;
      using Arg = _Arg;
      template<typename T> using Operation = _Operation<T>;
      using Tape = _Tape;
      static constexpr size_t primalValueOffset = _primalValueOffset;
      static constexpr size_t constantValueOffset = _constantValueOffset;

      using Real = typename Tape::Real;
      using Identifier = typename Tape::Identifier;
      using PassiveReal = typename Tape::PassiveReal;

      /// Unmodified offsets since there is just one sub tree
      using ArgConstructor = ConstructStaticContextLogic<Arg, Tape, primalValueOffset, constantValueOffset>;
      using ArgMod = typename ArgConstructor::ResultType;

      using ResultType = UnaryExpression<OpReal, ArgMod, Operation>;

      static ResultType construct(Real* primalVector, Identifier const* const identifiers, PassiveReal const* const constantData) {
        return ResultType(ArgConstructor::construct(primalVector, identifiers, constantData));
      }
  };
#endif
}
