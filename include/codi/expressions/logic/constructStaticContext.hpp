#pragma once

#include "../../config.h"
#include "../../aux/macros.h"
#include "../../traits/expressionTraits.hpp"
#include "../binaryExpression.hpp"
#include "../expressionInterface.hpp"
#include "../static/staticContextActiveType.hpp"
#include "../static/staticContextConstantExpression.hpp"
#include "../unaryExpression.hpp"
#include "nodeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Rhs, typename _Tape, size_t _primalValueOffset, size_t _constantValueOffset, typename = void>
  struct ConstructStaticContextLogic {
    public:

      using Rhs = _Rhs;
      using Tape = _Tape;
      static constexpr size_t primalValueOffset = _primalValueOffset;
      static constexpr size_t constantValueOffset = _constantValueOffset;

      using Real = typename Tape::Real;
      using Identifier = typename Tape::Identifier;
      using PassiveReal = typename Tape::PassiveReal;

      using ResultType = DECLARE_DEFAULT(_Rhs, TEMPLATE(ExpressionInterface<double, ANY>));

      static ResultType construct(Real* primalVector, Identifier const* const identifiers, PassiveReal const* const constantVector);
  };

  template<typename _Rhs, typename _Tape, size_t _primalValueOffset, size_t _constantValueOffset>
  struct ConstructStaticContextLogic<_Rhs, _Tape, _primalValueOffset, _constantValueOffset, enableIfLhsExpression<_Rhs>> {
    public:

      using Rhs = _Rhs;
      using Tape = _Tape;
      static constexpr size_t primalValueOffset = _primalValueOffset;
      static constexpr size_t constantValueOffset = _constantValueOffset;

      using Real = typename Tape::Real;
      using Identifier = typename Tape::Identifier;
      using PassiveReal = typename Tape::PassiveReal;

      using ResultType = StaticContextActiveType<Tape, primalValueOffset>;

      static ResultType construct(Real* primalVector, Identifier const* const identifiers, PassiveReal const* const constantVector) {
        CODI_UNUSED(constantVector);

        return ResultType(primalVector, identifiers);
      }
  };

  template<typename _Rhs, typename _Tape, size_t _primalValueOffset, size_t _constantValueOffset>
  struct ConstructStaticContextLogic<_Rhs, _Tape, _primalValueOffset, _constantValueOffset, enableIfConstantExpression<_Rhs>> {
    public:

      using Rhs = _Rhs;
      using Tape = _Tape;
      static constexpr size_t primalValueOffset = _primalValueOffset;
      static constexpr size_t constantValueOffset = _constantValueOffset;

      using Real = typename Tape::Real;
      using Identifier = typename Tape::Identifier;
      using PassiveReal = typename Tape::PassiveReal;

      using ResultType = StaticContextConstantExpression<PassiveReal, constantValueOffset>;

      static ResultType construct(Real* primalVector, Identifier const* const identifiers, PassiveReal const* const constantVector) {
        CODI_UNUSED(primalVector, identifiers);

        return ResultType(constantVector);
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

      using ArgAConstructor = ConstructStaticContextLogic<ArgA, Tape, primalValueOffset, constantValueOffset>;
      using ArgAMod = typename ArgAConstructor::ResultType;

      static size_t constexpr primalValueOffsetArgB = primalValueOffset + NumberOfActiveTypeArguments<ArgA>::value;
      static size_t constexpr constantValueOffsetArgB = constantValueOffset + NumberOfConstantTypeArguments<ArgA>::value;
      using ArgBConstructor = ConstructStaticContextLogic<ArgB, Tape, primalValueOffsetArgB, constantValueOffsetArgB>;
      using ArgBMod = typename ArgBConstructor::ResultType;

      using ResultType = BinaryExpression<OpReal, ArgAMod, ArgBMod, Operation>;

      static ResultType construct(Real* primalVector, Identifier const* const identifiers, PassiveReal const* const constantVector) {
        return ResultType(
              ArgAConstructor::construct(primalVector, identifiers, constantVector),
              ArgBConstructor::construct(primalVector, identifiers, constantVector));
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

      using ArgConstructor = ConstructStaticContextLogic<Arg, Tape, primalValueOffset, constantValueOffset>;
      using ArgMod = typename ArgConstructor::ResultType;

      using ResultType = UnaryExpression<OpReal, ArgMod, Operation>;

      static ResultType construct(Real* primalVector, Identifier const* const identifiers, PassiveReal const* const constantVector) {
        return ResultType(ArgConstructor::construct(primalVector, identifiers, constantVector));
      }
  };
}
