#pragma once

#include "../aux/macros.h"
#include "../config.h"
#include "expressionInterface.hpp"
#include "logic/compileTimeTraversalLogic.hpp"
#include "logic/nodeInterface.hpp"
#include "logic/traversalLogic.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real>
  struct BinaryOperation {
    public:

      using Real = DECLARE_DEFAULT(_Real, double);

      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB);

      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientA(ArgA const& argA, ArgB const& argB, Real const& result);

      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientB(ArgA const& argA, ArgB const& argB, Real const& result);
  };


  template<typename _Real, typename _ArgA, typename _ArgB, template<typename> class _Operation>
  struct BinaryExpression : public ExpressionInterface<_Real, BinaryExpression<_Real, _ArgA, _ArgB, _Operation> > {
    public:
      using Real = DECLARE_DEFAULT(_Real, double);
      using ArgA = DECLARE_DEFAULT(_ArgA, TEMPLATE(ExpressionInterface<double, ANY>));
      using ArgB = DECLARE_DEFAULT(_ArgB, TEMPLATE(ExpressionInterface<double, ANY>));
      using Operation = DECLARE_DEFAULT(TEMPLATE(_Operation<Real>), TEMPLATE(BinaryOperation<Real>));

      static bool constexpr EndPoint = false;
      using StoreAs = BinaryExpression;

    private:

      typename ArgA::StoreAs argA;
      typename ArgB::StoreAs argB;
      Real result;

    public:

      explicit BinaryExpression(ExpressionInterface<Real, ArgA> const& argA, ExpressionInterface<Real, ArgB> const& argB) :
        argA(argA.cast()),
        argB(argB.cast()),
        result(Operation::primal(this->argA.getValue(), this->argB.getValue())) {}


      /****************************************************************************
       * Section: Implementation of ExpressionInterface functions
       */

      CODI_INLINE Real const& getValue() const {
        return result;
      }

      template<size_t argNumber>
      CODI_INLINE Real getJacobian() const {
        if(0 == argNumber) {
          return Operation::gradientA(argA.getValue(), argB.getValue(), result);
        } else {
          return Operation::gradientB(argA.getValue(), argB.getValue(), result);
        }
      }

      /****************************************************************************
       * Section: Implementation of NodeInterface functions
       */

      template<typename Logic, typename ... Args>
      CODI_INLINE void forEachLink(TraversalLogic<Logic>& logic, Args&& ... args) const {
        logic.cast().template link<0>(argA, *this, std::forward<Args>(args)...);
        logic.cast().template link<1>(argB, *this, std::forward<Args>(args)...);
      }

      template<typename CompileTimeLogic, typename ... Args>
      CODI_INLINE static typename CompileTimeLogic::ResultType constexpr forEachLinkConstExpr(Args&& ... args) {
        return CompileTimeLogic::reduce(
              CompileTimeLogic::template link<0, ArgA, BinaryExpression>(std::forward<Args>(args)...),
              CompileTimeLogic::template link<1, ArgB, BinaryExpression>(std::forward<Args>(args)...));
      }
  };
}
