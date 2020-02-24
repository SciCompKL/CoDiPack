#pragma once

#include "../aux/macros.h"
#include "../config.h"
#include "expressionInterface.hpp"
#include "logic/nodeInterface.hpp"

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
      using Operation = DECLARE_DEFAULT(_Operation, BinaryOperation);

    private:

      typename ArgA::StoreAs argA;
      typename ArgB::StoreAs argB;
      Real result;

    public:

      explicit BinaryExpression(ExpressionInterface<Real, ArgA> const& argA, ExpressionInterface<Real, ArgB> const& argB) :
        argA(argA.cast()),
        argB(argB.cast()),
        result(Operation<Real>::primal(this->argA.getValue(), this->argB.getValue())) {}


      /****************************************************************************
       * Section: Implementation of ExpressionInterface functions
       */

      CODI_INLINE Real const& getValue() const {
        return result;
      }

      template<size_t argNumber>
      CODI_INLINE Real getJacobian() const {
        if(0 == linkNumber) {
          return Operation<Real>::gradientA(argA.getValue(), argB.getValue(), result);
        } else {
          return Operation<Real>::gradientB(argA.getValue(), argB.getValue(), result);
        }
      }
  };
}
