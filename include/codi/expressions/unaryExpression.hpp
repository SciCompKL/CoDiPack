#pragma once

#include "../aux/macros.h"
#include "../config.h"
#include "expressionInterface.hpp"
#include "logic/nodeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real>
  struct UnaryOperation {

      using Real = DECLARE_DEFAULT(_Real, double);

      template<typename Arg>
      static CODI_INLINE Real primal(const Arg& arg);

      template<typename Arg>
      static CODI_INLINE Real gradient(const Arg& arg, const Real& result);
  };


  template<typename _Real, typename _Arg, template<typename> class _Operation>
  struct UnaryExpression : public ExpressionInterface<_Real, UnaryExpression<_Real, _Arg, _Operation> > {
    public:
      using Real = DECLARE_DEFAULT(_Real, double);
      using Arg = DECLARE_DEFAULT(_Arg, TEMPLATE(ExpressionInterface<double, ANY>));
      using Operation = DECLARE_DEFAULT(_Operation, UnaryOperation);

    private:

      typename Arg::StoreAs arg;
      Real result;

    public:

      explicit UnaryExpression(const ExpressionInterface<Real, Arg>& arg) :
        arg(arg.cast()),
        result(Operation<Real>::primal(this->arg.getValue())) {}


      /****************************************************************************
       * Section: Implementation of ExpressionInterface functions
       */

      CODI_INLINE const Real& getValue() const {
        return result;
      }

      template<size_t argNumber>
      CODI_INLINE Real getJacobian() const {
        return Operation<Real>::gradient(arg.getValue(), result);
      }
  };
}
