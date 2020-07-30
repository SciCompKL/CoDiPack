#pragma once

#include "../aux/macros.hpp"
#include "../config.h"
#include "expressionInterface.hpp"
#include "logic/compileTimeTraversalLogic.hpp"
#include "logic/nodeInterface.hpp"
#include "logic/traversalLogic.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real>
  struct UnaryOperation {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg);

      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result);
  };


  template<typename _Real, typename _Arg, template<typename> class _Operation>
  struct UnaryExpression : public ExpressionInterface<_Real, UnaryExpression<_Real, _Arg, _Operation> > {
    public:
      using Real = CODI_DECLARE_DEFAULT(_Real, double);
      using Arg = CODI_DECLARE_DEFAULT(_Arg, CODI_TEMPLATE(ExpressionInterface<double, CODI_ANY>));
      using Operation = CODI_DECLARE_DEFAULT(CODI_TEMPLATE(_Operation<Real>), CODI_TEMPLATE(UnaryOperation<Real>));

      static bool constexpr EndPoint = false;
      using StoreAs = UnaryExpression;

    private:

      typename Arg::StoreAs arg;
      Real result;

    public:

      template<typename RealArg>
      explicit UnaryExpression(const ExpressionInterface<RealArg, Arg>& arg) :
        arg(arg.cast()),
        result(Operation::primal(this->arg.getValue())) {}


      /****************************************************************************
       * Section: Implementation of ExpressionInterface functions
       */

      CODI_INLINE Real const& getValue() const {
        return result;
      }

      template<size_t argNumber>
      CODI_INLINE Real getJacobian() const {
        return Operation::gradient(arg.getValue(), result);
      }

      /****************************************************************************
       * Section: Implementation of NodeInterface functions
       */

      template<typename Logic, typename ... Args>
      CODI_INLINE void forEachLink(TraversalLogic<Logic>& logic, Args&& ... args) const {
        logic.cast().template link<0>(arg, *this, std::forward<Args>(args)...);
      }

      template<typename CompileTimeLogic, typename ... Args>
      CODI_INLINE static typename CompileTimeLogic::ResultType constexpr forEachLinkConst(Args&& ... args) {
        return CompileTimeLogic::template link<0, Arg, UnaryExpression>(std::forward<Args>(args)...);
      }
  };
}
