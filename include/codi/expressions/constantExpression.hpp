#pragma once

#include "../aux/macros.hpp"
#include "../config.h"
#include "expressionInterface.hpp"
#include "logic/nodeInterface.hpp"
#include "logic/traversalLogic.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real>
  struct ConstantExpression : public ExpressionInterface<_Real, ConstantExpression<_Real>> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      using StoreAs = ConstantExpression;

      static bool constexpr EndPoint = true;

    private:
      Real primalValue;

    public:

      CODI_INLINE ConstantExpression(Real const& v) : primalValue(v) {}

      CODI_INLINE Real const& getValue() const {
        return primalValue;
      }

      template<size_t argNumber>
      CODI_INLINE Real getJacobian() const {
        return Real();
      }

      template<typename Logic, typename ... Args>
      CODI_INLINE void forEachLink(TraversalLogic<Logic>& logic, Args&& ... args) const {
        CODI_UNUSED(logic, args...);
      }

      template<typename Logic, typename ... Args>
      CODI_INLINE static typename Logic::ResultType constexpr forEachLinkConstExpr(Args&& ... CODI_UNUSED_ARG(args)) {
        return Logic::NeutralElement;
      }
  };
}
