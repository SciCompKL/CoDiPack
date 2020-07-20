#pragma once

#include "../../aux/macros.h"
#include "../../config.h"
#include "../../tapes/interfaces/gradientAccessTapeInterface.hpp"
#include "../../traits/realTraits.hpp"
#include "../expressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, size_t _offset>
  struct StaticContextConstantExpression : public ExpressionInterface<_Real, StaticContextConstantExpression<_Real, _offset>> {
    public:

      using Real = DECLARE_DEFAULT(_Real, double);
      static size_t constexpr offset = DECLARE_DEFAULT(_offset, 0);

    private:

        Real const primal;
    public:

      CODI_INLINE StaticContextConstantExpression(Real const* const primalVector) :
        primal(primalVector[offset])
      {}

      /*******************************************************************************
       * Section: ExpressionInterface implementation
       */

      using StoreAs = StaticContextConstantExpression;

      CODI_INLINE Real const getValue() const {
        return primal;
      }

      template<size_t argNumber>
      CODI_INLINE Real getJacobian() const {
        return Real();
      }

      /*******************************************************************************
       * Section: NodeInterface implementation
       *
       */

      static bool constexpr EndPoint = true;

      template<typename Logic, typename ... Args>
      CODI_INLINE void forEachLink(TraversalLogic<Logic>& logic, Args&& ... args) const {
        CODI_UNUSED(logic, args...);
      }

      template<typename Logic, typename ... Args>
      CODI_INLINE static typename Logic::ResultType constexpr forEachLinkConstExpr(Args&& ... args) {
        CODI_UNUSED(args...);
      }

    private:
      StaticContextConstantExpression& operator=(StaticContextConstantExpression const&) = delete;
  };
}
