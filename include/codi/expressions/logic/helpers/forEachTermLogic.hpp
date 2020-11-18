#pragma once

#include <type_traits>
#include <utility>

#include "../../../aux/macros.hpp"
#include "../../../config.h"
#include "../../../traits/expressionTraits.hpp"
#include "../traversalLogic.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Impl>
  struct ForEachTermLogic : public TraversalLogic<_Impl> {
    public:

      using Impl = CODI_DECLARE_DEFAULT(_Impl, CODI_TEMPLATE(TraversalLogic<CODI_ANY>));

      /*******************************************************************************
       * Section: Methods the child class can overwrite
       *
       * Description: TODO
       *
       */

      /// TODO
      template<typename Node, typename ... Args>
      void handleActive(Node const& node, Args&& ... args) {
        CODI_UNUSED(node, args...);
      }

      /// TODO
      template<typename Node, typename ... Args>
      void handleConstant(Node const& node, Args&& ... args) {
        CODI_UNUSED(node, args...);
      }

      /*******************************************************************************
       * Section: Implementation of the handling switch
       *
       * Description: TODO
       */

      template<typename Node, typename ... Args>
      CODI_INLINE enableIfLhsExpression<Node> term(Node const& node, Args&& ... args) {
        cast().handleActive(node, std::forward<Args>(args)...);
      }
      template<typename Node, typename ... Args>
      CODI_INLINE enableIfConstantExpression<Node> term(Node const& node, Args&& ... args) {
        cast().handleConstant(node, std::forward<Args>(args)...);
      }
      using TraversalLogic<Impl>::term;

    private:

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

  };
}
