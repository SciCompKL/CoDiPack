#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../tapes/interfaces/gradientAccessTapeInterface.hpp"
#include "../../traits/realTraits.hpp"
#include "../expressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Replacement type of ConstantExpression types in ConstructStaticContext.
   *
   * See ConstructStaticContext for detailed information.
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam _Real    Original primal value of the statement/expression.
   * @tparam _offset  The offset into the primal vector during construction.
   */
  template<typename _Real, size_t _offset>
  struct StaticContextConstantExpression : public ExpressionInterface<_Real, StaticContextConstantExpression<_Real, _offset>> {
      // TODO: Delete and use ConstantExpression instead.
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double); ///< See StaticContextConstantExpression
      static size_t constexpr offset = CODI_DECLARE_DEFAULT(_offset, 0); ///< See StaticContextConstantExpression

    private:

        Real const primal;
    public:

      /// Constructor
      CODI_INLINE StaticContextConstantExpression(Real const* const primalVector) :
        primal(primalVector[offset])
      {}

      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      using StoreAs = StaticContextConstantExpression; ///< codi::ExpressionInterface::StoreAs

      /// \copydoc codi::ExpressionInterface::getValue()
      CODI_INLINE Real const getValue() const {
        return primal;
      }

      /// \copydoc codi::ExpressionInterface::getJacobian()
      template<size_t argNumber>
      CODI_INLINE Real getJacobian() const {
        return Real();
      }

      /// @}
      /*******************************************************************************/
      /// @name Implementation of NodeInterface
      /// @{

      static bool constexpr EndPoint = true; ///< \copydoc codi::NodeInterface::EndPoint

      /// \copydoc codi::NodeInterface::forEachLink()
      template<typename Logic, typename ... Args>
      CODI_INLINE void forEachLink(TraversalLogic<Logic>& logic, Args&& ... args) const {
        CODI_UNUSED(logic, args...);
      }

      /// \copydoc codi::NodeInterface::forEachLinkConstExpr()
      template<typename Logic, typename ... Args>
      CODI_INLINE static typename Logic::ResultType constexpr forEachLinkConstExpr(Args&& ... CODI_UNUSED_ARG(args)) {
        return Logic::NeutralElement;
      }

      /// @}

    private:
      StaticContextConstantExpression& operator=(StaticContextConstantExpression const&) = delete;
  };
}
