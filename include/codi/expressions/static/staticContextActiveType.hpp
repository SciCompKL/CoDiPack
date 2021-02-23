#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../tapes/interfaces/gradientAccessTapeInterface.hpp"
#include "../../tapes/interfaces/internalStatementRecordingInterface.hpp"
#include "../../traits/realTraits.hpp"
#include "../expressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Replacement type of LhsExpressionInterface types in ConstructStaticContext.
   *
   * See ConstructStaticContext for detailed information.
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam _Tape  The tape that create the original expression.
   */
  template<typename _Tape>
  struct StaticContextActiveType : public ExpressionInterface<typename _Tape::Real, StaticContextActiveType<_Tape>> {
    public:

      using Tape = CODI_DECLARE_DEFAULT(_Tape, CODI_TEMPLATE(
                                      CODI_UNION<InternalStatementRecordingInterface<int>,
                                                 GradientAccessTapeInterface<double, int>>)); ///< See StaticContextActiveType

      using Real = typename Tape::Real;             ///< See TapeTypesInterface.
      using Identifier = typename Tape::Identifier; ///< See TapeTypesInterface.

    private:

      Real const primal;
      Identifier const identifier;

    public:

      /// Constructor
      CODI_INLINE StaticContextActiveType(Real const& primal, Identifier const& identifier) :
        primal(primal),
        identifier(identifier)
      {}

      /*******************************************************************************/
      /// @name Partial implementation of LhsExpressionInterface
      /// @{

      /// \copydoc codi::LhsExpressionInterface::getIdentifier() const
      CODI_INLINE Identifier const& getIdentifier() const {
        return identifier;
      }

      /// @}
      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      using StoreAs = StaticContextActiveType; ///< \copydoc codi::ExpressionInterface::EndPoint

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

      /// \copydoc codi::NodeInterface::forEachLink
      template<typename Logic, typename ... Args>
      CODI_INLINE void forEachLink(TraversalLogic<Logic>& logic, Args&& ... args) const {
        CODI_UNUSED(logic, args...);
      }

      /// \copydoc codi::NodeInterface::forEachLinkConstExpr
      template<typename Logic, typename ... Args>
      CODI_INLINE static typename Logic::ResultType constexpr forEachLinkConstExpr(Args&& ... CODI_UNUSED_ARG(args)) {
        return Logic::NeutralElement;
      }

      /// @}

    private:
      StaticContextActiveType& operator=(StaticContextActiveType const&) = delete;
  };
}
