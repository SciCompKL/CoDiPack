#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../tapes/interfaces/gradientAccessTapeInterface.hpp"
#include "../../tapes/interfaces/internalStatementRecordingInterface.hpp"
#include "../../traits/realTraits.hpp"
#include "../expressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Tape>
  struct StaticContextActiveType : public ExpressionInterface<typename _Tape::Real, StaticContextActiveType<_Tape>> {
    public:

      using Tape = CODI_DECLARE_DEFAULT(_Tape, CODI_TEMPLATE(CODI_UNION<InternalStatementRecordingInterface<int>, GradientAccessTapeInterface<double, int>>));

      using Real = typename Tape::Real;
      using Identifier = typename Tape::Identifier;

    private:

      Real const primal;
      Identifier const identifier;

    public:

      CODI_INLINE StaticContextActiveType(Real const& primal, Identifier const& identifier) :
        primal(primal),
        identifier(identifier)
      {}

      CODI_INLINE Identifier const& getIdentifier() const {
        return identifier;
      }

      /*******************************************************************************
       * Section: ExpressionInterface implementation
       */

      using StoreAs = StaticContextActiveType;

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
      CODI_INLINE static typename Logic::ResultType constexpr forEachLinkConstExpr(Args&& ... CODI_UNUSED_ARG(args)) {
        return Logic::NeutralElement;
      }

    private:
      StaticContextActiveType& operator=(StaticContextActiveType const&) = delete;
  };
}
