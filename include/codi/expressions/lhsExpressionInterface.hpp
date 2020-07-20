#pragma once

#include "../aux/macros.h"
#include "../config.h"
#include "../tapes/interfaces/internalStatementRecordingInterface.hpp"
#include "../traits/expressionTraits.hpp"
#include "../traits/realTraits.hpp"
#include "expressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Gradient, typename _Tape, typename _Impl>
  struct LhsExpressionInterface : public ExpressionInterface<_Real, _Impl> {
    public:

      using Real = DECLARE_DEFAULT(_Real, double);
      using Gradient = DECLARE_DEFAULT(_Gradient, Real);
      using Tape = DECLARE_DEFAULT(_Tape, TEMPLATE(InternalStatementRecordingInterface<ANY>));
      using Impl = DECLARE_DEFAULT(_Impl, LhsExpressionInterface);

      using Identifier = typename Tape::Identifier;
      using PassiveReal = PassiveRealType<Real>;

      static bool constexpr EndPoint = true;

      /*******************************************************************************
       * Section: Start of interface definition
       *
       */

      Real const& value() const;
      Real& value();

      Identifier& getIdentifier();
      Identifier const& getIdentifier() const;

      static Tape& getGlobalTape();

      /*******************************************************************************
       * Section: Start of general implementation
       *
       */

      using ExpressionInterface<Real, Impl>::cast;
      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

      CODI_INLINE Gradient& gradient() {
        return Impl::getGlobalTape().gradient(cast().getIdentifier());
      }

      CODI_INLINE Gradient const& gradient() const {
        return const_cast<Tape const&>(Impl::getGlobalTape()).gradient(cast().getIdentifier());
      }

      CODI_INLINE Gradient getGradient() const {
        return cast().gradient();
      }

      CODI_INLINE void setGradient(Gradient const& g) {
        cast().gradient() = g;
      }

      CODI_INLINE Real const& getValue() const {
        return cast().value();
      }

      CODI_INLINE void setValue(Real const& v) {
        cast().value() = v;
      }

      CODI_INLINE Impl& operator=(PassiveReal const& rhs){
        Impl::getGlobalTape().store(cast(), rhs);
        return cast();
      }

      template<typename Rhs>
      CODI_INLINE Impl& operator=(ExpressionInterface<Real, Rhs> const& rhs){
        Impl::getGlobalTape().store(cast(), rhs.cast());
        return cast();
      }

      CODI_INLINE Impl& operator=(LhsExpressionInterface const& rhs) {
        Impl::getGlobalTape().store(cast(), rhs);
        return cast();
      }

      template<typename Logic, typename ... Args>
      CODI_INLINE void forEachLink(TraversalLogic<Logic>& logic, Args&& ... args) const {
        CODI_UNUSED(logic, args...);
      }

      template<typename Logic, typename ... Args>
      CODI_INLINE static typename Logic::ResultType constexpr forEachLinkConstExpr(Args&& ... args) {
        CODI_UNUSED(args...);
      }

    protected:

      CODI_INLINE void init() {
        Impl::getGlobalTape().initIdentifier(cast().value(), cast().getIdentifier());
      }

      CODI_INLINE void destroy() {
        Impl::getGlobalTape().destroyIdentifier(cast().value(), cast().getIdentifier());
      }
  };

  template<typename _Type>
  struct RealTraits<_Type, enableIfLhsExpression<_Type>> {
    public:

      using Type = DECLARE_DEFAULT(
                      _Type,
                      TEMPLATE(LhsExpressionInterface<double, double, InternalStatementRecordingInterface<ANY>, _Type>)
                    );
      using Real = typename Type::Real;

      using PassiveReal = PassiveRealType<Real>;

      static int constexpr MaxDerivativeOrder = 1 + RealTraits<Real>::MaxDerivativeOrder;

      static CODI_INLINE PassiveReal const& getPassiveValue(Type const& v) {
        return RealTraits<Real>::getPassiveValue(v.getValue());
      }

      static CODI_INLINE bool isTotalZero(Type const& v) {
        // TODO: Specialize for forward and reverse tapes.
        return Real() == v.getValue() && Type::Gradient() == v.getGradient();
      }
  };
}
