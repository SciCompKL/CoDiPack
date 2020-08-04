#pragma once

#include <iostream>

#include "../aux/macros.hpp"
#include "../config.h"
#include "../tapes/interfaces/gradientAccessTapeInterface.hpp"
#include "../tapes/interfaces/internalStatementRecordingInterface.hpp"
#include "../traits/expressionTraits.hpp"
#include "../traits/realTraits.hpp"
#include "expressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Gradient, typename _Tape, typename _Impl>
  struct LhsExpressionInterface : public ExpressionInterface<_Real, _Impl> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);
      using Gradient = CODI_DECLARE_DEFAULT(_Gradient, Real);
      using Tape = CODI_DECLARE_DEFAULT(_Tape, CODI_TEMPLATE(CODI_UNION<InternalStatementRecordingInterface<int>, GradientAccessTapeInterface<double, int>>));
      using Impl = CODI_DECLARE_DEFAULT(_Impl, LhsExpressionInterface);

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
      CODI_INLINE static typename Logic::ResultType constexpr forEachLinkConstExpr(Args&& ... CODI_UNUSED_ARG(args)) {
        return Logic::NeutralElement;
      }

    protected:

      CODI_INLINE void init() {
        Impl::getGlobalTape().initIdentifier(cast().value(), cast().getIdentifier());
      }

      CODI_INLINE void destroy() {
        Impl::getGlobalTape().destroyIdentifier(cast().value(), cast().getIdentifier());
      }
  };

  template<typename Real, typename Gradient, typename Tape, typename Impl>
  std::ostream& operator<<(std::ostream& out, LhsExpressionInterface<Real, Gradient, Tape, Impl> const& v) {
    return out << v.cast().value();
  }

  template<typename _Type>
  struct RealTraits<_Type, enableIfLhsExpression<_Type>> {
    public:

      using Type = CODI_DECLARE_DEFAULT(
                      _Type,
                      CODI_TEMPLATE(LhsExpressionInterface<double, double, InternalStatementRecordingInterface<CODI_ANY>, _Type>)
                    );
      using Real = typename Type::Real;

      using PassiveReal = PassiveRealType<Real>;

      static int constexpr MaxDerivativeOrder = 1 + RealTraits<Real>::MaxDerivativeOrder;

      static CODI_INLINE PassiveReal const& getPassiveValue(Type const& v) {
        return RealTraits<Real>::getPassiveValue(v.getValue());
      }

      static CODI_INLINE bool isTotalZero(Type const& v) {
        // TODO: Specialize for forward and reverse tapes.
        return Real() == v.getValue() && typename Type::Gradient() == v.getGradient();
      }
  };
}
