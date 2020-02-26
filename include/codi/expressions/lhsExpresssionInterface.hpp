#pragma once

#include "../aux/macros.h"
#include "../config.h"
#include "../tapes/interfaces/internalExpressionTapeInterface.hpp"
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
      using Tape = DECLARE_DEFAULT(_Tape, TEMPLATE(InternalExpressionTapeInterface<ANY>));
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

      CODI_INLINE void initBase() {
        Impl::getGlobalTape().initIdentifier(cast().value(), cast().getIdentifier());
      }

      CODI_INLINE ~LhsExpressionInterface() {
        Impl::getGlobalTape().destroyIdentifier(cast().value(), cast().getIdentifier());
      }

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
      CODI_INLINE static typename Logic::ResultType constexpr forEachLinkConst(Args&& ... args) {
        CODI_UNUSED(args...);
      }
  };

  template<typename _Type>
  struct RealTraits<_Type, enableIfLhsExpression<typename _Type::Real, typename _Type::Gradient, typename _Type::Tape, _Type>> {
    public:

      using Type = DECLARE_DEFAULT(
                      _Type,
                      TEMPLATE(LhsExpressionInterface<double, double, InternalExpressionTapeInterface<ANY>, _Type>)
                    );
      using Real = typename Type::Real;

      using PassiveReal = PassiveRealType<Real>;

      static int constexpr MaxDerivativeOrder = 1 + RealTraits<Real>::MaxDerivativeOrder;

      static CODI_INLINE PassiveReal const& getPassiveValue(Type const& v) {
        return RealTraits<Real>::getPassiveValue(v.getValue());
      }
  };
}
