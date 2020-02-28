#pragma once

#include "../aux/macros.h"
#include "../config.h"
#include "../tapes/interfaces/gradientAccessTapeInterface.hpp"
#include "../traits/realTraits.hpp"
#include "assignmentOperators.hpp"
#include "incerementOperators.hpp"
#include "lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Tape>
  class ActiveType : public LhsExpressionInterface<typename _Tape::Real, typename _Tape::Gradient, _Tape, ActiveType<_Tape> >,
                     public AssignmentOperators<_Tape, ActiveType<_Tape>>,
                     public IncrementOperators<_Tape, ActiveType<_Tape>> {
    public:

      using Tape = DECLARE_DEFAULT(_Tape, TEMPLATE(GradientAccessTapeInterface<double, int>));

      using StoreAs = ActiveType const&;

      using Real = typename Tape::Real;
      using PassiveReal = PassiveRealType<Real>;
      using Identifier = typename Tape::Identifier;
      using Gradient = typename Tape::Gradient;

  private:

      Real primalValue;
      Identifier identifier;

      static Tape globalTape;

  public:

    ActiveType() : primalValue(), identifier() {
      this->initBase();
    }

    ActiveType(ActiveType<Tape> const& v) : primalValue(), identifier() {
      this->initBase();
      this->getGlobalTape().store(*this, v);
    }

    CODI_INLINE ActiveType(PassiveReal const& value) : primalValue(value), identifier() {
      this->initBase();
    }

    template<class Rhs>
    CODI_INLINE ActiveType(ExpressionInterface<Real, Rhs> const& rhs) : primalValue(), identifier() {
      this->initBase();
      this->getGlobalTape().store(*this, rhs.cast());
    }

    ActiveType<Tape>& operator=(ActiveType<Tape > const& v) {
      static_cast<LhsExpressionInterface<Real, Gradient, Tape, ActiveType>&>(*this) = v;
      return *this;
    }
    using LhsExpressionInterface<Real, Gradient, Tape, ActiveType>::operator=;

    CODI_INLINE Identifier& getIdentifier() {
      return identifier;
    }

    CODI_INLINE Identifier const& getIdentifier() const {
      return identifier;
    }

    CODI_INLINE Real& value() {
      return primalValue;
    }

    CODI_INLINE Real const& value() const {
      return primalValue;
    }

    static CODI_INLINE Tape& getGlobalTape() {
      return globalTape;
    }
  };

  template<typename Tape>
  Tape ActiveType<Tape>::globalTape = Tape();
}
