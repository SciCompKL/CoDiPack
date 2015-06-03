/**
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing, TU Kaiserslautern
 *
 * This file is part of CoDiPack.
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
 * License, or (at your option) any later version.
 *
 * CoDiPack is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU
 * General Public License along with CoDiPack.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * Authors: TODO
 */

#pragma once

#include "expressions.h"
#include "typeTraits.hpp"
#include "expressionTraits.hpp"
#include <iostream>

namespace codi {

  /*
   * Just a helper class to remove the globalTape from the ActiveReal implementation.
   */
  template<typename Tape>
  class GlobalActiveRealData  {
    public:
      static Tape globalTape;
  };

  template<typename Real, typename Tape>
  class ActiveReal : public Expression<Real, ActiveReal<Real, Tape> > {
  public:
    typedef Real RealType;
    typedef Tape TapeType;
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;
    typedef typename Tape::GradientData GradientData;

  private:
    Real value;
    GradientData gradientData;

  public:

    inline ActiveReal() : value() {
      getGlobalTape().initGradientData(value, gradientData);
    }

    inline ActiveReal(const PassiveReal& value) : value(value) {
      getGlobalTape().initGradientData(this->value, gradientData);
    }

    inline ActiveReal(const Real& value, const Real& gradient) : value(value) {
      getGlobalTape().initGradientData(this->value, gradientData);
      getGlobalTape().setGradient(gradientData, gradient);
    }

    template<class R>
    inline ActiveReal(const Expression<Real, R>& rhs) {
      getGlobalTape().store(value, gradientData, rhs.cast());
    }

    inline ActiveReal(const ActiveReal<Real, Tape>& v) {
      getGlobalTape().store(value, gradientData, v);
    }

    inline ~ActiveReal() {
      getGlobalTape().destroyGradientData(value, gradientData);
    }

    inline void calcGradient(Real& gradient) const {
      getGlobalTape().pushJacobi(gradient, value, gradientData);
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      getGlobalTape().pushJacobi(gradient, multiplier, value, gradientData);
    }

    inline GradientData& getGradientData() {
      return gradientData;
    }

    inline const GradientData& getGradientData() const {
      return gradientData;
    }

    inline Real& gradient() {
      return getGlobalTape().gradient(gradientData);
    }


    inline Real getGradient() const {
      return getGlobalTape().getGradient(gradientData);
    }

    inline void setGradient(const Real& gradient) {
      getGlobalTape().setGradient(gradientData, gradient);
    }

    inline Real& getValue() {
      return value;
    }

    inline const Real& getValue() const {
      return value;
    }

    inline void setValue(const Real& value) {
      this->value = value;
    }

    inline ActiveReal<Real, Tape>& operator=(const PassiveReal& rhs){
      getGlobalTape().store(value, gradientData, rhs);
      return *this;
    }

    template<class R>
    inline ActiveReal<Real, Tape>& operator=(const Expression<Real, R>& rhs){
      getGlobalTape().store(value, gradientData, rhs.cast());
      return *this;
    }

    inline ActiveReal<Real, Tape>& operator=(const ActiveReal<Real, Tape>& rhs) {
      getGlobalTape().store(value, gradientData, rhs);
      return *this;
    }

    template<class R>
    inline ActiveReal<Real, Tape>& operator+=(const Expression<Real, R>& rhs) {
      return *this = (*this + rhs);
    }
    template<class R>
    inline ActiveReal<Real, Tape>& operator-=(const Expression<Real, R>& rhs) {
      return *this = (*this - rhs);
    }
    template<class R>
    inline ActiveReal<Real, Tape>& operator*=(const Expression<Real, R>& rhs) {
      return *this = (*this * rhs);
    }
    template<class R>
    inline ActiveReal<Real, Tape>& operator/=(const Expression<Real, R>& rhs) {
      return *this = (*this / rhs);
    }

    inline ActiveReal<Real, Tape>& operator+=(const PassiveReal& rhs) {
      // Optimization of code: If jacobies would be stored an identity operation is produced on the tape
      value += rhs;
      return *this;
    }
    inline ActiveReal<Real, Tape>& operator-=(const PassiveReal& rhs) {
      // Optimization of code: If jacobies would be stored an identity operation is produced on the tape
      value -= rhs;
      return *this;
    }
    inline ActiveReal<Real, Tape>& operator*=(const PassiveReal& rhs) {
      return *this = (*this * rhs);
    }
    inline ActiveReal<Real, Tape>& operator/=(const PassiveReal& rhs) {
      return *this = (*this / rhs);
    }

    inline ActiveReal<Real, Tape> operator++() {
      return *this = *this + 1.0;
    }
    inline ActiveReal<Real, Tape> operator++(int) {
      ActiveReal<Real, Tape> r(*this);
      *this = *this + 1.0;
      return r;
    }
    inline ActiveReal<Real, Tape> operator--() {
      return *this = *this - 1.0;
    }
    inline ActiveReal<Real, Tape> operator--(int) {
      ActiveReal<Real, Tape> r(*this);
      *this = *this - 1.0;
      return r;
    }

    static inline Tape& getGlobalTape() {
      return GlobalActiveRealData<Tape>::globalTape;
    }
  };

  template<typename Real, typename Tape>
  class TypeTraits<ActiveReal<Real, Tape> > {
    public:
      typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

      static const typename TypeTraits<Real>::PassiveReal getBaseValue(const ActiveReal<Real, Tape>& t) {
        return TypeTraits<Real>::getBaseValue(t.getValue());
      }
  };

  template<typename Real, typename Tape>
  struct ExpressionTraits<ActiveReal<Real, Tape> >  {
    static const size_t maxActiveVariables = 1;
  };

  template<typename Tape>
  Tape GlobalActiveRealData<Tape>::globalTape;

  template<typename Real, class R>
  std::ostream& operator<<(std::ostream& os, const Expression<Real, R>& rhs){
    os << rhs.getValue();
    return os;
  }
  template<typename Real, typename Tape>
  std::istream& operator>>(std::istream& os, ActiveReal<Real, Tape>& rhs){
    Real temp;
    os >> temp;
    rhs.setValue(temp);
    return os;
  }
}
