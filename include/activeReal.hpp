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

namespace codi {
  template<typename Real, typename Tape>
  class ActiveReal : public Expression<Real, ActiveReal<Real, Tape> > {
  public:
    typedef Real RealType;
    typedef Tape TapeType;
    typedef typename Tape::GradientData GradientData;

    static Tape globalTape;

  private:
    Real value;
    GradientData gradientData;

  public:

    inline ActiveReal() : value() {
      globalTape.initGradientData(value, gradientData);
    }

    inline ActiveReal(const Real& value) : value(value) {
      globalTape.initGradientData(this->value, gradientData);
    }

    inline ActiveReal(const Real& value, const Real& gradient) : value(value) {
      globalTape.initGradientData(this->value, gradientData);
      globalTape.setGradient(gradientData, gradient);
    }

    template<class R>
    inline ActiveReal(const Expression<Real, R>& rhs) {
      globalTape.store(value, gradientData, rhs.cast());
    }

    inline ActiveReal(const ActiveReal<Real, Tape>& v) {
      globalTape.store(value, gradientData, v);
    }

    inline ~ActiveReal() {
      globalTape.destroyGradientData(value, gradientData);
    }

    inline void calcGradient(Real& gradient) const {
      globalTape.pushJacobi(gradient, value, gradientData);
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      globalTape.pushJacobi(gradient, multiplier, value, gradientData);
    }

    inline GradientData& getGradientData() {
      return gradientData;
    }

    inline const GradientData& getGradientData() const {
      return gradientData;
    }

    inline Real getGradient() const {
      return globalTape.getGradient(gradientData);
    }

    inline Real& getGradient() {
      return globalTape.getGradient(gradientData);
    }

    inline void setGradient(const Real& gradient) {
      globalTape.setGradient(gradientData, gradient);
    }

    inline Real getValue() const {
      return value;
    }

    inline void setValue(const Real& value) {
      this->value = value;
    }

    inline ActiveReal<Real, Tape>& operator=(const Real& rhs){
      globalTape.store(value, gradientData, rhs);
      return *this;
    }

    template<class R>
    inline ActiveReal<Real, Tape>& operator=(const Expression<Real, R>& rhs){
      globalTape.store(value, gradientData, rhs.cast());
      return *this;
    }

    inline ActiveReal<Real, Tape>& operator=(const ActiveReal<Real, Tape>& rhs) {
      globalTape.store(value, gradientData, rhs);
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

    inline ActiveReal<Real, Tape>& operator+=(const Real& rhs) {
      // Optimization of code: If jacobies would be stored an identity operation is produced on the tape
      value += rhs;
      return *this;
    }
    inline ActiveReal<Real, Tape>& operator-=(const Real& rhs) {
      // Optimization of code: If jacobies would be stored an identity operation is produced on the tape
      value -= rhs;
      return *this;
    }
    inline ActiveReal<Real, Tape>& operator*=(const Real& rhs) {
      return *this = (*this * rhs);
    }
    inline ActiveReal<Real, Tape>& operator/=(const Real& rhs) {
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
  };

  template<typename Real, typename Tape>
  Tape ActiveReal<Real, Tape>::globalTape;
}
