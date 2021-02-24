/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the
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
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *     Max Sagebaum
 *     Tim Albring
 *     Johannes Blühdorn
 */

#pragma once

#include "activeReal.hpp"
#include "tapes/tapeTraits.hpp"
#include "tools/direction.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  template<typename Real, typename = void>
  class Atomic {
    public:
      Atomic() = delete;
  };

  template<typename _Real>
  class Atomic<_Real, typename std::enable_if<std::is_same<_Real, float>::value || std::is_same<_Real, double>::value>::type> {
    public:
      using Real = _Real;

    private:
      Real value;

      CODI_INLINE void setValue(const Real& newValue) {

        #pragma omp atomic write
        this->value = newValue;
      }

      CODI_INLINE void addValue(const Real& increment) {

        #pragma omp atomic update
        this->value += increment;
      }

      CODI_INLINE Real getValue() const {
        double result;

        #pragma omp atomic read
        result = this->value;

        return result;
      }

    public:
      CODI_INLINE Atomic() : value() {}

      CODI_INLINE Atomic(const Atomic& other) {
        setValue(other.getValue());
      }

      CODI_INLINE Atomic(const Real& other) {
        setValue(other);
      }

      CODI_INLINE Atomic& operator = (const Atomic& other) {
        return operator = (other.getValue());
      }

      CODI_INLINE Atomic& operator = (const Real& other) {
        setValue(other);
        return *this;
      }

      CODI_INLINE Atomic& operator += (const Atomic& other) {
        return operator += (other.getValue());
      }

      CODI_INLINE Atomic& operator += (const Real& other) {
        addValue(other);
        return *this;
      }

      CODI_INLINE operator Real () const {
        return getValue();
      }

      CODI_INLINE bool isTotalZero() const {
        return codi::isTotalZero(getValue());
      }
  };

  template<typename _Tape>
  class Atomic<ActiveReal<_Tape>, enableIfForwardTape<_Tape>> : public ActiveReal<_Tape> {
    public:
      using Tape = _Tape;
      using Real = ActiveReal<Tape>;
      using NestedReal = typename Real::Real;
      using GradientValue = typename Real::GradientValue;

    private:
      CODI_INLINE void internalSetValue(const Real& newValue) {

        Atomic<NestedReal>* atomicPrimal = reinterpret_cast<Atomic<NestedReal>*>(&this->value());
        Atomic<GradientValue>* atomicGradient = reinterpret_cast<Atomic<GradientValue>*>(&this->gradient());

        *atomicPrimal = newValue.value();
        *atomicGradient = newValue.gradient();
      }

      CODI_INLINE void internalAddValue(const Real& increment) {

        Atomic<NestedReal>* atomicPrimal = reinterpret_cast<Atomic<NestedReal>*>(&this->value());
        Atomic<GradientValue>* atomicGradient = reinterpret_cast<Atomic<GradientValue>*>(&this->gradient());

        *atomicPrimal += increment.value();
        *atomicGradient += increment.gradient();
      }

      CODI_INLINE Real internalGetValue() const {
        Real result;

        const Atomic<NestedReal>* atomicPrimal = reinterpret_cast<const Atomic<NestedReal>*>(&this->value());
        const Atomic<GradientValue>* atomicGradient = reinterpret_cast<const Atomic<GradientValue>*>(&this->gradient());

        result.value() = *atomicPrimal;
        result.gradient() = *atomicGradient;

        return result;
      }

    public:
      CODI_INLINE Atomic() : Real() {}

      CODI_INLINE Atomic(const Atomic& other) : Real() {
        internalSetValue(other.internalGetValue());
      }

      CODI_INLINE Atomic(const Real& other) : Real() {
        internalSetValue(other);
      }

      CODI_INLINE Atomic& operator = (const Atomic& other) {
        return operator = (other.internalGetValue());
      }

      CODI_INLINE Atomic& operator = (const Real& other) {
        internalSetValue(other);
        return *this;
      }

      CODI_INLINE Atomic& operator += (const Atomic& other) {
        return operator += (other.internalGetValue());
      }

      CODI_INLINE Atomic& operator += (const Real& other) {
        internalAddValue(other);
        return *this;
      }

      CODI_INLINE operator Real () const {
        return internalGetValue();
      }

      CODI_INLINE bool isTotalZero() const {
        return codi::isTotalZero(internalGetValue());
      }
  };

  template<typename Real>
  struct InternalRemoveAtomic {
      using type = Real;
  };

  template<typename Real>
  struct InternalRemoveAtomic<Atomic<Real>> {
      using type = Real;
  };

  template<typename Real>
  using RemoveAtomic = typename InternalRemoveAtomic<Real>::type;

}
