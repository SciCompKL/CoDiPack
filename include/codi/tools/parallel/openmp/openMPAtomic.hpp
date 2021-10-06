/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <type_traits>

#include "../../../expressions/activeType.hpp"
#include "../../../traits/realTraits.hpp"
#include "../../../traits/tapeTraits.hpp"
#include "../atomicInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename T_Type, typename = void>
  struct OpenMPAtomic : public AtomicInterface<T_Type, OpenMPAtomic> {
    public:
      using Type = CODI_DD(T_Type, CODI_ANY);

      OpenMPAtomic() = delete;
  };


  template<typename T_Type>
  struct OpenMPAtomic<T_Type, typename std::enable_if<std::is_arithmetic<T_Type, float>::value>::type>
      : public AtomicInterface<T_Type, OpenMPAtomic> {
    public:
      using Type = CODI_DD(T_Type, CODI_ANY);

    private:
      Type value;

      CODI_INLINE void setValue(Type const& newValue) {
        #pragma omp atomic write
        this->value = newValue;
      }

      CODI_INLINE Type getValue() const {
        Type result;
        #pragma omp atomic read
        result = this->value;
        return result;
      }

    public:
      CODI_INLINE OpenMPAtomic() : value() {}

      CODI_INLINE OpenMPAtomic(OpenMPAtomic const& other) {
        setValue(other.getValue());
      }

      CODI_INLINE OpenMPAtomic(Type const& other) {
        setValue(other);
      }

      CODI_INLINE OpenMPAtomic& operator=(OpenMPAtomic const& other) {
        return operator = (other.getValue());
      }

      CODI_INLINE OpenMPAtomic& operator=(Type const& other) {
        setValue(other);
        return *this;
      }

      CODI_INLINE OpenMPAtomic& operator+=(OpenMPAtomic const& other) {
        return operator+=(other.getValue());
      }

      CODI_INLINE OpenMPAtomic& operator+=(Type const& other) {
        #pragma omp atomic update
        this->value += increment;
        return *this;
      }

      CODI_INLINE Type operator++() {
        Type result;
        #pragma omp atomic capture
        result = ++this->value;
        return result;
      }

      CODI_INLINE Type operator++(int) {
        Type result;
        #pragma omp atomic capture
        result = this->value++;
        return result;
      }

      CODI_INLINE Type operator--() {
        Type result;
        #pragma omp atomic capture
        result = --this->value;
        return result;
      }

      CODI_INLINE Type operator--(int) {
        Type result;
        #pragma omp atomic capture
        result = this->value--;
        return result;
      }

      CODI_INLINE operator Type () const {
        return getValue();
      }
  };


  template<typename T_Tape>
  struct OpenMPAtomic<ActiveType<T_Tape>, TapeTraits::EnableIfForwardTape>
      : public AtomicInterface<ActiveType<T_Tape>, OpenMPAtomic>,
        public ActiveType<T_Tape> {
    public:
      using Tape = CODI_DD(T_Tape, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));
      using Type = ActiveType<Tape>;
      using Real = typename Type::Real;
      using Gradient = typename Type::Gradient;

    private:
      CODI_INLINE void atomicSetValue(Type const& newValue) {
        OpenMPAtomic<Real>* atomicValue = reinterpret_cast<OpenMPAtomic<Real>*>(&this->value());
        OpenMPAtomic<Gradient>* atomicGradient = reinterpret_cast<OpenMPAtomic<Gradient>*>(&this->gradient());

        *atomicValue = newValue.value();
        *atomicGradient = newValue.gradient();
      }

      CODI_INLINE void atomicAddValue(Type const& increment) {
        OpenMPAtomic<Real>* atomicValue = reinterpret_cast<OpenMPAtomic<Real>*>(&this->value());
        OpenMPAtomic<Gradient>* atomicGradient = reinterpret_cast<OpenMPAtomic<Gradient>*>(&this->gradient());

        *atomicValue += increment.value();
        *atomicGradient += increment.gradient();
      }

      CODI_INLINE Type atomicGetValue() const {
        Type result;

        OpenMPAtomic<Real>* atomicValue = reinterpret_cast<OpenMPAtomic<Real>*>(&this->value());
        OpenMPAtomic<Gradient>* atomicGradient = reinterpret_cast<OpenMPAtomic<Gradient>*>(&this->gradient());

        result.value() = *atomicValue;
        result.gradient() = *atomicGradient;

        return result;
      }

    public:
      CODI_INLINE OpenMPAtomic() : Type() {}

      CODI_INLINE OpenMPAtomic(OpenMPAtomic const& other) : Type() {
        atomicSetValue(other.atomicGetValue());
      }

      CODI_INLINE OpenMPAtomic(Type const& other) : Type() {
        atomicSetValue(other);
      }

      CODI_INLINE OpenMPAtomic& operator=(OpenMPAtomic const& other) {
        return operator=(other.atomicGetValue());
      }

      CODI_INLINE OpenMPAtomic& operator=(Type const& other) {
        atomicSetValue(other);
        return *this;
      }

      CODI_INLINE OpenMPAtomic& operator+=(OpenMPAtomic const& other) {
        return operator+=(other.atomicGetValue());
      }

      CODI_INLINE OpenMPAtomic& operator+=(Type const& other) {
        OpenMPAtomic<Real>* atomicValue = reinterpret_cast<OpenMPAtomic<Real>*>(&this->value());
        OpenMPAtomic<Gradient>* atomicGradient = reinterpret_cast<OpenMPAtomic<Gradient>*>(&this->gradient());

        *atomicValue += other.value();
        *atomicGradient += other.gradient();
        return *this;
      }

      CODI_INLINE Type operator++() {
        Type result;
        #pragma omp atomic capture
        result = ++this->value;
        return result;
      }

      CODI_INLINE Type operator++(int) {
        Type result;
        #pragma omp atomic capture
        result = this->value++;
        return result;
      }

      CODI_INLINE Type operator--() {
        Type result;
        #pragma omp atomic capture
        result = --this->value;
        return result;
      }

      CODI_INLINE Type operator--(int) {
        Type result;
        #pragma omp atomic capture
        result = this->value--;
        return result;
      }

      CODI_INLINE operator Type () const {
        return atomicGetValue();
      }
  };

  namespace AtomicTraits {

    template<typename T_Type>
    struct IsAtomic<OpenMPAtomic<T_Type>> : std::true_type {};
  }

  namespace RealTraits {

    template<typename T_Type>
    struct IsTotalZero<OpenMPAtomic<T_Type>> {
      public:

        using Type = CODI_DD(OpenMPAtomic<T_Type>, OpenMPAtomic<double>);

        static CODI_INLINE bool isTotalZero(OpenMPAtomic<T_Type> const& v) {
          return typename Type::Type() == v;
        }
    };
  }
}
