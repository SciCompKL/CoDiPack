/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <type_traits>

#include "../../../expressions/activeType.hpp"
#include "../../../traits/atomicTraits.hpp"
#include "../../../traits/realTraits.hpp"
#include "../../../traits/tapeTraits.hpp"
#include "../atomicInterface.hpp"
#include "macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Atomic implementation for OpenMP.
   *
   * OpenMP atomics are disabled for all types by default. Atomics for arithmetic types and forward CoDiPack types are
   * enabled by specializations.
   *
   * See also AtomicInterface.
   *
   * @tparam T_Type    The underlying data type.
   * @tparam T_Sfinae  Additional SFNIAE parameter for enable-if constructs.
   */
  template<typename T_Type, typename T_Sfinae = void>
  struct OpenMPAtomicImpl : public AtomicInterface<T_Type, OpenMPAtomicImpl<T_Type, T_Sfinae>> {
    public:
      using Type = T_Type;  ///< See OpenMPAtomicImpl.

      OpenMPAtomicImpl() = delete;  ///< Constructor is deleted, will throw errors for unspecialized instantiations.
  };

#ifndef DOXYGEN_DISABLE

  // Specialization for arithmetic types.
  template<typename T_Type>
  struct OpenMPAtomicImpl<T_Type, typename std::enable_if<std::is_arithmetic<T_Type>::value>::type>
      : public AtomicInterface<
            T_Type, OpenMPAtomicImpl<T_Type, typename std::enable_if<std::is_arithmetic<T_Type>::value>::type>> {
    public:
      using Type = T_Type;
      using Base =
          AtomicInterface<T_Type,
                          OpenMPAtomicImpl<T_Type, typename std::enable_if<std::is_arithmetic<T_Type>::value>::type>>;

    private:
      Type value;

      CODI_INLINE void setValue(Type const& newValue) {
        CODI_OMP_ATOMIC(write)
        this->value = newValue;
      }

      CODI_INLINE Type getValue() const {
        Type result;
        CODI_OMP_ATOMIC(read)
        result = this->value;
        return result;
      }

    public:
      CODI_INLINE OpenMPAtomicImpl() : Base(), value() {}

      CODI_INLINE OpenMPAtomicImpl(OpenMPAtomicImpl const& other) : Base() {
        setValue(other.getValue());
      }

      CODI_INLINE OpenMPAtomicImpl(Type const& other) : Base() {
        setValue(other);
      }

      CODI_INLINE OpenMPAtomicImpl& operator=(OpenMPAtomicImpl const& other) {
        return operator=(other.getValue());
      }

      CODI_INLINE OpenMPAtomicImpl& operator=(Type const& other) {
        setValue(other);
        return *this;
      }

      CODI_INLINE Type operator+=(OpenMPAtomicImpl const& other) {
        return operator+=(other.getValue());
      }

      CODI_INLINE Type operator+=(Type const& other) {
        Type result;
        CODI_OMP_ATOMIC(capture) {
          this->value += other;
          result = this->value;
        }
        return result;
      }

      CODI_INLINE Type operator++() {
        Type result;
        CODI_OMP_ATOMIC(capture)
        result = ++this->value;
        return result;
      }

      CODI_INLINE Type operator++(int) {
        Type result;
        CODI_OMP_ATOMIC(capture)
        result = this->value++;
        return result;
      }

      CODI_INLINE Type operator--() {
        Type result;
        CODI_OMP_ATOMIC(capture)
        result = --this->value;
        return result;
      }

      CODI_INLINE Type operator--(int) {
        Type result;
        CODI_OMP_ATOMIC(capture)
        result = this->value--;
        return result;
      }

      CODI_INLINE operator Type() const {
        return getValue();
      }
  };

  // Specialization for forward CoDiPack types. Acts on value and gradient with individual atomic operations.
  template<typename T_Tape>
  struct OpenMPAtomicImpl<ActiveType<T_Tape>, TapeTraits::EnableIfForwardTape<T_Tape>>
      : public AtomicInterface<ActiveType<T_Tape>,
                               OpenMPAtomicImpl<ActiveType<T_Tape>, TapeTraits::EnableIfForwardTape<T_Tape>>>,
        public ActiveType<T_Tape> {
    public:
      using Tape = CODI_DD(T_Tape, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));
      using Base = AtomicInterface<ActiveType<T_Tape>,
                                   OpenMPAtomicImpl<ActiveType<T_Tape>, TapeTraits::EnableIfForwardTape<T_Tape>>>;
      using Type = ActiveType<Tape>;
      using Real = typename Type::Real;
      using Gradient = typename Type::Gradient;

    private:
      CODI_INLINE void atomicSetValue(Type const& newValue) {
        OpenMPAtomicImpl<Real>* atomicValue = reinterpret_cast<OpenMPAtomicImpl<Real>*>(&this->value());
        OpenMPAtomicImpl<Gradient>* atomicGradient = reinterpret_cast<OpenMPAtomicImpl<Gradient>*>(&this->gradient());

        *atomicValue = newValue.value();
        *atomicGradient = newValue.gradient();
      }

      CODI_INLINE Type atomicGetValue() const {
        Type result;

        OpenMPAtomicImpl<Real> const* atomicValue = reinterpret_cast<OpenMPAtomicImpl<Real> const*>(&this->value());
        OpenMPAtomicImpl<Gradient> const* atomicGradient =
            reinterpret_cast<OpenMPAtomicImpl<Gradient> const*>(&this->gradient());

        result.value() = *atomicValue;
        result.gradient() = *atomicGradient;

        return result;
      }

    public:
      CODI_INLINE OpenMPAtomicImpl() : Base(), Type() {}

      CODI_INLINE OpenMPAtomicImpl(OpenMPAtomicImpl const& other) : Base(), Type() {
        atomicSetValue(other.atomicGetValue());
      }

      CODI_INLINE OpenMPAtomicImpl(Type const& other) : Base(), Type() {
        atomicSetValue(other);
      }

      CODI_INLINE OpenMPAtomicImpl& operator=(OpenMPAtomicImpl const& other) {
        return operator=(other.atomicGetValue());
      }

      CODI_INLINE OpenMPAtomicImpl& operator=(Type const& other) {
        atomicSetValue(other);
        return *this;
      }

      CODI_INLINE OpenMPAtomicImpl& operator+=(OpenMPAtomicImpl const& other) {
        return operator+=(other.atomicGetValue());
      }

      CODI_INLINE OpenMPAtomicImpl& operator+=(Type const& other) {
        OpenMPAtomicImpl<Real>* atomicValue = reinterpret_cast<OpenMPAtomicImpl<Real>*>(&this->value());
        OpenMPAtomicImpl<Gradient>* atomicGradient = reinterpret_cast<OpenMPAtomicImpl<Gradient>*>(&this->gradient());

        *atomicValue += other.value();
        *atomicGradient += other.gradient();
        return *this;
      }

      CODI_INLINE operator Type() const {
        return atomicGetValue();
      }
  };

#endif

  /// Wrapper for atomics for OpenMP.
  /// @tparam Type  An arithmetic type or CoDiPack forward type.
  template<typename Type>
  using OpenMPAtomic = OpenMPAtomicImpl<Type>;

  /// Declare OpenMPAtomic to be atomic in terms of AtomicTraits.
  template<typename T_Type>
  struct AtomicTraits::IsAtomic<OpenMPAtomic<T_Type>> : std::true_type {};

#ifndef DOXYGEN_DISABLE
  // Specialize IsTotalZero for OpenMPAtomic on arithmetic types.
  template<typename T_Type>
  struct RealTraits::IsTotalZero<
      OpenMPAtomicImpl<T_Type, typename std::enable_if<std::is_arithmetic<T_Type>::value>::type>> {
    public:

      using Type =
          CODI_DD(CODI_T(OpenMPAtomicImpl<T_Type, typename std::enable_if<std::is_arithmetic<T_Type>::value>::type>),
                  OpenMPAtomic<double>);

      static CODI_INLINE bool isTotalZero(Type const& v) {
        return typename Type::Type() == v;
      }
  };
#endif
}
