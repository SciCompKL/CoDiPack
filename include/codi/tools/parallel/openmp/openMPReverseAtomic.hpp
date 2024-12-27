/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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
#include "../reverseAtomicInterface.hpp"
#include "macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Reverse atomic implementation for OpenMP.
   *
   * OpenMP reverse atomics are disabled for all types by default. Reverse atomics for arithmetic types and forward
   * CoDiPack types are enabled by specializations.
   *
   * See also ReverseAtomicInterface.
   *
   * @tparam T_Type    The underlying data type.
   * @tparam T_Sfinae  Additional SFNIAE parameter for enable-if constructs.
   */
  template<typename T_Type, typename T_Sfinae = void>
  struct OpenMPReverseAtomicImpl : public ReverseAtomicInterface<T_Type, OpenMPReverseAtomicImpl<T_Type, T_Sfinae>> {
    public:
      using Type = T_Type;  ///< See OpenMPReverseAtomicImpl.

      OpenMPReverseAtomicImpl() = delete;  ///< Constructor is deleted, will throw errors for unspecialized
                                           ///< instantiations.
  };

#ifndef DOXYGEN_DISABLE

  // Specialization for arithmetic types.
  template<typename T_Type>
  struct OpenMPReverseAtomicImpl<T_Type, typename std::enable_if<std::is_arithmetic<T_Type>::value>::type>
      : public ReverseAtomicInterface<
            T_Type, OpenMPReverseAtomicImpl<T_Type, typename std::enable_if<std::is_arithmetic<T_Type>::value>::type>> {
    public:
      using Type = T_Type;
      using Base = ReverseAtomicInterface<
          T_Type, OpenMPReverseAtomicImpl<T_Type, typename std::enable_if<std::is_arithmetic<T_Type>::value>::type>>;

    private:
      Type value;

    public:
      CODI_INLINE OpenMPReverseAtomicImpl() : Base(), value() {}

      CODI_INLINE OpenMPReverseAtomicImpl(OpenMPReverseAtomicImpl const& other) : Base(), value(other.value) {}

      CODI_INLINE OpenMPReverseAtomicImpl(Type const& other) : Base(), value(other) {}

      CODI_INLINE OpenMPReverseAtomicImpl& operator=(OpenMPReverseAtomicImpl const& other) {
        return operator=(other.value);
      }

      CODI_INLINE OpenMPReverseAtomicImpl& operator=(Type const& other) {
        this->value = other;
        return *this;
      }

      CODI_INLINE void operator+=(OpenMPReverseAtomicImpl const& other) {
        operator+=(other.value);
      }

      CODI_INLINE void operator+=(Type const& other) {
        CODI_OMP_ATOMIC(update)
        this->value += other;
      }

      CODI_INLINE operator Type() const {
        return value;
      }
  };

  // Specialization for forward CoDiPack types. Acts on value and gradient with individual atomic operations.
  template<typename T_Type>
  struct OpenMPReverseAtomicImpl<T_Type, TapeTraits::EnableIfForwardTape<typename T_Type::Tape>>
      : public ReverseAtomicInterface<
            T_Type, OpenMPReverseAtomicImpl<T_Type, TapeTraits::EnableIfForwardTape<typename T_Type::Tape>>>,
        public T_Type {
    public:
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);
      using Base = ReverseAtomicInterface<
          T_Type, OpenMPReverseAtomicImpl<T_Type, TapeTraits::EnableIfForwardTape<typename T_Type::Tape>>>;
      using Tape = typename Type::Tape;
      using Real = typename Type::Real;
      using Gradient = typename Type::Gradient;

      CODI_INLINE OpenMPReverseAtomicImpl() : Base(), Type() {}

      CODI_INLINE OpenMPReverseAtomicImpl(OpenMPReverseAtomicImpl const& other) : Base(), Type(other) {}

      CODI_INLINE OpenMPReverseAtomicImpl(Type const& other) : Base(), Type(other) {}

      CODI_INLINE OpenMPReverseAtomicImpl& operator=(OpenMPReverseAtomicImpl const& other) {
        return operator=(static_cast<Type const&>(other));
      }

      CODI_INLINE OpenMPReverseAtomicImpl& operator=(Type const& other) {
        Type::operator=(other);
        return *this;
      }

      CODI_INLINE void operator+=(OpenMPReverseAtomicImpl const& other) {
        operator+=(static_cast<Type const&>(other));
      }

      CODI_INLINE void operator+=(Type const& other) {
        OpenMPReverseAtomicImpl<Real>* atomicValue = reinterpret_cast<OpenMPReverseAtomicImpl<Real>*>(&this->value());
        OpenMPReverseAtomicImpl<Gradient>* atomicGradient =
            reinterpret_cast<OpenMPReverseAtomicImpl<Gradient>*>(&this->gradient());

        *atomicValue += other.value();
        *atomicGradient += other.gradient();
      }

      CODI_INLINE operator Type() const {
        return static_cast<Type>(*this);
      }
  };

#endif

  /// Wrapper for reverse atomics for OpenMP.
  /// @tparam Type  An arithmetic type or CoDiPack forward type.
  template<typename Type>
  using OpenMPReverseAtomic = OpenMPReverseAtomicImpl<Type>;

  /// Declare OpenMPAtomic to be atomic in terms of AtomicTraits.
  template<typename T_Type>
  struct AtomicTraits::IsAtomic<OpenMPReverseAtomic<T_Type>> : std::true_type {};

#ifndef DOXYGEN_DISABLE
  // Specialize IsTotalZero for OpenMPReverseAtomic on arithmetic types.
  template<typename T_Type>
  struct RealTraits::IsTotalZero<
      OpenMPReverseAtomicImpl<T_Type, typename std::enable_if<std::is_arithmetic<T_Type>::value>::type>> {
    public:

      using Type = CODI_DD(
          CODI_T(OpenMPReverseAtomicImpl<T_Type, typename std::enable_if<std::is_arithmetic<T_Type>::value>::type>),
          OpenMPReverseAtomic<double>);

      static CODI_INLINE bool isTotalZero(Type const& v) {
        return typename Type::Type() == v;
      }
  };
#endif
}
