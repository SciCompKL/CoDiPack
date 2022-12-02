/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

#include <complex>

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../../misc/self.hpp"
#include "../../traits/computationTraits.hpp"
#include "../activeType.hpp"
#include "../aggregate/aggregatedActiveType.hpp"
#include "../binaryExpression.hpp"
#include "../expressionInterface.hpp"
#include "../expressionMemberOperations.hpp"
#include "../real/binaryOperators.hpp"
#include "../real/unaryOperators.hpp"
#include "../static/staticContextAggregateActiveType.hpp"
#include "../unaryExpression.hpp"
#include "allOperators.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   *  Define member operations for complex type expressions.
   *
   *  @tparam T_Real See ExpressionMemberOperations.
   *  @tparam T_Impl See ExpressionMemberOperations.
   */
  template<typename T_Real, typename T_Impl>
  struct ExpressionMemberOperations<
      T_Real, T_Impl,
      typename std::enable_if<std::is_same<T_Real, std::complex<typename T_Real::value_type>>::value>::type> {
      using Real = CODI_DD(T_Real, std::complex<double>);                     ///< See ExpressionMemberOperations.
      using Impl = CODI_DD(T_Impl, CODI_T(ExpressionInterface<Real, void>));  ///< See ExpressionMemberOperations.

      using InnerType = typename Real::value_type;  ///< Inner type of the complex.

      /// Expression returned by real.
      using ExpressionComplexReal = UnaryExpression<InnerType, Impl, OperationComplexReal>;

      /// real member function for complex.
      ExpressionComplexReal real() const {
        return ExpressionComplexReal(cast());
      }

      /// Expression returned by imag.
      using ExpressionComplexImag = UnaryExpression<InnerType, Impl, OperationComplexImag>;

      /// imag member function for complex.
      ExpressionComplexImag imag() const {
        return ExpressionComplexImag(cast());
      }

    private:
      CODI_INLINE Impl const& cast() const {
        return static_cast<Impl const&>(*this);
      }
  };

  namespace RealTraits {

    /// Specialize real traits for std::complex.
    template<typename T_InnerReal>
    struct AggregatedTypeTraits<std::complex<T_InnerReal>>
        : public ArrayAggregatedTypeTraitsBase<std::complex<T_InnerReal>, T_InnerReal,
                                               std::complex<RealTraits::Real<T_InnerReal>>, 2> {};
  }

  /// Implementation of AggregatedActiveType for complex types.
  /// @tparam T_InnerActive The active type which is used in the inner store.
  /// @tparam T_Impl The implementation if it extends this class, otherwise nothing.
  template<typename T_InnerActive, typename T_Impl = Self>
  struct ActiveComplex : public AggregatedActiveType<std::complex<typename T_InnerActive::Real>, T_InnerActive,
                                                     ReturnSelf<T_Impl, ActiveComplex<T_InnerActive, T_Impl>>>,
                         public AssignmentOperators<typename T_InnerActive::Real, false,
                                                    ReturnSelf<T_Impl, ActiveComplex<T_InnerActive, T_Impl>>>,
                         public AssignmentOperators<std::complex<typename T_InnerActive::Real>, false,
                                                    ReturnSelf<T_Impl, ActiveComplex<T_InnerActive, T_Impl>>> {
    public:

      using InnerActive = CODI_DD(T_InnerActive, CODI_T(ActiveType<CODI_ANY>));        ///< See ActiveComplex.
      using Impl = CODI_DD(CODI_T(ReturnSelf<T_Impl, ActiveComplex>), ActiveComplex);  ///< See ActiveComplex.

      using InnerReal = typename InnerActive::Real;                  ///< Inner real value.
      using Real = std::complex<InnerReal>;                          ///< Complex value with the inner real.
      using PassiverInnerReal = RealTraits::PassiveReal<InnerReal>;  ///< Passive inner value.

      using Base = AggregatedActiveType<Real, InnerActive, Impl>;  ///< Abbreviation for the base class.
      using AssignReal = AssignmentOperators<InnerReal, false, ReturnSelf<Impl, ActiveComplex>>;  ///< Abbreviation for
                                                                                                  ///< real assignments.
      using AssignComplex =
          AssignmentOperators<std::complex<InnerReal>, false, ReturnSelf<Impl, ActiveComplex>>;  ///< Abbreviations for
                                                                                                 ///< complex
                                                                                                 ///< assignments.

      using value_type = InnerActive;  ///< std::complex interface.

      using Base::Base;       ///< Use Base constructors.
      using Base::operator=;  ///< Use Base assign operators.

      using AssignReal::operator+=;     ///< Use real operator +=
      using AssignComplex::operator+=;  ///< Use real operator +=
      using AssignReal::operator-=;     ///< Use real operator -=
      using AssignComplex::operator-=;  ///< Use real operator -=
      using AssignReal::operator*=;     ///< Use real operator *=
      using AssignComplex::operator*=;  ///< Use real operator *=
      using AssignReal::operator/=;     ///< Use real operator /=
      using AssignComplex::operator/=;  ///< Use real operator /=

      /// Constructor.
      template<typename ArgR>
      CODI_INLINE ActiveComplex(ExpressionInterface<InnerReal, ArgR> const& argR) : Base() {
        Base::arrayValue[0] = argR;
      }

      /// Constructor.
      CODI_INLINE ActiveComplex(PassiverInnerReal const& argR) : Base() {
        Base::arrayValue[0] = argR;
      }

      /// Constructor.
      template<typename ArgR, typename ArgI>
      CODI_INLINE ActiveComplex(ExpressionInterface<InnerReal, ArgR> const& argR,
                                ExpressionInterface<InnerReal, ArgI> const& argI)
          : Base() {
        Base::arrayValue[0] = argR;
        Base::arrayValue[1] = argI;
      }

      /// Constructor.
      template<typename ArgI>
      CODI_INLINE ActiveComplex(PassiverInnerReal const& argR, ExpressionInterface<InnerReal, ArgI> const& argI)
          : Base() {
        Base::arrayValue[0] = argR;
        Base::arrayValue[1] = argI;
      }

      /// Constructor.
      template<typename ArgR>
      CODI_INLINE ActiveComplex(ExpressionInterface<InnerReal, ArgR> const& argR, PassiverInnerReal const& argI)
          : Base() {
        Base::arrayValue[0] = argR;
        Base::arrayValue[1] = argI;
      }

      /// Constructor.
      CODI_INLINE ActiveComplex(PassiverInnerReal const& argR, PassiverInnerReal const& argI) : Base() {
        Base::arrayValue[0] = argR;
        Base::arrayValue[1] = argI;
      }

      CODI_INLINE ~ActiveComplex() = default;  ///< Destructor
  };

  /// Specialization of ActiveResultImpl for std::complex.
  template<typename T_InnerReal, typename T_Tape>
  struct ExpressionTraits::ActiveResultImpl<std::complex<T_InnerReal>, T_Tape, false> {
      using InnerReal = CODI_DD(T_InnerReal, CODI_ANY);  ///< See ActiveResultImpl.
      using Tape = CODI_DD(T_Tape, CODI_ANY);            ///< See ActiveResultImpl.

      /// Active result of the inner type.
      using InnerActiveResult = ExpressionTraits::ActiveResult<InnerReal, Tape, false>;

      /// The resulting active type of an expression.
#if CODI_SpecializeStdComplex
      using ActiveResult = std::complex<InnerActiveResult>;
#else
      using ActiveResult = ActiveComplex<InnerActiveResult>;
#endif
  };

  /// Specialization of ActiveResultImpl for std::complex in a static context.
  template<typename T_InnerReal, typename T_Tape>
  struct ExpressionTraits::ActiveResultImpl<std::complex<T_InnerReal>, T_Tape, true> {
      using InnerReal = CODI_DD(T_InnerReal, CODI_ANY);  ///< See ActiveResultImpl.
      using Tape = CODI_DD(T_Tape, CODI_ANY);            ///< See ActiveResultImpl.

      /// Active result of the inner type.
      using InnerActiveResult = ExpressionTraits::ActiveResult<InnerReal, Tape, true>;

      /// The resulting active type of an expression.
      using ActiveResult = StaticAggregatedActiveType<std::complex<InnerReal>, InnerActiveResult>;
  };
}

#if CODI_SpecializeStdComplex
namespace std {

  /*******************************************************************************/
  /// @name Specialization of std::complex for active types.
  /// @{

  /// Specialization of std::complex for CoDiPack types.
  /// @tparam T_Tape  The tape that manages all expressions created with this type.
  template<typename T_Tape>
  class complex<codi::ActiveType<T_Tape>>
      : public codi::ActiveComplex<codi::ActiveType<T_Tape>, complex<codi::ActiveType<T_Tape>>> {
    public:

      /// See complex.
      using Tape = CODI_DD(T_Tape, CODI_T(codi::FullTapeInterface<double, double, int, codi::EmptyPosition>));

      using Real = complex<typename Tape::Real>;     ///< See ActiveComplex::Real.
      using InnerActive = codi::ActiveType<Tape>;    ///< See ActiveComplex::InnerActive.
      using InnerReal = typename InnerActive::Real;  ///< See ActiveComplex::InnerReal.

      /// Base class abbreviation.
      using Base = codi::ActiveComplex<codi::ActiveType<Tape>, complex<codi::ActiveType<T_Tape>>>;

      using Base::Base;                                           ///< Use Base constructors.
      CODI_INLINE complex(complex const& value) : Base(value) {}  ///< Constructor.

      CODI_INLINE ~complex() = default;  ///< Destructor

      using Base::operator=;  ///< Use Base assign operators.

      /// Assign operator.
      CODI_INLINE complex& operator=(complex const& value) {
        Base::store(value);

        return *this;
      }
  };

  /// @}
  /*******************************************************************************/
  /// @name Binary  operators and standard math library binary functions
  /// @{

  #define FUNCTION operator+
  #define OPERATION_LOGIC codi::OperationAdd
  #include "binaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION operator-
  #define OPERATION_LOGIC codi::OperationSubstract
  #include "binaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION operator*
  #define OPERATION_LOGIC codi::OperationMultiply
  #include "binaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION operator/
  #define OPERATION_LOGIC codi::OperationDivide
  #include "binaryComplexToComplexStdSpecialization.tpp"

  // polar needs to be created in the codi namespace.

  #define FUNCTION pow
  #define OPERATION_LOGIC codi::OperationPow
  #include "binaryComplexToComplexStdSpecialization.tpp"

  /// @}
  /*******************************************************************************/
  /// @name Unary operators and standard math library unary functions
  /// @{

  /// Function overload for operator +
  template<typename Tape>
  CODI_INLINE complex<codi::ActiveType<Tape>> const& operator+(complex<codi::ActiveType<Tape>> const& arg) {
    return arg;
  }

  #define FUNCTION operator-
  #define OPERATION_LOGIC codi::OperationUnaryMinus
  #include "unaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION abs
  #define OPERATION_LOGIC codi::OperationComplexAbs
  #include "unaryComplexToRealStdSpecialization.tpp"

  #define FUNCTION acos
  #define OPERATION_LOGIC codi::OperationAcos
  #include "unaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION acosh
  #define OPERATION_LOGIC codi::OperationAcosh
  #include "unaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION arg
  #define OPERATION_LOGIC codi::OperationComplexArg
  #include "unaryComplexToRealStdSpecialization.tpp"

  #define FUNCTION asin
  #define OPERATION_LOGIC codi::OperationAsin
  #include "unaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION asinh
  #define OPERATION_LOGIC codi::OperationAsinh
  #include "unaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION atan
  #define OPERATION_LOGIC codi::OperationAtan
  #include "unaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION atanh
  #define OPERATION_LOGIC codi::OperationAtanh
  #include "unaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION conj
  #define OPERATION_LOGIC codi::OperationComplexConj
  #include "unaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION cos
  #define OPERATION_LOGIC codi::OperationCos
  #include "unaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION cosh
  #define OPERATION_LOGIC codi::OperationCosh
  #include "unaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION exp
  #define OPERATION_LOGIC codi::OperationExp
  #include "unaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION imag
  #define OPERATION_LOGIC codi::OperationComplexImag
  #include "unaryComplexToRealStdSpecialization.tpp"

  #define FUNCTION log
  #define OPERATION_LOGIC codi::OperationLog
  #include "unaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION log10
  #define OPERATION_LOGIC codi::OperationLog10
  #include "unaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION norm
  #define OPERATION_LOGIC codi::OperationComplexNorm
  #include "unaryComplexToRealStdSpecialization.tpp"

  #define FUNCTION proj
  #define OPERATION_LOGIC codi::OperationComplexProj
  #include "unaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION real
  #define OPERATION_LOGIC codi::OperationComplexReal
  #include "unaryComplexToRealStdSpecialization.tpp"

  #define FUNCTION sin
  #define OPERATION_LOGIC codi::OperationSin
  #include "unaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION sinh
  #define OPERATION_LOGIC codi::OperationSinh
  #include "unaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION sqrt
  #define OPERATION_LOGIC codi::OperationSqrt
  #include "unaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION tan
  #define OPERATION_LOGIC codi::OperationTan
  #include "unaryComplexToComplexStdSpecialization.tpp"

  #define FUNCTION tanh
  #define OPERATION_LOGIC codi::OperationTanh
  #include "unaryComplexToComplexStdSpecialization.tpp"

  /// @}
}
#endif
