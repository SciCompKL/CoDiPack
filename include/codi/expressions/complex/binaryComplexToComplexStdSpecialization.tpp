/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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

// clang-format off
/*
 * This file defines function overloads of the kind:
 *
 * namespace std {
 *   codi::ComputeExpression<complex<Real>, ...> FUNCTION(complex<codi::ActiveType> const&, complex<codi::ActiveType> const&);
 *   codi::ComputeExpression<complex<Real>, ...> FUNCTION(complex<codi::ActiveType> const&, codi::ActiveType const&);
 *   codi::ComputeExpression<complex<Real>, ...> FUNCTION(complex<codi::ActiveType> const&, PassiveReal const&);
 *   codi::ComputeExpression<complex<Real>, ...> FUNCTION(codi::ActiveType const&,          complex<codi::ActiveType> const&);
 *   codi::ComputeExpression<complex<Real>, ...> FUNCTION(PassiveReal const&, complex<codi::ActiveType> const&);
 * }
 *
 * In order to include this file the user has to define the preprocessor macros OPERATION_LOGIC and FUNCTION.
 * OPERATION_LOGIC contains the name of the operation logic class. FUNCTION represents the normal name of that function
 * e.g. 'operator -' or 'sin'.
 *
 * The defines OPERATION_LOGIC and FUNCTION will be undefined at the end of this template.
 *
 * Prior to including this file, the user has to implement the operation's primal and derivative logic according to
 * BinaryOpInterface.
 */
// clang-format on

#ifndef OPERATION_LOGIC
  #error Please define a name for the binary expression.
#endif
#ifndef FUNCTION
  #error Please define the primal function representation.
#endif

// Create a correct include environment for viewing and programming in an IDE.
#ifndef FUNCTION
  #include <complex>

  #include "../../config.h"
  #include "../../misc/macros.hpp"
  #include "../activeType.hpp"
  #include "../computeExpression.hpp"
  #include "../constantExpression.hpp"
  #include "realToComplexCast.hpp"

  #define OPERATION_LOGIC codi::BinaryOperation
  #define FUNCTION func

namespace std {
#endif

  // Define (complex, complex) bindings

  /// Function overload for FUNCTION(complex, complex).
  template<typename Tape>
  CODI_INLINE auto FUNCTION(complex<codi::ActiveType<Tape>> const& argA, complex<codi::ActiveType<Tape>> const& argB) {
    return codi::ComputeExpression<complex<typename Tape::Real>, OPERATION_LOGIC, complex<codi::ActiveType<Tape>>,
                                   complex<codi::ActiveType<Tape>>>(argA, argB);
  }

  // Define (complex, real) bindings

  /// Function overload for FUNCTION(complex, real).
  template<typename Tape>
  CODI_INLINE auto FUNCTION(complex<codi::ActiveType<Tape>> const& argA, codi::ActiveType<Tape> const& argB) {
    return codi::ComputeExpression<complex<typename Tape::Real>, OPERATION_LOGIC, complex<codi::ActiveType<Tape>>,
                                   codi::RealToComplexCast<typename Tape::Real, codi::ActiveType<Tape>>>(
        argA, codi::RealToComplexCast<typename Tape::Real, codi::ActiveType<Tape>>(argB));
  }

  /// Function overload for FUNCTION(complex, passive real).
  template<typename Tape>
  CODI_INLINE auto FUNCTION(complex<codi::ActiveType<Tape>> const& argA, typename Tape::PassiveReal const& argB) {
    return codi::ComputeExpression<complex<typename Tape::Real>, OPERATION_LOGIC, complex<codi::ActiveType<Tape>>,
                                   codi::ConstantExpression<typename Tape::PassiveReal>>(
        argA, codi::ConstantExpression<typename Tape::PassiveReal>(argB));
  }

  // Define (real, complex) bindings

  /// Function overload for FUNCTION(real, complex).
  template<typename Tape>
  CODI_INLINE auto FUNCTION(codi::ActiveType<Tape> const& argA, complex<codi::ActiveType<Tape>> const& argB) {
    return codi::ComputeExpression<complex<typename Tape::Real>, OPERATION_LOGIC,
                                   codi::RealToComplexCast<typename Tape::Real, codi::ActiveType<Tape>>,
                                   complex<codi::ActiveType<Tape>>>(
        codi::RealToComplexCast<typename Tape::Real, codi::ActiveType<Tape>>(argA), argB);
  }

  /// Function overload for FUNCTION(passive real, complex).
  template<typename Tape>
  CODI_INLINE auto FUNCTION(typename Tape::PassiveReal const& argA, complex<codi::ActiveType<Tape>> const& argB) {
    return codi::ComputeExpression<complex<typename Tape::Real>, OPERATION_LOGIC,
                                   codi::ConstantExpression<typename Tape::PassiveReal>,
                                   complex<codi::ActiveType<Tape>>>(
        codi::ConstantExpression<typename Tape::PassiveReal>(argA), argB);
  }

// Create a correct include environment for viewing and programming in an IDE.
#ifndef FUNCTION
}
#endif

#undef FUNCTION
#undef OPERATION_LOGIC
