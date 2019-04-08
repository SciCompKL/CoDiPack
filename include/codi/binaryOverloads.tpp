/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 *          Prof. Robin Hogan, (Univ. of Reading).
 *
 * Originally based on Adept 1.0 (http://www.met.rdg.ac.uk/clouds/adept/)
 * released under GPL 3.0 (Copyright (C) 2012-2013 Robin Hogan and the University of Reading).
 */

/*
 * In order to include this file the user has to define the preprocessor macros OPERATION_LOGIC and FUNCTION.
 * OPERATION_LOGIC contains the name of the operation logic class. FUNCTION represents the normal name of that function
 * e.g. 'operator -' or 'sin'.
 *
 * The defines OPERATION_LOGIC and FUNCTION will be undefined at the end of this template.
 *
 * Prior to including this file, the user has to implement the operation's primal and derivative logic according to BinaryOpInterface.
 */

#ifndef OPERATION_LOGIC
  #error Please define a name for the binary expression.
#endif
#ifndef FUNCTION
  #error Please define the primal function representation.
#endif

#include "macros.h"

/**
 * @brief Overload for FUNCTION with the CoDiPack expressions.
 *
 * @param[in] a  The first argument of the operation.
 * @param[in] b  The second argument of the operation.
 *
 * @return BinaryOp11 instanciated for OPERATION_LOGIC.
 *
 * @tparam Real  The real type used in the active types.
 * @tparam    A  The expression for the first argument of the function
 * @tparam    B  The expression for the second argument of the function
 */
template <typename Real, class A, class B>
CODI_INLINE BinaryOp11<Real, A, B, OPERATION_LOGIC> FUNCTION(const Expression<Real, A>& a, const Expression<Real, B>& b) {
  return BinaryOp11<Real, A, B, OPERATION_LOGIC>(a.cast(), b.cast());
}
/**
 * @brief Overload for FUNCTION with the CoDiPack expressions.
 *
 * @param[in] a  The first argument of the operation.
 * @param[in] b  The second argument of the operation.
 *
 * @return BinaryOp10 instanciated for OPERATION_LOGIC.
 *
 * @tparam Real  The real type used in the active types.
 * @tparam    A  The expression for the first argument of the function
 */
template <typename Real, class A>
CODI_INLINE BinaryOp10<Real, A, OPERATION_LOGIC> FUNCTION(const Expression<Real, A>& a, const typename TypeTraits<Real>::PassiveReal& b) {
  return BinaryOp10<Real, A, OPERATION_LOGIC>(a.cast(), b);
}
/**
 * @brief Overload for FUNCTION with the CoDiPack expressions.
 *
 * @param[in] a  The first argument of the operation.
 * @param[in] b  The second argument of the operation.
 *
 * @return BinaryOp01 instanciated for OPERATION_LOGIC.
 *
 * @tparam Real  The real type used in the active types.
 * @tparam    B  The expression for the second argument of the function
 */
template <typename Real, class B>
CODI_INLINE BinaryOp01<Real, B, OPERATION_LOGIC> FUNCTION(const typename TypeTraits<Real>::PassiveReal& a, const Expression<Real, B>& b) {
  return BinaryOp01<Real, B, OPERATION_LOGIC>(a, b.cast());
}

#undef FUNCTION
#undef OPERATION_LOGIC
