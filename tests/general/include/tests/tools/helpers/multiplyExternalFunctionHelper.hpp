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

#include <codi.hpp>

template<typename T_Number, typename = void>
struct MultiplyExternalFunctionHelper {
  public:
    using Number = CODI_DECLARE_DEFAULT(T_Number, codi::ActiveType<CODI_ANY>);

    static Number create(Number const& x1, Number const& x2, bool passive) {
      codi::CODI_UNUSED(passive);

      return x1 * x2;
    }
};

template<typename T_Number>
struct MultiplyExternalFunctionHelper<T_Number, codi::TapeTraits::EnableIfReverseTape<typename T_Number::Tape>> {
  public:
    using Number = CODI_DECLARE_DEFAULT(T_Number, codi::ActiveType<CODI_ANY>);

    using Real = typename Number::Real;

    static void funcCall(Number& w, Number const& x1, Number const& x2) {
      w = x1 * x2;
    }

    static void func_primal(Real const* x, size_t m, Real* y, size_t n, codi::ExternalFunctionUserData* d) {
      codi::CODI_UNUSED(m, n, d);

      y[0] = x[0] * x[1];
    }

    static void func_reverse(Real const* x, Real* x_b, size_t m, Real const* y, Real const* y_b, size_t n,
                             codi::ExternalFunctionUserData* d) {
      codi::CODI_UNUSED(m, n, y, d);

      x_b[0] = x[1] * y_b[0];
      x_b[1] = x[0] * y_b[0];
    }

    static void func_forward(Real const* x, Real const* x_d, size_t m, Real* y, Real* y_d, size_t n,
                             codi::ExternalFunctionUserData* d) {
      codi::CODI_UNUSED(m, n, y, d);

      y[0] = x[0] * x[1];
      y_d[0] = x[1] * x_d[0] + x_d[1] * x[0];
    }

    static Number create(Number const& x1, Number const& x2, bool primalFuncUsesADType) {
      static codi::ExternalFunctionHelper<Number> eh(primalFuncUsesADType);

      Number w;

      eh.addInput(x1);
      eh.addInput(x2);
      eh.addOutput(w);

      if (primalFuncUsesADType) {
        eh.callPrimalFuncWithADType(MultiplyExternalFunctionHelper::funcCall, w, x1, x2);
      } else {
        eh.callPrimalFunc(MultiplyExternalFunctionHelper::func_primal);
      }

      eh.addToTape(MultiplyExternalFunctionHelper::func_reverse, MultiplyExternalFunctionHelper::func_forward,
                   MultiplyExternalFunctionHelper::func_primal);

      return w;
    }
};
