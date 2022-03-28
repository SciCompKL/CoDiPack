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

#include <vector>

#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../misc/macros.hpp"
#include "../../tapes/interfaces/fullTapeInterface.hpp"
#include "../../tapes/misc/vectorAccessInterface.hpp"
#include "../../traits/tapeTraits.hpp"
#include "../data/externalFunctionUserData.hpp"
#include "externalFunctionHelper.hpp"

extern int enzyme_dup;
extern int enzyme_out;
extern int enzyme_const;

template<typename... Args>
void __enzyme_autodiff(void*, Args...);

template<typename... Args>
void __enzyme_fwddiff(void*, Args...);

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Helper class for the implementation of an external function with Enzyme in CoDiPack.
   *
   * The class helps the user to create derivative functions with Enzyme and add them to the tape. See
   * ExternalFunctionHelper for the general configuration options and procedures. This class only supports the operation
   * mode 1 "Implemented primal function". In this implementation, this mode requires just one call to
   * callAndAddToTape() which calls the primal function and adds it to the tape.
   *
   * The procedure of pushing a differentiated function with Enzyme is as follows.
   * 1. All function inputs and outputs are specified.
   * 2. Call callAndAddToTape method with the primal function.
   *
   * An example is:
   * \snippet examples/Example_24_Enzyme_external_function_helper.cpp Enzyme-differentiated function
   *
   * This helper provides an additional overload where the inputs and outputs can be directly added:
   * \snippet examples/Example_24_Enzyme_external_function_helper.cpp Enzyme-differentiated function - short
   *
   * The full example can be found at \ref Example_24_Enzyme_external_function_helper.
   *
   * @tparam T_Type  The CoDiPack type that is used outside of the external function.
   */
  template<typename T_Type>
  struct EnzymeExternalFunctionHelper : public ExternalFunctionHelper<T_Type> {
    public:

      /// See ExternalFunctionHelper.
      using Type = CODI_DD(T_Type, CODI_T(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

      using Base = ExternalFunctionHelper<Type>;  ///< Base class abbreviation.

      using Real = typename Type::Real;              ///< See LhsExpressionInterface.
      using Identifier = typename Type::Identifier;  ///< See LhsExpressionInterface.

      /// See LhsExpressionInterface.
      using Tape = CODI_DD(typename Type::Tape, CODI_T(FullTapeInterface<double, double, int, CODI_ANY>));

      using PrimalFunc = typename Base::PrimalFunc;  ///< See ExternalFunctionHelper.

    public:

      /// Constructor
      EnzymeExternalFunctionHelper() : Base(false) {}

    private:
      // Hide base class functionality.
      using Base::addToTape;
      using Base::callPrimalFunc;
      using Base::callPrimalFuncWithADType;

    public:
      // Implement enzyme specific functionality.

      /// Calls the primal function and adds the Enzyme generated functions to the tape.
      ///
      /// Is a combination of ExternalFunctionHelper::callPrimalFunc and ExternalFunctionHelper::addToTape.
      template<PrimalFunc func>
      void callAndAddToTape() {
        Base::callPrimalFunc(func);
        Base::addToTape(enzymeDiff_b<func>, enzymeDiff_d<func>, func);
      }

      /// Adds all inputs in x and outputs in y to the external function and then calls callAndAddToTape().
      template<PrimalFunc func>
      void callAndAddToTape(Type const* x, size_t m, Type* y, size_t n) {
        for (size_t i = 0; i < m; i += 1) {
          Base::addInput(x[i]);
        }

        for (size_t i = 0; i < n; i += 1) {
          Base::addOutput(y[i]);
        }

        callAndAddToTape<func>();
      }

    private:
      template<PrimalFunc func>
      static void enzymeDiff_b(Real const* x, Real* x_b, size_t m, Real const* y, Real const* y_b, size_t n,
                               ExternalFunctionUserData* d) {
        // clang-format off
        __enzyme_autodiff(
            (void*) func,
            enzyme_dup, x, x_b,
            enzyme_const, m,
            enzyme_dup, y, y_b,
            enzyme_const, n,
            enzyme_const, d);
        // clang-format on
      }

      template<PrimalFunc func>
      static void enzymeDiff_d(Real const* x, Real const* x_d, size_t m, Real* y, Real* y_d, size_t n,
                               ExternalFunctionUserData* d) {
        // clang-format off
        __enzyme_fwddiff(
            (void*) func,
            enzyme_dup, x, x_d,
            enzyme_const, m,
            enzyme_dup, y, y_d,
            enzyme_const, n,
            enzyme_const, d);
        // clang-format on
      }
  };
}
