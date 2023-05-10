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

#include <algorithm>
#include <functional>
#include <type_traits>

#include "../../expressions/activeType.hpp"
#include "../../misc/exceptions.hpp"
#include "../../misc/macros.hpp"
#include "../../misc/memberStore.hpp"
#include "statementEvaluatorInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Only stores the function handle for the reverse evaluation.
   *
   * Uses the StatementEvaluatorTapeInterface.
   *
   * @tparam T_Real  The computation type of a tape, usually chosen as ActiveType::Real.
   */
  template<typename T_Real>
  struct ReverseStatementEvaluator : public StatementEvaluatorInterface<T_Real> {
      using Real = CODI_DD(T_Real, double);  ///< See ReverseStatementEvaluator.

    public:

      /*******************************************************************************/
      /// @name StatementEvaluatorInterface implementation
      /// @{

      using Handle = void*;  ///< Function pointer to the reverse evaluation.

      /// Throws CODI_EXCEPTION on call.
      template<typename Tape, typename... Args>
      static Real callForward(Handle const& h, Args&&... args) {
        CODI_UNUSED(h, args...);

        CODI_EXCEPTION("ReverseStatementEvaluator does not support forward evaluation calls.");

        return Real();
      }

      /// Throws CODI_EXCEPTION on call.
      template<typename Tape, typename... Args>
      static Real callPrimal(Handle const& h, Args&&... args) {
        CODI_UNUSED(h, args...);

        CODI_EXCEPTION("ReverseStatementEvaluator does not support primal evaluation calls.");

        return Real();
      }

      /// \copydoc StatementEvaluatorInterface::callReverse
      template<typename Tape, typename... Args>
      static void callReverse(Handle const& h, Args&&... args) {
        HandleTyped<Tape> func = (HandleTyped<Tape>)h;

        func(std::forward<Args>(args)...);
      }

      /// \copydoc StatementEvaluatorInterface::createHandle
      template<typename Tape, typename Generator, typename Expr>
      static Handle createHandle() {
        return (Handle*)Generator::template statementEvaluateReverse<Expr>;
      }

      /// @}

    protected:

      /// Full reverse function type.
      template<typename Tape>
      using HandleTyped = decltype(&Tape::template statementEvaluateReverse<ActiveType<Tape>>);
  };
}
