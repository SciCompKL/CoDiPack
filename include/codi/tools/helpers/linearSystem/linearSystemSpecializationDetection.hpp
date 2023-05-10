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

#include <vector>

#include "../../../config.h"
#include "../../../expressions/lhsExpressionInterface.hpp"
#include "../../../misc/exceptions.hpp"
#include "../../../traits/tapeTraits.hpp"
#include "../../data/direction.hpp"
#include "linearSystemFlags.hpp"
#include "linearSystemInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * Detects if certain methods in a LinearSystemInterface have been specialized.
   *
   * @tparam T_LinearSystem  Implementation of LinearSystemInterface.
   */
  template<typename T_LinearSystem>
  struct LinearSystemSpecializationDetection {
      using LinearSystem =
          CODI_DD(T_LinearSystem,
                  CODI_T(LinearSystemInterface<LinearSystemInterfaceTypes>));  ///< See LinearSystemOverloadDetection.

      using Type = CODI_DD(typename LinearSystem::Type,
                           CODI_DEFAULT_LHS_EXPRESSION);  ///< See LinearSystemInterfaceTypes.

      using Matrix = typename LinearSystem::Matrix;                      ///< See LinearSystemInterfaceTypes.
      using MatrixReal = typename LinearSystem::MatrixReal;              ///< See LinearSystemInterfaceTypes.
      using MatrixIdentifier = typename LinearSystem::MatrixIdentifier;  ///< See LinearSystemInterfaceTypes.
      using Vector = typename LinearSystem::Vector;                      ///< See LinearSystemInterfaceTypes.
      using VectorReal = typename LinearSystem::VectorReal;              ///< See LinearSystemInterfaceTypes.
      using VectorIdentifier = typename LinearSystem::VectorIdentifier;  ///< See LinearSystemInterfaceTypes.

    private:

      using Real = typename Type::Real;
      using Identifier = typename Type::Identifier;

      using Interface = LinearSystemInterface<typename LinearSystem::InterfaceTypes>;

      static void dyadicProxy(Identifier&, Real&, Real&) {}

    public:

      // gcc issues false warning for comparison with NULL.
#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Waddress"
#endif

      /// Checks if iterateDyadic is specialized in LinearSystem.
      CODI_INLINE static bool IsDyadicImplemented() {
        return static_cast<void (Interface::*)(decltype(dyadicProxy), MatrixIdentifier*, VectorReal*, VectorReal*)>(
                   &LinearSystem::iterateDyadic) !=
               static_cast<void (Interface::*)(decltype(dyadicProxy), MatrixIdentifier*, VectorReal*, VectorReal*)>(
                   &Interface::iterateDyadic);
      }

      /// Checks if transposeMatrix is specialized in LinearSystem.
      CODI_INLINE static bool IsTransposeImplemented() {
        return &LinearSystem::transposeMatrix != &Interface::transposeMatrix;
      }

      /// Checks if subtractMultiply is specialized in LinearSystem.
      CODI_INLINE static bool IsSubtractMultiplyImplemented() {
        return &LinearSystem::subtractMultiply != &Interface::subtractMultiply;
      }

      /// Checks if solveSystemPrimal is specialized in LinearSystem.
      CODI_INLINE static bool IsSolvePrimalImplemented() {
        return &LinearSystem::solveSystemPrimal != &Interface::solveSystemPrimal;
      }

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif

      /// True if all functions for the reverse mode support are specialized.
      CODI_INLINE static bool SupportsReverseMode() {
        return IsDyadicImplemented() && IsTransposeImplemented();
      }

      /// True if all functions for the forward mode support are specialized.
      CODI_INLINE static bool SupportsForwardMode() {
        return IsSubtractMultiplyImplemented();
      }
  };
}
