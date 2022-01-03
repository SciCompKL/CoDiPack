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

#include "../../misc/macros.hpp"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Allow for a direct access to the gradient information computed by the tape.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * The gradient information is usually accessed via the helper functions of the ActiveType, for example
   * \code{.cpp}
   *   ActiveType<Tape> w = 1.0;
   *
   *   w.gradient() = 100.0;
   *
   *   std::cout << "Gradient of w: " << w.getGradient() << std::endl;
   * \endcode
   *
   * These helper function are shortcuts to the functions provided in this interface, but the functions here can also be
   * used to obtain the sensitivity information of a variable that is no longer present, for example
   * (documentation/examples/gradientAccessTapeInterface.cpp):
   * \snippet examples/gradientAccessTapeInterface.cpp Gradient Access
   *
   * Implementation hint: Size of the gradient vector should be checked before access.
   *
   * @tparam T_Gradient    The gradient type of a tape, usually chosen as ActiveType::Gradient.
   * @tparam T_Identifier  The adjoint/tangent identification of a tape, usually chosen as ActiveType::Identifier.
   */
  template<typename T_Gradient, typename T_Identifier>
  struct GradientAccessTapeInterface {
    public:

      using Gradient = CODI_DD(T_Gradient, double);   ///< See GradientAccessTapeInterface.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See GradientAccessTapeInterface.

      /*******************************************************************************/
      /// @name Interface definition

      void setGradient(Identifier const& identifier, Gradient const& gradient);  ///< Set the gradient.
      Gradient const& getGradient(Identifier const& identifier) const;           ///< Get the gradient.

      Gradient& gradient(Identifier const& identifier);              ///< Reference access to gradient.
      Gradient const& gradient(Identifier const& identifier) const;  ///< Constant reference access to gradient.
  };
}
