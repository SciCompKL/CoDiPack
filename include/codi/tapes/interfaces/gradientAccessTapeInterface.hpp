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

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../misc/tapeParameters.hpp"

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

      /**
       * @brief Set the gradient.
       *
       * Automatic adjoints management involves bounds checking, resizing, and locking, see AdjointsManagement for
       * details.
       */
      void setGradient(Identifier const& identifier, Gradient const& gradient,
                       AdjointsManagement adjointsManagement = AdjointsManagement::Automatic);

      /**
       * @brief Set the gradient.
       *
       * Automatic adjoints management involves bounds checking and locking. If no adjoint variable with the given
       * identifier exists, a reference to adjoints[0] is returned. See AdjointsManagement for details.
       */
      Gradient const& getGradient(Identifier const& identifier,
                                  AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) const;

      /**
       * @brief Reference access to gradient.
       *
       * Automatic adjoints management involves bounds checking, resizing, and locking, see AdjointsManagement for
       * details.
       */
      Gradient& gradient(Identifier const& identifier,
                         AdjointsManagement adjointsManagement = AdjointsManagement::Automatic);

      /**
       * @brief Constant reference access to gradient.
       *
       * Automatic adjoints management involves bounds checking and locking. If no adjoint variable with the given
       * identifier exists, a reference to adjoints[0] is returned. See AdjointsManagement for details.
       */
      Gradient const& gradient(Identifier const& identifier,
                               AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) const;
  };
}
