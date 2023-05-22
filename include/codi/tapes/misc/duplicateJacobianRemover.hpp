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

#include <array>

#include "../../config.h"
#include "../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Combines entries of Jacobians with the same identifier.
   *
   * This class is used in the storing process of the Jacobians for an expression. For each pushData, it checks if a
   * Jacobian with the same identifier has already been pushed. If so, then it combines these Jacobians.
   *
   * This behavior can be enabled with `-DCODI_RemoveDuplicateJacobianArguments=1`. See JacobianBaseTape::pushJacobians
   * for details.
   *
   * @tparam T_Real  The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_Identifier  The adjoint/tangent identifier type of a tape, usually chosen as ActiveType::Identifier.
   */
  template<typename T_Real, typename T_Identifier>
  struct DuplicateJacobianRemover {
    public:

      using Real = CODI_DD(T_Real, double);           ///< See DuplicateJacobianRemover.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See DuplicateJacobianRemover.
      using ArgumentSize = Config::ArgumentSize;      ///< Definition of ArgumentSize type.

    private:
      std::array<Identifier, Config::MaxArgumentSize> indices;
      std::array<Real, Config::MaxArgumentSize> jacobians;
      ArgumentSize size;

    public:

      /// Constructor
      DuplicateJacobianRemover() = default;

      /// For all added items, check if one matches the identifier. If yes combine, if no append.
      CODI_INLINE void pushData(Real const& jacobian, Identifier const& index) {
        bool found = false;
        ArgumentSize pos;
        for (pos = 0; pos < size; pos += 1) {
          if (indices[pos] == index) {
            found = true;
            break;
          }
        }

        if (!found) {
          size += 1;
          indices[pos] = index;
          jacobians[pos] = jacobian;
        } else {
          jacobians[pos] += jacobian;
        }
      }

      /// Add the data to the provided vector. Resets the internal data for a new statement push.
      /// @tparam Vec  DataInterface with Chunk2<Real, Identifier> as data.
      template<typename Vec>
      CODI_INLINE void storeData(Vec& vec) {
        for (ArgumentSize pos = 0; pos < size; pos += 1) {
          vec.pushData(jacobians[pos], indices[pos]);
        }

        // Reset the data for the next statement.
        size = 0;
      }
  };
}
