/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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
#pragma once

#include <cstddef>

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "adjointVectorAccess.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Implementation of VectorAccessInterface for adjoint and primal vectors.
   *
   * Both vectors are used as is, they are assumed to have correct sizes. No bounds checking is performed.
   *
   * Inherits from AdjointVectorAccess and overwrites all methods specific to the primals.
   *
   * @tparam T_Real           The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_Identifier     The adjoint/tangent identification of a tape, usually chosen as ActiveType::Identifier.
   * @tparam T_AdjointVector  Either a pointer type, for example Gradient*, or a reference to something that can be
   *                          accessed like a gradient array, for example std::vector<Gradient>&.
   */
  template<typename T_Real, typename T_Identifier, typename T_AdjointVector>
  struct PrimalAdjointVectorAccess : public AdjointVectorAccess<T_Real, T_Identifier, T_AdjointVector> {
      using Real = CODI_DD(T_Real, double);           ///< See PrimalAdjointVectorAccess.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See PrimalAdjointVectorAccess.
      using Gradient = CODI_DD(AdjointVectorTraits::Gradient<T_AdjointVector>, double);  ///< Adjoint vector entry type.
      using AdjointVector = CODI_DD(T_AdjointVector, double*);  ///< See PrimalAdjointVectorAccess.

      using Base = AdjointVectorAccess<Real, Identifier, AdjointVector>;  ///< Base class abbreviation.

    private:

      Real* primalVector;  ///< Pointer to the primal vector.

    public:

      /// Constructor. See interface documentation for details about the vectors.
      PrimalAdjointVectorAccess(AdjointVector adjointVector, Real* primalVector)
          : Base(adjointVector), primalVector(primalVector) {}

      /*******************************************************************************/
      /// @name Misc

      /// \copydoc codi::VectorAccessInterface::clone
      VectorAccessInterface<Real, Identifier>* clone() const {
        return new PrimalAdjointVectorAccess(this->adjointVector, this->primalVector);
      }

      /*******************************************************************************/
      /// @name Primal access

      /// \copydoc VectorAccessInterface::setPrimal
      void setPrimal(Identifier const& index, Real const& primal) {
        primalVector[index] = primal;
      }

      /// \copydoc VectorAccessInterface::getPrimal
      Real getPrimal(Identifier const& index) {
        return primalVector[index];
      }

      /// \copydoc VectorAccessInterface::hasPrimals <br>
      /// Implementation: Always returns true.
      bool hasPrimals() {
        return true;
      }
  };
}
