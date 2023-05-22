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
#include <cstddef>

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../../tools/data/direction.hpp"
#include "../../traits/realTraits.hpp"
#include "vectorAccessInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Implementation of VectorAccessInterface for adjoint vectors.
   *
   * The adjoint vector is used as is, it is assumed to have the correct size. No bounds checking is performed.
   *
   * @tparam T_Real        The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_Identifier  The adjoint/tangent identification of a tape, usually chosen as ActiveType::Identifier.
   * @tparam T_Gradient    The gradient type of a tape, usually chosen as ActiveType::Gradient.
   */
  template<typename T_Real, typename T_Identifier, typename T_Gradient>
  struct AdjointVectorAccess : public VectorAccessInterface<T_Real, T_Identifier> {
      using Real = CODI_DD(T_Real, double);           ///< See AdjointVectorAccess.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See AdjointVectorAccess.
      using Gradient = CODI_DD(T_Gradient, double);   ///< See AdjointVectorAccess.

    protected:

      Gradient* adjointVector;  ///< Pointer to the gradient vector.

    private:

      Gradient lhs;  ///< Temporary storage for indirect adjoint or tangent updates.

      std::array<Real, GradientTraits::dim<Gradient>()> buffer;  ///< Temporary storage for getAdjointVec.

    public:

      /// Constructor. See interface documentation for details about the adjoint vector.
      AdjointVectorAccess(Gradient* adjointVector) : adjointVector(adjointVector), lhs(), buffer() {}

      /*******************************************************************************/
      /// @name Misc

      /// \copydoc codi::VectorAccessInterface::getVectorSize
      size_t getVectorSize() const {
        return GradientTraits::dim<Gradient>();
      }

      /// \copydoc codi::VectorAccessInterface::isLhsZero
      bool isLhsZero() {
        return RealTraits::isTotalZero(lhs);
      }

      /// \copydoc codi::VectorAccessInterface::clone
      VectorAccessInterface<Real, Identifier>* clone() const {
        return new AdjointVectorAccess(this->adjointVector);
      }

      /*******************************************************************************/
      /// @name Indirect adjoint access

      /// \copydoc codi::VectorAccessInterface::setLhsAdjoint
      void setLhsAdjoint(Identifier const& index) {
        lhs = adjointVector[index];
        adjointVector[index] = Gradient();
      }

      /// \copydoc codi::VectorAccessInterface::updateAdjointWithLhs
      void updateAdjointWithLhs(Identifier const& index, Real const& jacobian) {
        adjointVector[index] += jacobian * lhs;
      }

      /*******************************************************************************/
      /// @name Indirect tangent access

      /// \copydoc codi::VectorAccessInterface::setLhsTangent
      void setLhsTangent(Identifier const& index) {
        adjointVector[index] = lhs;
        lhs = Gradient();
      }

      /// \copydoc codi::VectorAccessInterface::updateTangentWithLhs
      void updateTangentWithLhs(Identifier const& index, Real const& jacobian) {
        lhs += jacobian * adjointVector[index];
      }

      /*******************************************************************************/
      /// @name Direct adjoint access

      /// \copydoc codi::VectorAccessInterface::resetAdjoint
      void resetAdjoint(Identifier const& index, size_t dim) {
        GradientTraits::at(adjointVector[index], dim) = typename GradientTraits::Real<Gradient>();
      }

      /// \copydoc codi::VectorAccessInterface::resetAdjointVec
      void resetAdjointVec(Identifier const& index) {
        adjointVector[index] = Gradient();
      }

      /// \copydoc codi::VectorAccessInterface::getAdjoint
      Real getAdjoint(Identifier const& index, size_t dim) {
        CODI_UNUSED(dim);

        return (Real)GradientTraits::at(adjointVector[index], dim);
      }

      /// \copydoc codi::VectorAccessInterface::getAdjointVec
      void getAdjointVec(Identifier const& index, Real* const vec) {
        for (size_t i = 0; i < getVectorSize(); ++i) {
          vec[i] = (Real)GradientTraits::at(adjointVector[index], i);
        }
      }

      /// \copydoc codi::VectorAccessInterface::getAdjointVec
      Real const* getAdjointVec(Identifier const& index) {
        getAdjointVec(index, buffer.data());
        return buffer.data();
      }

      /// \copydoc codi::VectorAccessInterface::updateAdjoint
      void updateAdjoint(Identifier const& index, size_t dim, Real const& adjoint) {
        GradientTraits::at(adjointVector[index], dim) += adjoint;
      }

      /// \copydoc codi::VectorAccessInterface::updateAdjointVec
      void updateAdjointVec(Identifier const& index, Real const* const vec) {
        for (size_t i = 0; i < getVectorSize(); ++i) {
          GradientTraits::at(adjointVector[index], i) += vec[i];
        }
      }

      /*******************************************************************************/
      /// @name Primal access

      /// \copydoc codi::VectorAccessInterface::setPrimal <br><br>
      /// Implementation: Not implemented, empty function.
      void setPrimal(Identifier const& index, Real const& primal) {
        CODI_UNUSED(index, primal);
      }

      /// \copydoc codi::VectorAccessInterface::getPrimal <br><br>
      /// Implementation: Not implemented, returns zero.
      Real getPrimal(Identifier const& index) {
        CODI_UNUSED(index);

        return Real();
      }

      /// \copydoc codi::VectorAccessInterface::setPrimal <br><br>
      /// Implementation: Always returns false.
      bool hasPrimals() {
        return false;
      }
  };
}
