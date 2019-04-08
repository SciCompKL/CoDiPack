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
 */

#pragma once

#include <cstddef>

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   *
   * @tparam Real  The floating point type that is used in the tape for the computation.
   * @tparam Index The identifier that is used in the tape for the identification of variables.
   */
  template<typename Real, typename Index>
  struct AdjointInterface {

      /**
       * Destructor for the interface.
       */
      virtual ~AdjointInterface() {}

      /**
       * @brief Get the vector size of an adjoint value.
       * @return The vector size of an adjoint value.
       */
      virtual size_t getVectorSize() const = 0;

      /**
       * @brief Set the adjoint value at the position and dimension to zero.
       *
       * @param[in] index  The position for the adjoint.
       * @param[in]   dim  The dimension in the vector.
       */
      virtual void resetAdjoint(const Index index, const size_t dim) = 0;

      /**
       * @brief Set the adjoint vector at the position to zero.
       * @param[in] index  The position for the adjoint.
       */
      virtual void resetAdjointVec(const Index index) = 0;

      /**
       * @brief Get the adjoint value at the specified position and dimension.
       *
       * If the adjoint vector is a vector of vectors the result is:
       *
       *  adjoint[index][dim]
       *
       * @param[in] index  The position for the adjoint
       * @param[in]   dim  The dimension in the vector.
       * @return The adjoint value at the position with the dimension.
       */
      virtual Real getAdjoint(const Index index, const size_t dim) = 0;

      /**
       * @brief Get the adjoint vector at the specified position.
       *
       * If the adjoint vector is a vector of vectors the result is:
       *
       * for(size_t i = 0; i < dim; ++i) {
       *   vec[i] = adjoint[index][i];
       * }
       *
       * where dim is the result from getVectorSize()
       *
       * @param[in] index  The position for the adjoint
       * @param[out]  vec  The vector for the storage of the data.
       */
      virtual void getAdjointVec(const Index index, Real* vec) = 0;

      /**
       * @brief Update the adjoint value at the specified position and dimension.
       *
       * If the adjoint vector is a vector of vectors the update is:
       *
       *  adjoint[index][dim] += adjoint;
       *
       * @param[in]   index  The position for the adjoint
       * @param[in]     dim  The dimension in the vector.
       * @param[in] adjoint  The update for the adjoint value.
       */
      virtual void updateAdjoint(const Index index, const size_t dim, const Real adjoint) = 0;

      /**
       * @brief Update the adjoint vector at the specified position.
       *
       * If the adjoint vector is a vector of vectors the update is:
       *
       * for(size_t i = 0; i < dim; ++i) {
       *   adjoint[index][i] += vec[i];
       * }
       *
       * where dim is the result from getVectorSize()
       *
       * @param[in]   index  The position for the adjoint
       * @param[in]     vec  The update for the adjoint vector.
       */
      virtual void updateAdjointVec(const Index index, const Real* vec) = 0;

      /**
       * @brief The adjoint target for the adjoint of the left hand side of an equation.
       *
       * The function needs to be used together with updateJacobiAdjoint. For the statement
       *
       *  w = h(x)
       *
       * the adjoint update
       *
       * x_b += jac * w_b
       * w_b = 0.0
       *
       * needs to be performed where jac = dh/dx. This function call identifies w_b by the index and stores it internally.
       * A call to resetAdjointVec with the same index can then be used to reset w_b to zero. With the function call to
       * updateJacobiAdjoint the multiplication jac * w_b is performed and the adjoint identified with the index is
       * updated.
       *
       * @param[in] index  The index of the adjoint value that is stored and reset to zero.
       */
      virtual void setLhsAdjoint(const Index index) = 0;

      /**
       * @brief Updates the target adjoint with the prior specified lhs multiplied with the jacobi.
       *
       * See also the documentation of setLhsAdjoint
       *
       * @param[in]  index  The index of the adjoint value that receives the update.
       * @param[in] jacobi  The jacobi value that is multiplied with the lhs adjoint.
       */
      virtual void updateJacobiAdjoint(const Index index, Real jacobi) = 0;

      /**
       * @brief The tangent target for the tangent of the left hand side of an equation.
       *
       * The function needs to be used together with updateJacobiTangent. For the statement
       *
       *  w = h(x)
       *
       * the tangent update
       *
       * w_d += jac * x_d
       *
       * needs to be performed where jac = dh/dx. This function call identifies w_d by the index and sets it to the currently
       * accumulated value. With the function call to updateJacobiAdjoint the multiplication jac * x_d is performed and needs
       * to be called before the call to this function.
       *
       * @param[in] index  The index of the tangent value that is set to the current accumulated value.
       */
      virtual void setLhsTangent(const Index index) = 0;

      /**
       * @brief Updates the lhs tangent with the jacobian multiplied with target tangent value.
       *
       * See also the documentation of setLhsTangent
       *
       * @param[in]  index  The index of the tangent value that is used for the update.
       * @param[in] jacobi  The jacobi value that is multiplied with the tangnet value defined by index.
       */
      virtual void updateJacobiTangent(const Index index, Real jacobi) = 0;

      /**
       * @brief Some tapes need to revert the primal values in the primal value vector to the old value
       * for output variables.
       *
       * If the tape needs this behaviour can be checked with Tape::RequiresPrimalReset. The value required
       * here is returned on a registerExtFunctionOutput call.
       *
       * The function can also be used to enable a primal evaluation for primal value tapes.
       *
       * @param[in]  index  The index of the primal value that needs to be reverted.
       * @param[in] primal  The primal value that is set.
       */
      virtual void setPrimal(const Index index, Real primal) = 0;

      /**
       * @brief Get primal value support for the primal evaluation of primal value tapes.
       *
       * @param[in]  index  The index of the primal value
       */
      virtual Real getPrimal(const Index index) = 0;
  };
}
