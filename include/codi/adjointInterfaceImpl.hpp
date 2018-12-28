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

#include "macros.h"
#include "tools/direction.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief The implementation assumes that each element in the adjoint vector consists only of one entry.
   *
   * Nearly everything of the base interface is implemented only the methods setPrimal and getPrimal are left out.
   *
   * @tparam          Real  The primal value of the CoDiPack type.
   * @tparam  GradientData  The identifier the CoDiPack type.
   * @tparam GradientValue  The adjoint value for the current evaluation. This type needs to support additions and
   *                        multiplications.
   */
  template<typename Real, typename GradientData, typename GradientValue>
  struct AdjointInterfaceImplBase : public AdjointInterface<Real, GradientData> {
      GradientValue* adjointVector; /**< The vector for the adjoint data.*/

      GradientValue lhs; /**< The stored value for the inplace updates. */

      /**
       * @brief Create a new instance.
       *
       * The vector is used for all operations.
       *
       * @param[in] adjointVector  The adjoint vector on which all the operations are evaluated.
       */
      explicit AdjointInterfaceImplBase(GradientValue* adjointVector) :
        adjointVector(adjointVector),
        lhs() {}

      /**
       * @brief Get the vector size of an adjoint value.
       * @return The vector size of an adjoint value.
       */
      size_t getVectorSize() const {
        return 1;
      }

      /**
       * @brief Set the adjoint value at the position and dimension to zero.
       *
       * @param[in] index  The position for the adjoint.
       * @param[in]   dim  The dimension in the vector.
       */
      void resetAdjoint(const GradientData index, const size_t dim) {
        CODI_UNUSED(dim);

        adjointVector[arrayAccess(index)] = GradientValue();
      }

      /**
       * @brief Set the adjoint vector at the position to zero.
       * @param[in] index  The position for the adjoint.
       */
      void resetAdjointVec(const GradientData index) {
        adjointVector[arrayAccess(index)] = GradientValue();
      }

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
      Real getAdjoint(const GradientData index, const size_t dim) {
        CODI_UNUSED(dim);

        return (Real) adjointVector[arrayAccess(index)];
      }

      /**
       * @brief Get the adjoint vector at the specified position.
       *
       * If the adjoint vector is a vector of vectors the result is:
       *
       * for(size_t i = 0; i < dim; ++i) {
       *   vec[dim] = adjoint[index][dim];
       * }
       *
       * where dim is the result from getVectorSize()
       *
       * @param[in] index  The position for the adjoint
       * @param[out]  vec  The vector for the storage of the data.
       */
      void getAdjointVec(const GradientData index, Real* vec) {
        *vec = (Real)adjointVector[arrayAccess(index)];
      }

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
      virtual void updateAdjoint(const GradientData index, const size_t dim, const Real adjoint) {
        CODI_UNUSED(dim);

        adjointVector[arrayAccess(index)] += adjoint;
      }

      /**
       * @brief Update the adjoint vector at the specified position.
       *
       * If the adjoint vector is a vector of vectors the update is:
       *
       * for(size_t i = 0; i < dim; ++i) {
       *   adjoint[index][dim] += vec[dim];
       * }
       *
       * where dim is the result from getVectorSize()
       *
       * @param[in] index  The position for the adjoint
       * @param[in]   vec  The update for the adjoint value.
       */
      virtual void updateAdjointVec(const GradientData index, const Real* vec) {
        adjointVector[arrayAccess(index)] += *vec;
      }

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
      void setLhsAdjoint(const GradientData index) {
        lhs = adjointVector[arrayAccess(index)];
      }

      /**
       * @brief Updates the target adjoint with the prior specified lhs multiplied with the jacobi.
       *
       * See also the documentation of setLhsAdjoint
       *
       * @param[in]  index  The index of the adjoint value that receives the update.
       * @param[in] jacobi  The jacobi value that is multiplied with the lhs adjoint.
       */
      void updateJacobiAdjoint(const GradientData index, Real jacobi) {
        adjointVector[arrayAccess(index)] += jacobi * lhs;
      }

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
      void setLhsTangent(const GradientData index) {
        adjointVector[arrayAccess(index)] = lhs;
        lhs = GradientValue();
      }

      /**
       * @brief Updates the lhs tangent with the jacobian multiplied with target tangent value.
       *
       * See also the documentation of setLhsTangent
       *
       * @param[in]  index  The index of the tangent value that is used for the update.
       * @param[in] jacobi  The jacobi value that is multiplied with the tangnet value defined by index.
       */
      void updateJacobiTangent(const GradientData index, Real jacobi) {
        lhs +=  jacobi * adjointVector[arrayAccess(index)];
      }
  };

  /**
   * @brief Specialization for the the codi::Direction structure.
   *
   * @tparam         Real  The primal value of the CoDiPack type.
   * @tparam GradientData  The identifier the CoDiPack type.
   * @tparam      RealDir  The type for the entries of the vectors. This type needs to support addition and multiplication
   *                       operations.
   * @tparam       vecDim  The dimension of the vector
   */
  template<typename Real, typename GradientData, typename RealDir, size_t vecDim>
  struct AdjointInterfaceImplBase <Real, GradientData, Direction<RealDir, vecDim> > : public AdjointInterface<Real, GradientData> {
      Direction<RealDir, vecDim>* adjointVector; /**< The vector for the adjoint data.*/

      Direction<RealDir, vecDim> lhs; /**< The stored value for the inplace updates. */

      /**
       * @brief Create a new instance.
       *
       * The vector is used for all operations.
       * @param[in] adjointVector  The adjoint vector on which all the operations are evaluated.
       */
      explicit AdjointInterfaceImplBase(Direction<RealDir, vecDim>* adjointVector) :
        adjointVector(adjointVector) {}

      /**
       * @brief Get the vector size of an adjoint value.
       * @return The vector size of an adjoint value.
       */
      size_t getVectorSize() const {
        return vecDim;
      }

      /**
       * @brief Set the adjoint value at the position and dimension to zero.
       *
       * @param[in] index  The position for the adjoint.
       * @param[in]   dim  The dimension in the vector.
       */
      void resetAdjoint(const GradientData index, const size_t dim) {
        adjointVector[arrayAccess(index)][dim] = RealDir();
      }

      /**
       * @brief Set the adjoint vector at the position to zero.
       * @param[in] index  The position for the adjoint.
       */
      void resetAdjointVec(const GradientData index) {
        adjointVector[arrayAccess(index)] = Direction<RealDir, vecDim>();
      }

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
      Real getAdjoint(const GradientData index, const size_t dim) {
        return (Real) adjointVector[arrayAccess(index)][dim];
      }

      /**
       * @brief Get the adjoint vector at the specified position.
       *
       * If the adjoint vector is a vector of vectors the result is:
       *
       * for(size_t i = 0; i < dim; ++i) {
       *   vec[dim] = adjoint[index][dim];
       * }
       *
       * where dim is the result from getVectorSize()
       *
       * @param[in] index  The position for the adjoint
       * @param[out]  vec  The vector for the storage of the data.
       */
      void getAdjointVec(const GradientData index, Real* vec) {
        for(size_t i = 0; i < vecDim; ++i) {
          vec[i] = (Real)adjointVector[arrayAccess(index)][i];
        }
      }

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
      virtual void updateAdjoint(const GradientData index, const size_t dim, const Real adjoint) {
        adjointVector[arrayAccess(index)][dim] += adjoint;
      }

      /**
       * @brief Update the adjoint vector at the specified position.
       *
       * If the adjoint vector is a vector of vectors the update is:
       *
       * for(size_t i = 0; i < dim; ++i) {
       *   adjoint[index][dim] += vec[dim];
       * }
       *
       * where dim is the result from getVectorSize()
       *
       * @param[in] index  The position for the adjoint
       * @param[in]   vec  The update for the adjoint value.
       */
      virtual void updateAdjointVec(const GradientData index, const Real* vec) {
        for(size_t i = 0; i < vecDim; ++i) {
          adjointVector[arrayAccess(index)][i] += vec[i];
        }
      }

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
      void setLhsAdjoint(const GradientData index) {
        lhs = adjointVector[arrayAccess(index)];
      }

      /**
       * @brief Updates the target adjoint with the prior specified lhs multiplied with the jacobi.
       *
       * See also the documentation of setLhsAdjoint
       *
       * @param[in]  index  The index of the adjoint value that receives the update.
       * @param[in] jacobi  The jacobi value that is multiplied with the lhs adjoint.
       */
      void updateJacobiAdjoint(const GradientData index, Real jacobi) {
        adjointVector[arrayAccess(index)] += jacobi * lhs;
      }

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
      void setLhsTangent(const GradientData index) {
        adjointVector[arrayAccess(index)] = lhs;
        lhs = Direction<RealDir, vecDim>();
      }

      /**
       * @brief Updates the lhs tangent with the jacobian multiplied with target tangent value.
       *
       * See also the documentation of setLhsTangent
       *
       * @param[in]  index  The index of the tangent value that is used for the update.
       * @param[in] jacobi  The jacobi value that is multiplied with the tangnet value defined by index.
       */
      void updateJacobiTangent(const GradientData index, Real jacobi) {
        lhs +=  jacobi * adjointVector[arrayAccess(index)];
      }
  };

  /**
   * @brief The implementation for tapes that do not require a primal value reset.
   *
   * @tparam          Real  The primal value of the CoDiPack type.
   * @tparam  GradientData  The identifier the CoDiPack type.
   * @tparam GradientValue  The adjoint value for the current evaluation. This type needs to support additions and
   *                        multiplications.
   */
  template<typename Real, typename GradientData, typename GradientValue>
  struct AdjointInterfaceImpl final : public AdjointInterfaceImplBase<Real, GradientData, GradientValue> {

      /**
       * @brief Create a new instance.
       *
       * The vector is used for all operations.
       *
       * @param[in] adjointVector  The adjoint vector on which all the operations are evaluated.
       */
      explicit AdjointInterfaceImpl(GradientValue* adjointVector) :
        AdjointInterfaceImplBase<Real, GradientData, GradientValue>(adjointVector) {}

      /**
       * @brief Some tapes need to revert the primal values in the primal value vector to the old value
       * for output variables.
       *
       * If the tape needs this behaviour can be checked with Tape::RequiresPrimalReset. The value required
       * here is returned on a registerExtFunctionOutput call.
       *
       * @param[in]  index  The index of the primal value that needs to be reverted.
       * @param[in] primal  The primal value that is set.
       */
      virtual void setPrimal(const GradientData index, Real primal) {
        CODI_UNUSED(index);
        CODI_UNUSED(primal);

        // no primal handling required for the tape
      }

      /**
       * @brief Get primal value support for the primal evaluation of primal value tapes.
       *
       * @param[in]  index  The index of the primal value
       */
      Real getPrimal(const GradientData index) {
        CODI_UNUSED(index);

        // no primal handling required for the tape
        return Real();
      }
  };

  /**
   * @brief The implementation for tapes that require a primal value reset.
   *
   * @tparam          Real  The primal value of the CoDiPack type.
   * @tparam  GradientData  The identifier the CoDiPack type.
   * @tparam GradientValue  The adjoint value for the current evaluation. This type needs to support additions and
   *                        multiplications.
   */
  template<typename Real, typename GradientData, typename GradientValue>
  struct AdjointInterfacePrimalImpl final : public AdjointInterfaceImplBase<Real, GradientData, GradientValue> {

      Real* primalVector; /**< The vector for the primal values */

      /**
       * @brief Create a new instance.
       *
       * The adjoint vector is used for all adjoint operations.
       * The primal vector is used for all primal operations.
       *
       * @param[in] adjointVector  The adjoint vector on which all the operations are evaluated.
       * @param[in]  primalVector  The primal vector on which the primal reset is done.
       */
      explicit AdjointInterfacePrimalImpl(GradientValue* adjointVector, Real* primalVector) :
        AdjointInterfaceImplBase<Real, GradientData, GradientValue>(adjointVector),
        primalVector(primalVector) {}

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
      void setPrimal(const GradientData index, Real primal) {
        primalVector[index] = primal;
      }

      /**
       * @brief Get primal value support for the primal evaluation of primal value tapes.
       *
       * @param[in]  index  The index of the primal value
       */
      Real getPrimal(const GradientData index) {
        return primalVector[index];
      }
  };
}
