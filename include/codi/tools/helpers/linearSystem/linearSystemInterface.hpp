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
#include "../../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * Types for the LinearSystemInterface. These are passed as template parameters to other structures and therefore they
   * are combined in an extra definition.
   *
   * Implementations of LinearSystemInterface must also define their own set of types with this structure.
   */
  struct LinearSystemInterfaceTypes {
    public:

      using Type = CODI_ANY;  ///< A CoDiPack type or a floating point type.

      using Matrix = CODI_ANY;  ///< The Matrix with Type as the computation type (e.g. M<Type>).
      using MatrixReal =
          CODI_ANY;  ///< The Matrix with the Real type of the CoDiPack type Type (e.g. M<typename Type::Real>).
      using MatrixIdentifier = CODI_ANY;  ///< The Matrix with the Identifier type of the CoDiPack type Type (e.g.
                                          ///< M<typename Type::Identifier>).
      using Vector = CODI_ANY;            ///< The Vector with Type as the computation type (e.g. V<Type>).
      using VectorReal =
          CODI_ANY;  ///< The Vector with the Real type of the CoDiPack type Type (e.g. V<typename Type::Real>).
      using VectorIdentifier = CODI_ANY;  ///< The Vector with the Identifier type of the CoDiPack type Type (e.g.
                                          ///< V<typename Type::Identifier>).
  };

  /**
   * The interface defines all mandatory and optional functions that are required by the LinearSystemSolverHandler.
   * Implementations need to define all mandatory functions, the optional ones either depend on the use case and on
   * specifics of the algorithm. For more details see LinearSystemSolverHandler.
   *
   * If the interface has been specialized, it can be called with:
   * \code{.cpp}
   * codi::solveLinearSystem({Specialization}, A, rhs, sol, hints);
   * \endcode
   *  The hints parameter is optional, but can be used to improve the runtime and memory.
   *  The set of flags that can be passed as hints is defined in #LinearSystemSolverFlags.
   *  - ReverseEvaluation: Prepare for a reverse mode evaluation. Stores A_v_trans.
   *  - ForwardEvaluation: Prepare for a forward mode evaluation. Stores A_v.
   *  - PrimalEvaluation:  Prepare for a primal reevaluation. Stores A_v.
   *  - ProvidePrimalSolution: Read x_v before the system is solved and provide it to the solveSystem or
   *                           solveSystemPrimal methods (only during the primal computation, not in the external
   *                           function implementations).
   *  - RecomputePrimalInForwardEvaluation: In the AD forward mode also solve the primal linear system again.
   *
   *  See \ref Example_21_Special_handling_of_linear_system_solvers for an example with the Eigen implementation.
   *
   *  \subsection sec_mandatory Mandatory methods
   *   - #createMatrixReal, #createMatrixIdentifier, #deleteMatrixReal, #deleteMatrixIdentifier
   *   - #createVectorReal, #createVectorIdentifier, #deleteVectorReal, #deleteVectorIdentifier
   *   - #iterateMatrix with two and three arguments
   *   - #iterateVector with two to four arguments
   *   - #solveSystem
   *
   *   \subsection sec_optional  Optional methods
   *    - Reverse mode AD support:
   *      - #iterateDyadic
   *      - #transposeMatrix
   *    - Foward mode AD support:
   *      - #subtractMultiply
   *    - Other:
   *      - #solveSystemPrimal
   *
   * @tparam T_InterfaceTypes  The definition of LinearSystemInterfaceTypes for the implementation.
   */
  template<typename T_InterfaceTypes>
  struct LinearSystemInterface {
    public:

      using InterfaceTypes = CODI_DD(T_InterfaceTypes, LinearSystemInterfaceTypes);  ///< See LinearSystemInterface.

      using Type = typename InterfaceTypes::Type;  ///< See LinearSystemInterfaceTypes.

      using Matrix = typename InterfaceTypes::Matrix;                      ///< See LinearSystemInterfaceTypes.
      using MatrixReal = typename InterfaceTypes::MatrixReal;              ///< See LinearSystemInterfaceTypes.
      using MatrixIdentifier = typename InterfaceTypes::MatrixIdentifier;  ///< See LinearSystemInterfaceTypes.
      using Vector = typename InterfaceTypes::Vector;                      ///< See LinearSystemInterfaceTypes.
      using VectorReal = typename InterfaceTypes::VectorReal;              ///< See LinearSystemInterfaceTypes.
      using VectorIdentifier = typename InterfaceTypes::VectorIdentifier;  ///< See LinearSystemInterfaceTypes.

      /*******************************************************************************/
      /// @name Mandatory: Implementations for matrix and vector creation and deletion.
      /// @{

      /// Create a real matrix from an existing one. Values do not need to be copied.
      /// @tparam M  M is either Matrix or MatrixIdentifier.
      template<typename M>
      MatrixReal* createMatrixReal(M* mat);

      /// Create an identifier matrix from an existing one. Values do not need to be copied.
      /// @tparam M  M is either Matrix or MatrixIdentifier.
      template<typename M>
      MatrixIdentifier* createMatrixIdentifier(M* mat);

      /// Create a real vector from an existing one. Values do not need to be copied.
      /// @tparam V  V is either Vector or VectorIdentifier.
      template<typename V>
      VectorReal* createVectorReal(V* vec);

      /// Create an identifier vector from an existing one. Values do not need to be copied.
      /// @tparam V  V is either Vector or VectorIdentifier.
      template<typename V>
      VectorIdentifier* createVectorIdentifier(V* vec);

      /// Delete a real matrix.
      void deleteMatrixReal(MatrixReal* A_v);

      /// Delete an identifier matrix.
      void deleteMatrixIdentifier(MatrixIdentifier* A_id);

      /// Delete a real Vector.
      void deleteVectorReal(VectorReal* vec_v);

      /// Delete an identifier vector.
      void deleteVectorIdentifier(VectorIdentifier* vec_id);

      /// @}
      /*******************************************************************************/
      /// @name Mandatory: Implementations for matrix and vector iterations.
      /// @{

      /// Iterate over all elements in the matrices at the same time.
      template<typename Func, typename MatrixA, typename MatrixB>
      void iterateMatrix(Func func, MatrixA* matA, MatrixB* matB);

      /// Iterate over all elements in the matrices at the same time.
      template<typename Func, typename MatrixA, typename MatrixB, typename MatrixC>
      void iterateMatrix(Func func, MatrixA* matA, MatrixB* matB, MatrixC* matC);

      /// Iterate over all elements in the vectors at the same time.
      template<typename Func, typename VectorA, typename VectorB>
      void iterateVector(Func func, VectorA* vecA, VectorB* vecB);

      /// Iterate over all elements in the vectors at the same time.
      template<typename Func, typename VectorA, typename VectorB, typename VectorC>
      void iterateVector(Func func, VectorA* vecA, VectorB* vecB, VectorC* vecC);

      /// Iterate over all elements in the vectors at the same time.
      template<typename Func, typename VectorA, typename VectorB, typename VectorC, typename VectorD>
      void iterateVector(Func func, VectorA* vecA, VectorB* vecB, VectorC* vecC, VectorD* vecD);

      /// @}
      /*******************************************************************************/
      /// @name Mandatory: Implementations for the linear system solve.
      /// @{

      /// Solve the linear system with the real valued matrices and vectors.
      /// Solves Ax = b for x.
      void solveSystem(MatrixReal const* A, VectorReal const* b, VectorReal* x);

      /// @}
      /*******************************************************************************/
      /// @name Optional: Implementation for reverse mode AD.
      /// @{

      /// Iterate over all elements in \c mat_id \c and provide the elements in \c b_b \c and \c x_v \c .
      /// For element (i,j) func needs to be called with \c func(mat_id(i,j), b_b(i), x_v(j)) \c .
      /// For sparse matrices, only the elements of the sparsity pattern need to be considered.
      /// Used for e.g. the computation of the dyadic product \f$A=x_v*b_b^T\f$.
      template<typename Func>
      void iterateDyadic(Func func, MatrixIdentifier* mat_id, VectorReal* x_v, VectorReal* b_b) {
        CODI_UNUSED(func, mat_id, x_v, b_b);
      }

      /// Create a transposed matrix.
      MatrixReal* transposeMatrix(MatrixReal* A_v) {
        CODI_UNUSED(A_v);
        return NULL;
      }

      /// @}
      /*******************************************************************************/
      /// @name Optional: Implementation for forward mode AD.
      /// @{

      /// Computes t = b_d - A_d * x.
      void subtractMultiply(VectorReal* t, VectorReal const* b_d, MatrixReal const* A_d, VectorReal const* x) {
        CODI_UNUSED(t, b_d, A_d, x);
      }

      /// @}
      /*******************************************************************************/
      /// @name Optional: Implementations for specializations of the algorithm.
      /// @{

      /// Solve the linear system with the real valued matrices and vectors.
      /// Implementation that is called in the primal routine. If not specialized solveSystem is called.
      /// Solves Ax = b for x.
      void solveSystemPrimal(MatrixReal const* A, VectorReal const* b, VectorReal* x) {
        CODI_UNUSED(A, b, x);
      }

      /// @}
  };
}
