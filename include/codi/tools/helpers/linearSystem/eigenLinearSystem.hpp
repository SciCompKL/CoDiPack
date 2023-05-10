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

#if CODI_EnableEigen

  #include <Eigen/Eigen>
  #include <vector>

  #include "../../../config.h"
  #include "../../../misc/macros.hpp"
  #include "linearSystemHandler.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Eigen definition for the LinearSystemInterfaceTypes.
  template<typename T_Type, template<typename> class T_Matrix, template<typename> class T_Vector>
  struct EigenLinearSystemTypes {
    public:

      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See LinearSystemInterfaceTypes.

      using Matrix = CODI_DD(CODI_T(T_Matrix<Type>),
                             CODI_T(Eigen::Matrix<Type, 2, 2>));  ///< See LinearSystemInterfaceTypes.
      using Vector = CODI_DD(CODI_T(T_Vector<Type>),
                             CODI_T(Eigen::Matrix<Type, 2, 1>));  ///< See LinearSystemInterfaceTypes.

      using Real = typename Type::Real;              ///< See LhsExpressionInterface.
      using Identifier = typename Type::Identifier;  ///< See LhsExpressionInterface.

      using MatrixReal = CODI_DD(CODI_T(T_Matrix<Real>),
                                 CODI_T(Eigen::Matrix<Real, 2, 2>));  ///< See LinearSystemInterfaceTypes.
      using VectorReal = CODI_DD(CODI_T(T_Vector<Real>),
                                 CODI_T(Eigen::Matrix<Real, 2, 1>));  ///< See LinearSystemInterfaceTypes.

      using MatrixIdentifier = CODI_DD(CODI_T(T_Matrix<Identifier>),
                                       CODI_T(Eigen::Matrix<Identifier, 2, 2>));  ///< See LinearSystemInterfaceTypes.
      using VectorIdentifier = CODI_DD(CODI_T(T_Vector<Identifier>),
                                       CODI_T(Eigen::Matrix<Identifier, 2, 1>));  ///< See LinearSystemInterfaceTypes.
  };

  /// Eigen implementation of LinearSystemInterface. The only methods missing are
  /// solveSystem and solveSystemPrimal (optional). TODO: Link example
  template<typename T_Type, template<typename> class T_Matrix, template<typename> class T_Vector>
  struct EigenLinearSystem : public LinearSystemInterface<EigenLinearSystemTypes<T_Type, T_Matrix, T_Vector>> {
    public:
      using InterfaceTypes = EigenLinearSystemTypes<T_Type, T_Matrix, T_Vector>;  ///< See LinearSystemInterface.

      using Type = typename InterfaceTypes::Type;  ///< See LinearSystemInterfaceTypes.

      using Matrix = typename InterfaceTypes::Matrix;                      ///< See LinearSystemInterfaceTypes.
      using MatrixReal = typename InterfaceTypes::MatrixReal;              ///< See LinearSystemInterfaceTypes.
      using MatrixIdentifier = typename InterfaceTypes::MatrixIdentifier;  ///< See LinearSystemInterfaceTypes.
      using Vector = typename InterfaceTypes::Vector;                      ///< See LinearSystemInterfaceTypes.
      using VectorReal = typename InterfaceTypes::VectorReal;              ///< See LinearSystemInterfaceTypes.
      using VectorIdentifier = typename InterfaceTypes::VectorIdentifier;  ///< See LinearSystemInterfaceTypes.

      using Index = typename Matrix::Index;  ///< Index of an Eigen matrix.

      /*******************************************************************************/
      /// @name Mandatory: Implementations for matrix and vector creation and deletion.
      /// @{

      /// \copydoc codi::LinearSystemInterface::createMatrixReal
      template<typename M>
      MatrixReal* createMatrixReal(M* mat) {
        return new MatrixReal(mat->rows(), mat->cols());
      }

      /// \copydoc codi::LinearSystemInterface::createMatrixIdentifier
      template<typename M>
      MatrixIdentifier* createMatrixIdentifier(M* mat) {
        return new MatrixIdentifier(mat->rows(), mat->cols());
      }

      /// \copydoc codi::LinearSystemInterface::createVectorReal
      template<typename V>
      VectorReal* createVectorReal(V* vec) {
        return new VectorReal(vec->size());
      }

      /// \copydoc codi::LinearSystemInterface::createVectorIdentifier
      template<typename V>
      VectorIdentifier* createVectorIdentifier(V* vec) {
        return new VectorIdentifier(vec->size());
      }

      /// \copydoc codi::LinearSystemInterface::deleteMatrixReal
      void deleteMatrixReal(MatrixReal* A_v) {
        delete A_v;
      }

      /// \copydoc codi::LinearSystemInterface::deleteMatrixIdentifier
      void deleteMatrixIdentifier(MatrixIdentifier* A_id) {
        delete A_id;
      }

      /// \copydoc codi::LinearSystemInterface::deleteVectorReal
      void deleteVectorReal(VectorReal* vec_v) {
        delete vec_v;
      }

      /// \copydoc codi::LinearSystemInterface::deleteVectorIdentifier
      void deleteVectorIdentifier(VectorIdentifier* vec_id) {
        delete vec_id;
      }

      /// @}
      /*******************************************************************************/
      /// @name Mandatory: Implementations for matrix and vector iterations.
      /// @{

      /// \copydoc codi::LinearSystemInterface::iterateMatrix
      template<typename Func, typename MatrixA, typename MatrixB>
      void iterateMatrix(Func func, MatrixA* matA, MatrixB* matB) {
        Index rows = matA->rows();
        Index cols = matA->cols();

        for (int i = 0; i < rows; i += 1) {
          for (int j = 0; j < cols; j += 1) {
            func(matA->coeffRef(i, j), matB->coeffRef(i, j));
          }
        }
      }

      /// \copydoc codi::LinearSystemInterface::iterateMatrix
      template<typename Func, typename MatrixA, typename MatrixB, typename MatrixC>
      void iterateMatrix(Func func, MatrixA* matA, MatrixB* matB, MatrixC* matC) {
        Index rows = matA->rows();
        Index cols = matA->cols();

        for (int i = 0; i < rows; i += 1) {
          for (int j = 0; j < cols; j += 1) {
            func(matA->coeffRef(i, j), matB->coeffRef(i, j), matC->coeffRef(i, j));
          }
        }
      }

      /// \copydoc codi::LinearSystemInterface::iterateVector
      template<typename Func, typename VectorA, typename VectorB>
      void iterateVector(Func func, VectorA* vecA, VectorB* vecB) {
        Index size = vecA->size();

        for (int i = 0; i < size; i += 1) {
          func(vecA->coeffRef(i), vecB->coeffRef(i));
        }
      }

      /// \copydoc codi::LinearSystemInterface::iterateVector
      template<typename Func, typename VectorA, typename VectorB, typename VectorC>
      void iterateVector(Func func, VectorA* vecA, VectorB* vecB, VectorC* vecC) {
        Index size = vecA->size();

        for (int i = 0; i < size; i += 1) {
          func(vecA->coeffRef(i), vecB->coeffRef(i), vecC->coeffRef(i));
        }
      }

      /// \copydoc codi::LinearSystemInterface::iterateVector
      template<typename Func, typename VectorA, typename VectorB, typename VectorC, typename VectorD>
      void iterateVector(Func func, VectorA* vecA, VectorB* vecB, VectorC* vecC, VectorD* vecD) {
        Index size = vecA->size();

        for (int i = 0; i < size; i += 1) {
          func(vecA->coeffRef(i), vecB->coeffRef(i), vecC->coeffRef(i), vecD->coeffRef(i));
        }
      }

      /// @}
      /*******************************************************************************/
      /// @name Mandatory: Implementations for the linear system solve.
      /// @{

      /// \copydoc codi::LinearSystemInterface::solveSystem <br> Needs to be implemented by the user.
      void solveSystem(MatrixReal const* A, VectorReal const* b, VectorReal* x);

      /// @}
      /*******************************************************************************/
      /// @name Implementation for reverse mode AD.
      /// @{

      /// \copydoc codi::LinearSystemInterface::transposeMatrix
      MatrixReal* transposeMatrix(MatrixReal* A_v) {
        MatrixReal* A_v_trans = new MatrixReal(A_v->rows(), A_v->cols());
        *A_v_trans = A_v->transpose();

        return A_v_trans;
      }

      /// \copydoc codi::LinearSystemInterface::iterateDyadic
      template<typename Func>
      void iterateDyadic(Func func, MatrixIdentifier* mat_id, VectorReal* x_v, VectorReal* b_b) {
        Index rows = mat_id->rows();
        Index cols = mat_id->cols();

        for (int i = 0; i < rows; i += 1) {
          for (int j = 0; j < cols; j += 1) {
            func(mat_id->coeffRef(i, j), x_v->coeffRef(j), b_b->coeffRef(i));
          }
        }
      }

      /// @}
      /*******************************************************************************/
      /// @name Implementation for forward mode AD.
      /// @{

      /// \copydoc codi::LinearSystemInterface::subtractMultiply
      void subtractMultiply(VectorReal* t, VectorReal const* b_d, MatrixReal const* A_d, VectorReal const* x) {
        *t = *b_d - *A_d * *x;
      }

      /// @}
  };

  /// Eigen implementation of LinearSystemInterface for sparse matrices. The only methods missing are
  /// solveSystem and solveSystemPrimal (optional).  TODO: Link example
  template<typename T_Type, template<typename> class T_Matrix, template<typename> class T_Vector>
  struct SparseEigenLinearSystem : public EigenLinearSystem<T_Type, T_Matrix, T_Vector> {
    public:
      using InterfaceTypes = EigenLinearSystemTypes<T_Type, T_Matrix, T_Vector>;  ///< See LinearSystemInterface.

      using Type = typename InterfaceTypes::Type;  ///< See LinearSystemInterfaceTypes.

      using Matrix = typename InterfaceTypes::Matrix;                      ///< See LinearSystemInterfaceTypes.
      using MatrixReal = typename InterfaceTypes::MatrixReal;              ///< See LinearSystemInterfaceTypes.
      using MatrixIdentifier = typename InterfaceTypes::MatrixIdentifier;  ///< See LinearSystemInterfaceTypes.
      using Vector = typename InterfaceTypes::Vector;                      ///< See LinearSystemInterfaceTypes.
      using VectorReal = typename InterfaceTypes::VectorReal;              ///< See LinearSystemInterfaceTypes.
      using VectorIdentifier = typename InterfaceTypes::VectorIdentifier;  ///< See LinearSystemInterfaceTypes.

      using Real = typename Type::Real;              ///< See LhsExpressionInterface.
      using Identifier = typename Type::Identifier;  ///< See LhsExpressionInterface.

      using Index = typename Matrix::Index;  ///< Index of an Eigen matrix.

    private:

      /// Clone a sparse matrix.
      template<typename R, typename T, typename M>
      R* cloneMatrix(M* mat) {
        R* r = new R(mat->rows(), mat->cols());

        std::vector<Eigen::Triplet<T>> entries(mat->nonZeros());

        Index outerSize = mat->outerSize();

        for (int k = 0; k < outerSize; ++k) {
          typename M::InnerIterator iter(*mat, k);
          while (iter) {
            entries.push_back(Eigen::Triplet<T>(iter.row(), iter.col(), T()));
            ++iter;
          }
        }
        r->setFromTriplets(entries.begin(), entries.end());

        return r;
      }

    public:

      /*******************************************************************************/
      /// @name Mandatory: Implementations for matrix and vector creation and deletion.
      /// @{

      /// \copydoc LinearSystemInterface::createMatrixReal
      template<typename M>
      MatrixReal* createMatrixReal(M* mat) {
        return cloneMatrix<MatrixReal, Real>(mat);
      }

      /// \copydoc LinearSystemInterface::createMatrixIdentifier
      template<typename M>
      MatrixIdentifier* createMatrixIdentifier(M* mat) {
        return cloneMatrix<MatrixIdentifier, Identifier>(mat);
      }

      /// @}
      /*******************************************************************************/
      /// @name Mandatory: Implementations for matrix and vector iterations.
      /// @{

      /// \copydoc LinearSystemInterface::iterateMatrix
      template<typename Func, typename MatrixA, typename MatrixB>
      void iterateMatrix(Func func, MatrixA* matA, MatrixB* matB) {
        Index outerSize = matA->outerSize();

        for (int k = 0; k < outerSize; ++k) {
          typename MatrixA::InnerIterator iterA(*matA, k);
          typename MatrixB::InnerIterator iterB(*matB, k);
          while (iterA) {
            func(iterA.valueRef(), iterB.valueRef());

            ++iterA;
            ++iterB;
          }
        }
      }

      /// \copydoc LinearSystemInterface::iterateMatrix
      template<typename Func, typename MatrixA, typename MatrixB, typename MatrixC>
      void iterateMatrix(Func func, MatrixA* matA, MatrixB* matB, MatrixC* matC) {
        Index outerSize = matA->outerSize();

        for (int k = 0; k < outerSize; ++k) {
          typename MatrixA::InnerIterator iterA(*matA, k);
          typename MatrixB::InnerIterator iterB(*matB, k);
          typename MatrixC::InnerIterator iterC(*matC, k);
          while (iterA) {
            func(iterA.valueRef(), iterB.valueRef(), iterC.valueRef());

            ++iterA;
            ++iterB;
            ++iterC;
          }
        }
      }

      /// @}
      /*******************************************************************************/
      /// @name Implementation for reverse mode AD.
      /// @{

      /// \copydoc LinearSystemInterface::iterateDyadic
      template<typename Func>
      void iterateDyadic(Func func, MatrixIdentifier* mat_id, VectorReal* x_v, VectorReal* b_b) {
        Index outerSize = mat_id->outerSize();

        for (int k = 0; k < outerSize; ++k) {
          typename MatrixIdentifier::InnerIterator iterId(*mat_id, k);
          while (iterId) {
            func(iterId.valueRef(), x_v->coeffRef(iterId.col()), b_b->coeffRef(iterId.row()));

            ++iterId;
          }
        }
      }

      /// @}
  };
}

#endif
