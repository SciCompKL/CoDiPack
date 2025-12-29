/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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

#include <codi/config.h>
#include <codi/misc/macros.hpp>
#include <codi/tools/lowlevelFunctions/generationHelperCoDiPack.hpp>
#include <codi/tools/lowlevelFunctions/eigenWrappers.hpp>
#include <codi/tools/lowlevelFunctions/lowLevelFunctionCreationUtilities.hpp>

/** \copydoc codi::Namespace */
namespace codi {

  /**
   *  Low level function for \f$R = A * B\f$ with
   *   - \f$ R \in \R^{n \times m} \f$
   *   - \f$ A \in \R^{n \times k} \f$
   *   - \f$ B \in \R^{k \times m} \f$
   *
   * @tparam eigenStore One of Eigen::StorageOptions.
   */
  template<Eigen::StorageOptions eigenStore, typename Type, typename ActiveType = Type>
  struct matrixMatrixMultiplication : public LowLevelFunction {
      [[in, store(active_B), size(n * k)]] Type* ACTIVE_ARG(A);
      [[in, store(active_A), size(k * m)]] Type* ACTIVE_ARG(B);
      [[out, size(n * m)]] Type* ACTIVE_ARG(R);

      [[storeType(uint8_t), order(1)]] int n;
      [[storeType(uint8_t), order(2)]] int k;
      [[storeType(uint8_t), order(3)]] int m;

      void primal() {
        mapEigen<eigenStore>(R, n, m) = mapEigen<eigenStore>(A, n, k) * mapEigen<eigenStore>(B, k, m);
      }

      void primal_activity() {
        if(active) {
          mapEigen<eigenStore>(R_i_out, n, m).setZero();
          if(active_A) {
            mapEigen<eigenStore>(R_i_out, n, m).colwise() += mapEigen<eigenStore>(A_i_in, n, k).rowwise().any();
          }
          if(active_B) {
            mapEigen<eigenStore>(R_i_out, n, m).rowwise() += mapEigen<eigenStore>(B_i_in, k, m).colwise().any();
          }
        }
      }

      void diff_A_fwd() {
        mapEigen<eigenStore>(R_d_out, n, m) += mapEigen<eigenStore>(A_d_in, n, k) * mapEigen<eigenStore>(B, k, m);
      }
      void diff_B_fwd() {
        mapEigen<eigenStore>(R_d_out, n, m) += mapEigen<eigenStore>(A, n, k) * mapEigen<eigenStore>(B_d_in, k, m);
      }

      void diff_A_rws() {
        mapEigen<eigenStore>(A_b_in, n, k) = mapEigen<eigenStore>(R_b_out, n, m) * mapEigen<eigenStore>(B, k, m).transpose();
      }

      void diff_B_rws() {
        mapEigen<eigenStore>(B_b_in, k, m) = mapEigen<eigenStore>(A, n, k).transpose() * mapEigen<eigenStore>(R_b_out, n, m);
      }
  };

    /**
  *  Low level function for \f$R = A * B\f$ with
  *   - \f$ R \in \R^{n \times m} \f$
  *   - \f$ A \in \R^{n \times k} \f$
  *   - \f$ B \in \R^{k \times m} \f$
  */
  template<typename Type>
  void matrixMatrixMultiplicationRowMajor(
          Type const* A,
          Type const* B,
          Type* R,
          int n,
          int k,
          int m
  ) {
    matrixMatrixMultiplication<Eigen::StorageOptions::RowMajor>(A, B, R, n, k, m);
  }

  /**
  *  Low level function for \f$R = A * B\f$ with
  *   - \f$ R \in \R^{n \times m} \f$
  *   - \f$ A \in \R^{n \times k} \f$
  *   - \f$ B \in \R^{k \times m} \f$
  */
  template<typename Type>
  void matrixMatrixMultiplicationColMajor(
          Type const* A,
          Type const* B,
          Type* R,
          int n,
          int k,
          int m
  ) {
    matrixMatrixMultiplication<Eigen::StorageOptions::ColMajor>(A, B, R, n, k, m);
  }
}
