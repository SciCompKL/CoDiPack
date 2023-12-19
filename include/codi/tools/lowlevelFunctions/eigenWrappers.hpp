#pragma once
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

#if CODI_EnableEigen

  #include <Eigen/Eigen>

  #include "../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Abbreviation for a mapped Eigen matrix.
  template<typename T, Eigen::StorageOptions store>
  using MapEigenMatrix = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, store>>;

  /// Abbreviation for a constant mapped Eigen matrix.
  template<typename T, Eigen::StorageOptions store>
  using MapEigenMatrixConst = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, store> const>;

  /// Abbreviation for a mapped Eigen vector.
  template<typename T>
  using MapEigenVector = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>;

  /// Abbreviation for a constant mapped Eigen vector.
  template<typename T>
  using MapEigenVectorConst = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1> const>;

  /// Create a mapped Eigen matrix with specified storing option.
  template<Eigen::StorageOptions store, typename T>
  MapEigenMatrix<T, store> mapEigen(T* p, int rows, int cols) {
    return MapEigenMatrix<T, store>(p, rows, cols);
  }

  /// Create a mapped Eigen matrix with specified storing option.
  template<Eigen::StorageOptions store, typename T>
  MapEigenMatrixConst<T, store> mapEigen(T const* p, int rows, int cols) {
    return MapEigenMatrixConst<T, store>(p, rows, cols);
  }

  /// Create a mapped Eigen matrix with a row major data layout.
  template<typename T>
  MapEigenMatrix<T, Eigen::StorageOptions::RowMajor> mapEigenRowMajor(T* p, int rows, int cols) {
    return mapEigen<Eigen::StorageOptions::RowMajor>(p, rows, cols);
  }

  /// Create a constant mapped Eigen matrix with a row major data layout.
  template<typename T>
  MapEigenMatrixConst<T, Eigen::StorageOptions::RowMajor> mapEigenRowMajor(T const* p, int rows, int cols) {
    return mapEigen<Eigen::StorageOptions::RowMajor>(p, rows, cols);
  }

  /// Create a mapped Eigen matrix with a column major data layout.
  template<typename T>
  MapEigenMatrix<T, Eigen::StorageOptions::ColMajor> mapEigenColMajor(T* p, int rows, int cols) {
    return mapEigen<Eigen::StorageOptions::ColMajor>(p, rows, cols);
  }

  /// Create a constant mapped Eigen matrix with a column major data layout.
  template<typename T>
  MapEigenMatrixConst<T, Eigen::StorageOptions::ColMajor> mapEigenColMajor(T const* p, int rows, int cols) {
    return mapEigen<Eigen::StorageOptions::ColMajor>(p, rows, cols);
  }

  /// Create a mapped Eigen vector.
  template<typename T>
  MapEigenVector<T> mapEigen(T* p, int size) {
    return MapEigenVector<T>(p, size);
  }

  /// Create a constant mapped Eigen vector.
  template<typename T>
  MapEigenVectorConst<T> mapEigen(T const* p, int size) {
    return MapEigenVectorConst<T>(p, size);
  }
}

#endif
