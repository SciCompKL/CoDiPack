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

#include <codi.hpp>

#include "../../../../testInterface.hpp"

struct TestMatrixMatrixMultiplication : public TestInterface {
  public:
    NAME("MatrixMatrixMultiplication")
    IN(1)
    OUT(1)
    POINTS(1) = {{2.0}};

    static int constexpr ITER = 5;

    template<typename Number>
    static void func(Number* x, Number* y) {
      Number A[4] = {1.0 * x[0], 2.0 * x[0], 3.0 * x[0], 4.0 * x[0]};
      Number B[4] = {0.4 * x[0], 0.3 * x[0], 0.2 * x[0], 0.1 * x[0]};
      Number C[4];

#if REVERSE_TAPE && CODI_EnableEigen
      codi::matrixMatrixMultiplicationRowMajor(A, B, C, 2, 2, 2);
#else
      for (int row = 0; row < 2; row += 1) {
        for (int col = 0; col < 2; col += 1) {
          C[row * 2 + col] = A[row * 2 + 0] * B[0 * 2 + col] + A[row * 2 + 1] * B[1 * 2 + col];
        }
      }
#endif

      y[0] = C[0] + C[1] + C[2] + C[3];
    }
};
