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

#include "../../../testInterface.hpp"
#include "baseLinearSystemSolverHandler.hpp"

struct TestEigenLinearSystemSolverHandler : public TestInterface {
  public:
    NAME("EigenLinearSystemSolverHandler")
    IN(6)
    OUT(2)
    POINTS(1) = {{1.0, 2.0, 3.0, 4.0, 20.0, 10.0}};

#if CODI_EnableEigen
    template<typename T>
    using Matrix = Eigen::Matrix<T, 2, 2>;
    template<typename T>
    using Vector = Eigen::Matrix<T, 2, 1>;

    template<typename Number>
    struct EigenLinearSystemTest : public codi::EigenLinearSystem<Number, Matrix, Vector> {
      public:

        using Base = codi::EigenLinearSystem<Number, Matrix, Vector>;
        using MatrixReal = typename Base::MatrixReal;
        using VectorReal = typename Base::VectorReal;

        void solveSystem(MatrixReal const* A, VectorReal const* b, VectorReal* x) {
          *x = A->colPivHouseholderQr().solve(*b);
        }
    };
#endif

    template<typename Number>
    static void func(Number* x, Number* y) {
#if CODI_EnableEigen
      Matrix<Number> A;
      A << x[0], x[1], x[2], x[3];
      Vector<Number> b = {x[4], x[5]};
      Vector<Number> sol;

      using Solver = EigenLinearSystemTest<Number>;
#else
      Number* A = &x[0];
      Number b[2] = {x[4], x[5]};
      Number sol[2];

      using Solver = int;
#endif

      BaseLinearSystemSolverHandler::func(Solver(), A, b, sol, b[0]);

      y[0] = sol[0];
      y[1] = sol[1];
    }
};
