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

struct BaseLinearSystemSolverHandler {
  public:

    template<typename Number>
    static void solveSystemDirect(Number* A, Number* b, Number* x) {
      Number det = 1.0 / (A[0] * A[3] - A[1] * A[2]);

      x[0] = (A[3] * b[0] - A[1] * b[1]) * det;
      x[1] = (-A[2] * b[0] + A[0] * b[1]) * det;
    }

    template<typename Solver, typename M, typename V, typename Number>
    static void func(Solver solver, M& A, V& b, V& sol, Number scale) {
      for (int i = 0; i < 2; i += 1) {
        sol[i] = scale * b[i];
      }

      for (int i = 0; i < 2; i += 1) {
        b[i] = sol[i] / scale;
      }

#if CODI_EnableEigen
      codi::solveLinearSystem(solver, A, b, sol);
#else
      (void)solver;
      solveSystemDirect(A, b, sol);
#endif
    }
};
