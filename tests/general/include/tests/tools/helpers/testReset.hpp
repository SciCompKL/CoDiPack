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
#include "../../../testInterface.hpp"

struct TestReset : public TestInterface {
  public:
    NAME("Reset")
    IN(2)
    OUT(1)
    POINTS(1) = {{1.0, 0.5}};

    template<typename Number>
    static void func(Number* x, Number* y) {
#if REVERSE_TAPE
      typename Number::Tape& tape = Number::getTape();
#endif

      Number a = x[0] * x[1];
      Number b = x[0] / sin(x[1]);
      Number c = b * a;

#if REVERSE_TAPE
      typename Number::Tape::Position pos = tape.getPosition();
#endif

      b = a * x[0];

#if REVERSE_TAPE
      tape.resetTo(pos);
#endif

      y[0] = c * a;
    }
};
