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

template<typename ActiveType>
void test(size_t nInputs, ActiveType* inputs, size_t nOutputs, ActiveType* outputs) {
  // process inputs

  std::cout << "# Active type computations" << std::endl;
  ActiveType a = 0.0, b = 0.0, c = 0.0, d = 0.0;

  for (size_t i = 0; i < nInputs; ++i) {
    a += sin(inputs[i]);
    b += cos(inputs[i]);
    c += 3.0 * inputs[i];
    d += inputs[i] * inputs[i];
  }

  // computations

  std::cout << "# Preaccumulation" << std::endl;
  codi::PreaccumulationHelper<ActiveType> ph;

  // forward mode preacc

  ph.start(a);

  ActiveType v = a * a;
  ActiveType w = v * cos(a);
  ActiveType x = v * w;
  ActiveType y = b + b * w;

  ph.finish(false, v, w, x, y);

  // reverse mode preacc

  ph.start(a, b, c, d);

  ActiveType z = a * b + c * d;

  ph.finish(false, z);

  // produce outputs

  std::cout << "# Assign outputs" << std::endl;
  for (size_t i = 0; i < nOutputs; ++i) {
    outputs[i] = exp(v * w / (i + 1)) + sin(i * (x + y)) + cos(z / (i + 1));
  }
}
