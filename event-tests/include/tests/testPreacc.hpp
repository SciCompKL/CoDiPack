/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 *  - SciComp, TU Kaiserslautern:
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

  ActiveType a = 0.0, b = 0.0, c = 0.0, d = 0.0;

  for (size_t i = 0; i < nInputs; ++i) {
    a += sin(inputs[i]);
    b += cos(inputs[i]);
    c += 3.0 * inputs[i];
    d += inputs[i] * inputs[i];
  }

  // computations

  // forward mode preacc

  codi::PreaccumulationHelper<ActiveType> ph;

  ph.start(a, b);

  ActiveType u = a * b;
  ActiveType v = u * cos(a);
  ActiveType w = u * v;
  ActiveType x = b + b * v;

  ph.finish(false, q, v, w, x);

  // reverse mode preacc

  ph.start(a, b, c, d);

  ActiveType y = a * b + c * d;
  ActiveType z = a * c + b * d;

  ph.finish(false, y, z);

  // produce outputs

  for (size_t i = 0; i < nOutputs; ++i) {
    outputs[i] = exp(u * v / (i + 1)) + sin(i * (w + x)) + cos(y * z / (i + 1));
  }

}
