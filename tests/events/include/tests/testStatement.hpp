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

  // active type

  ActiveType x = 0.0;    // passive
  ActiveType y = a;      // copy
  ActiveType z = a * b;  // expression

  x = c * d;  // expression
  y = b;      // copy
  z = 2.0;    // passive

  std::cout << "# Assign outputs" << std::endl;
  for (size_t i = 0; i < nInputs; ++i) {
    outputs[i] = sin(i * (x + y)) + cos(y * z / (i + 1));
  }

  // active type wrapper

  std::cout << "# Active type wrapper computations" << std::endl;
  codi::ActiveTypeWrapper<ActiveType> aRef(a.value(), a.getIdentifier());
  codi::ActiveTypeWrapper<ActiveType> bRef(b.value(), b.getIdentifier());

  codi::ActiveTypeWrapper<ActiveType> xRef(x.value(), x.getIdentifier());
  codi::ActiveTypeWrapper<ActiveType> yRef(y.value(), y.getIdentifier());
  codi::ActiveTypeWrapper<ActiveType> zRef(z.value(), z.getIdentifier());

  xRef = 0.0;          // passive
  yRef = aRef;         // copy
  zRef = aRef * bRef;  // expression

  std::cout << "# Assign outputs" << std::endl;
  for (size_t i = 0; i < nInputs; ++i) {
    outputs[i] = sin(i * (xRef + yRef)) + cos(yRef * zRef / (i + 1));
  }

  // active type copy
  std::cout << "# Immutable active type computations" << std::endl;
  codi::ImmutableActiveType<ActiveType> cCopy(c.value(), c.getIdentifier());
  codi::ImmutableActiveType<ActiveType> dCopy(d.value(), d.getIdentifier());

  x = cCopy * dCopy;  // expression
  y = cCopy;          // copy

  std::cout << "# Assign outputs" << std::endl;
  for (size_t i = 0; i < nInputs; ++i) {
    outputs[i] = sin(i * (x + y));
  }
}
