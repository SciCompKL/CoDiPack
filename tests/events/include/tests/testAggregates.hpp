/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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

#include <codi.hpp>

#if CODI_SpecializeStdComplex
template<typename T>
using TestComplex = std::complex<T>;
#else
template<typename T>
using TestComplex = codi::ActiveComplex<T>;
#endif

template<typename ActiveType>
void test(size_t nInputs, ActiveType* inputs, size_t nOutputs, ActiveType* outputs) {
  // process inputs

  std::cout << "# Active type computations" << std::endl;

  if (nInputs < 4) {
    std::cerr << "Test requires at least 4 inputs." << std::endl;
  }
  if (nOutputs < 4) {
    std::cerr << "Test requires at least 4 outputs." << std::endl;
  }

  TestComplex<ActiveType> a(inputs[0], inputs[1]);
  TestComplex<ActiveType> b(inputs[0], inputs[1]);

  // active type

  TestComplex<ActiveType> x = 0.0;    // passive
  TestComplex<ActiveType> y = a;      // copy
  TestComplex<ActiveType> z = a * b;  // expression

  x = a * b;  // expression
  y = b;      // copy
  z = 2.0;    // passive

  std::cout << "# Assign outputs" << std::endl;
  outputs[0] = a.real();
  outputs[1] = a.imag();
  outputs[2] = b.real();
  outputs[3] = b.imag();
}
