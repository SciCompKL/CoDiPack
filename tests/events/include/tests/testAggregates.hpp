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

  std::cout << "# Active type computations" << std::endl;

  if(nInputs < 4) {
    std::cerr << "Test requires at least 4 inputs." << std::endl;
  }
  if(nOutputs < 4) {
    std::cerr << "Test requires at least 4 outputs." << std::endl;
  }

  std::complex<ActiveType> a(inputs[0], inputs[1]);
  std::complex<ActiveType> b(inputs[0], inputs[1]);

  // active type

  std::complex<ActiveType> x = 0.0;    // passive
  std::complex<ActiveType> y = a;      // copy
  std::complex<ActiveType> z = a * b;  // expression

  x = a * b;  // expression
  y = b;      // copy
  z = 2.0;    // passive

  std::cout << "# Assign outputs" << std::endl;
  outputs[0] = a.real();
  outputs[1] = a.imag();
  outputs[2] = b.real();
  outputs[3] = b.imag();
}
