/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
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
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *     Max Sagebaum
 *     Tim Albring
 *     Johannes Blühdorn
 */
#include <codi.hpp>
#include <iostream>

codi::RealReverse func(const codi::RealReverse& x) {
  return x * x * x;
}

void call(bool reset, bool stats) {
  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();

  const size_t n = 5;
  double points[n] = {2.0, 2.1, 2.5, 3.0, -1.0};

  for(size_t i = 0; i < n; i += 1) {
    if(reset) {
      tape.reset();
    }
    codi::RealReverse x = points[i];

    tape.setActive();

    tape.registerInput(x);
    codi::RealReverse y = func(x);
    tape.registerOutput(y);

    tape.setPassive();
    y.setGradient(1.0);
    tape.evaluate();

    std::cout << "f(" << x.value() << ") = " << y << std::endl;
    std::cout << "df/dx(" << x.value() << ") = " << x.getGradient() << std::endl;

    if(stats) {
      tape.printStatistics();
    }
  }
}
int main(int nargs, char** args) {

  std::cout << "Recording tapes without reset:" << std::endl;
  call(false, true);

  std::cout << std::endl;
  std::cout << "Recording tapes without reset:" << std::endl;
  call(true, true);

  return 0;
}
