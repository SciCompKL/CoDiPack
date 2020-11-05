/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2020 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

void func(const codi::RealReverse* x, size_t l, codi::RealReverse* y) {
  y[0] = 0.0;
   y[1] = 1.0;
  for(size_t i = 0; i < l; ++i) {
    y[0] += x[i];
    y[1] *= x[i];
  }
}

int main(int nargs, char** args) {
  codi::RealReverse x[5];
  codi::RealReverse y[2];
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  x[3] = 4.0;
  x[4] = 5.0;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setActive();

  for(size_t i = 0; i < 5; ++i) {
    tape.registerInput(x[i]);
  }
  func(x, 5, y);
  tape.registerOutput(y[0]);
  tape.registerOutput(y[1]);

  tape.setPassive();
  std::cout << "f(1 .. 5) = (" << y[0] << ", " << y[1] << ")" << std::endl;

  y[0].setGradient(1.0);
  tape.evaluate();

  std::cout << "df_1/dx(1 .. 5) = (";
  for(size_t i = 0; i < 5; ++i) {
    if(0 != i) {
      std::cout << ", ";
    }
    std::cout << x[i].getGradient();
  }
  std::cout << ")" << std::endl;

  tape.clearAdjoints();
  y[1].setGradient(1.0);
  tape.evaluate();

  std::cout << "df_2/dx(1 .. 5) = (";
  for(size_t i = 0; i < 5; ++i) {
    if(0 != i) {
      std::cout << ", ";
    }
    std::cout << x[i].getGradient();
  }
  std::cout << ")" << std::endl;

  return 0;
}
