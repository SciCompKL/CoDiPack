/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */
#include <codi.hpp>
#include <iostream>

void func(const codi::RealForward* x, size_t l, codi::RealForward* y) {
  y[0] = 0.0;
   y[1] = 1.0;
  for(size_t i = 0; i < l; ++i) {
    y[0] += x[i];
    y[1] *= x[i];
  }
}

int main(int nargs, char** args) {
  codi::RealForward x[5];
  codi::RealForward y[2];
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  x[3] = 4.0;
  x[4] = 5.0;

  for(size_t i = 0; i < 5; ++i) {
    x[i].setGradient(1.0);

    func(x, 5, y);

    if(0 == i) {
      std::cout << "f(1 .. 5) = (" << y[0] << ", " << y[1] << ")" << std::endl;
    }
    std::cout << "df/dx_" << (i + 1) << " (1 .. 5) = (" << y[0].getGradient() << ", " << y[1].getGradient() << ")" << std::endl;

    x[i].setGradient(0.0);
  }

  return 0;
}
