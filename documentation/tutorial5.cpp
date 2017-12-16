/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2017 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
#include <sstream>

void func(const codi::RealReverse& x, codi::RealReverse& y) {
  y = 3.0*x*x*x*x + 5.0*x*x*x - 3.0*x*x + 2.0*x -4.0;
}

typedef codi::ReferenceActiveReal<codi::RealReverse> RefReal;

void funcRef(const codi::RealReverse& x, codi::RealReverse& y) {
  RefReal xRef = x;

  y = 3.0*xRef*xRef*xRef*xRef + 5.0*xRef*xRef*xRef - 3.0*xRef*xRef + 2.0*xRef -4.0;
}

int main(int nargs, char** args) {
  codi::RealReverse x = 3.14;
  codi::RealReverse y;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();

  std::cout << "Func with standard codi type." << std::endl;
  tape.setActive();

  tape.registerInput(x);
  func(x, y);
  tape.registerOutput(y);

  tape.setPassive();
  std::cout << "f(3.14) = (" << y << ")" << std::endl;

  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "df/dx = (" << x.getGradient() << ")" << std::endl;

  std::ostringstream run1;
  tape.printStatistics(run1);

  tape.reset();

  std::cout << "Func with reference codi type." << std::endl;
  tape.setActive();

  tape.registerInput(x);
  funcRef(x, y);
  tape.registerOutput(y);

  tape.setPassive();
  std::cout << "f(3.14) = (" << y << ")" << std::endl;

  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "df/dx = (" << x.getGradient() << ")" << std::endl;

  std::ostringstream run2;
  tape.printStatistics(run2);

  tape.reset();

  std::cout << std::endl;
  std::cout << "Statistics for the standard codi type:" << std::endl;
  std::cout << run1.str() << std::endl << std::endl;

  std::cout << "Statistics for the reference codi type:" << std::endl;
  std::cout << run2.str() << std::endl << std::endl;
}
