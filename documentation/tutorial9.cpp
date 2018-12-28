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

codi::RealReverse func(const codi::RealReverse& a, const codi::RealReverse& b) {
  return a * b;
}

void call(bool clear, bool stats) {
  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  codi::RealReverse a = 4.0;
  codi::RealReverse b = 3.0;

  // record with respect to a
  tape.setActive();

  tape.registerInput(a);
  codi::RealReverse y = func(a, b);
  tape.registerOutput(y);

  tape.setPassive();
  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "f(4.0, 3.0) = " << y << std::endl;
  std::cout << "df/da(4.0, 3.0) = " << a.getGradient() << std::endl;
  std::cout << "df/db(4.0, 3.0) = " << b.getGradient() << std::endl;

  if(stats) {
    tape.printStatistics();
  }
  if(clear) {
    tape.deactivateValue(a);
  }

  // record with respect to b
  tape.reset();
  tape.setActive();

  tape.registerInput(b);
  y = func(a, b);
  tape.registerOutput(y);

  tape.setPassive();
  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "f(4.0, 3.0) = " << y << std::endl;
  std::cout << "df/da(4.0, 3.0) = " << a.getGradient() << std::endl;
  std::cout << "df/db(4.0, 3.0) = " << b.getGradient() << std::endl;

  if(stats) {
    tape.printStatistics();
  }
}
int main(int nargs, char** args) {

  std::cout << "Recording tapes without clear:" << std::endl;
  call(false, false);

  std::cout << std::endl;
  std::cout << "Recording tapes with clear:" << std::endl;
  call(true, false);

  return 0;
}
