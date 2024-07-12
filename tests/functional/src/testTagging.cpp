/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2024 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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
#define CODI_IDE 0

#include <codi.hpp>
#include <fstream>
#include <iostream>

using Real = codi::RealReverseTag;
using Tape = typename Real::Tape;

Real func(const Real& x, const Real& y) {
  return x * y;
}

static void tagPropertyErrorCallback(double const& currentValue, double const& newValue, codi::TagFlags flag, void* userData) {
  std::ofstream* out = (std::ofstream*)userData;

  *out << "Property error '" << std::to_string(flag) << "' on value. current value: " << currentValue << " new value: " << newValue << "" << std::endl;
}

static void tagErrorCallback(int const& correctTag, int const& wrongTag, void* userData) {
  std::ofstream* out = (std::ofstream*)userData;

  *out << "Use of variable with bad tag '" << wrongTag << "', should be '" << correctTag << "'." << std::endl;
}


int main(int nargs, char** args) {

  std::ofstream out("run.out");

  Real x = 4.0;
  Real y = 3.0;
  Real z = 1.0;

  codi::PreaccumulationHelper<Real> ph;
  Tape& tape = Real::getTape();
  tape.setTagErrorCallback(tagErrorCallback, &out);
  tape.setTagPropertyErrorCallback(tagPropertyErrorCallback, &out);
  tape.setCurTag(42);
  tape.setActive();

  tape.registerInput(x);
  tape.registerInput(y);

  out << "Default test:" << std::endl;
  ph.start(x, y);
  Real w = func(x, y);
  ph.finish(false, w);
  w = w * z;

  out << "Input error test:" << std::endl;
  ph.start(x);
  w = func(x, y);
  ph.finish(false, w);
  w = w * z;

  out << "Output error test:" << std::endl;
  ph.start(x, y);
  w = func(x, y);
  ph.finish(false);
  w = w * z;

  out << "Do not use error:" << std::endl;
  tape.setTagPropertyOnVariable(x, codi::TagFlags::DoNotUse);
  w = func(x, y);
  tape.clearTagPropertiesOnVariable(x);

  out << "Do not change with same value:" << std::endl;
  tape.setTagPropertyOnVariable(w, codi::TagFlags::DoNotChange);
  w = func(x, y);

  out << "Do not change error test:" << std::endl;
  tape.setTagPropertyOnVariable(w, codi::TagFlags::DoNotChange);
  w = func(x, z);
  tape.clearTagPropertiesOnVariable(w);

  out << "Do not write error test:" << std::endl;
  tape.setTagPropertyOnVariable(w, codi::TagFlags::DoNotWrite);
  w = func(x, z);

  tape.registerOutput(y);

  tape.setPassive();
  tape.reset();

  out.close();

  return 0;
}
