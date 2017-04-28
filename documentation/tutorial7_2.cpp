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

#include <iostream>

#include <codi.hpp>

using namespace std;

typedef codi::RealForwardGen<double> t1s;
typedef codi::RealForwardGen<t1s>    t2s;
typedef codi::RealForwardGen<t2s>    t3s;
typedef codi::RealForwardGen<t3s>    t4s;
typedef codi::RealForwardGen<t4s>    t5s;
typedef codi::RealForwardGen<t5s>    t6s;

typedef codi::RealReverseGen<t5s>    r6s;

template<typename T>
T func(const T& x) {
  T t = x * x * x * x * x * x * x;
  return t * 3.0;
}

int main() {

  {
    t2s aFor2 = 2.0;
    // set all first order directions in order to get the 2. order derivative
    aFor2.value().gradient() = 1.0;
    aFor2.gradient().value() = 1.0;

    t2s cFor2 = func(aFor2);

    cout << "t0s:   " << cFor2.value().value() << std::endl;
    cout << "t1_1s: " << cFor2.value().gradient() << std::endl;
    cout << "t1_2s: " << cFor2.gradient().value() << std::endl;
    cout << "t2s:   " << cFor2.gradient().gradient() << std::endl;
  }

  {
    t6s aFor = 2.0;

    // set all first order directions in order to get the 6. order derivative
    aFor.value().value().value().value().value().gradient() = 1.0;
    aFor.value().value().value().value().gradient().value() = 1.0;
    aFor.value().value().value().gradient().value().value() = 1.0;
    aFor.value().value().gradient().value().value().value() = 1.0;
    aFor.value().gradient().value().value().value().value() = 1.0;
    aFor.gradient().value().value().value().value().value() = 1.0;

    t6s cFor = func(aFor);

    cout << "t0s: " << cFor << std::endl;
    cout << "t6s: " << cFor.gradient().gradient().gradient().gradient().gradient().gradient() << std::endl;
  }

  {
    r6s::TapeType& tape = r6s::getGlobalTape();
    r6s aRev = 2.0;
    // set all first order directions on the primal value
    aRev.value().value().value().value().value().gradient() = 1.0;
    aRev.value().value().value().value().gradient().value() = 1.0;
    aRev.value().value().value().gradient().value().value() = 1.0;
    aRev.value().value().gradient().value().value().value() = 1.0;
    aRev.value().gradient().value().value().value().value() = 1.0;

    tape.setActive();
    tape.registerInput(aRev);

    r6s cRev = func(aRev);

    tape.registerOutput(cRev);
    // set all first order directions on the adjoint value
    cRev.gradient().value().value().value().value().value() = 1.0;

    tape.setPassive();
    tape.evaluate();

    cout << "r0s: " << cRev << std::endl;
    cout << "r6s: " << aRev.gradient().gradient().gradient().gradient().gradient().gradient() << std::endl;
  }

  return 0;
}
