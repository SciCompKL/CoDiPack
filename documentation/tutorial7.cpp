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
    typedef codi::DerivativeHelper<t2s> DH;

    t2s aFor = 2.0;
    // set all first order directions in order to get the 2. order derivative
    DH::setDerivatives(aFor, 1, 1.0);

    t2s cFor = func(aFor);

    cout << "t0s:   " << DH::derivative(cFor, 0, 0) << std::endl;
    cout << "t1_1s: " << DH::derivative(cFor, 1, 0) << std::endl;
    cout << "t1_2s: " << DH::derivative(cFor, 1, 1) << std::endl;
    cout << "t2s:   " << DH::derivative(cFor, 2, 0) << std::endl;
  }

  {
    t6s aFor = 2.0;

    // set all first order directions in order to get the 6. order derivative
    typedef codi::DerivativeHelper<t6s> DH;
    DH::setDerivatives(aFor, 1, 1.0);

    t6s cFor = func(aFor);

    cout << "t0s: " << cFor << std::endl;
    cout << "t6s: " << DH::derivative(cFor, 6, 0) << std::endl;
  }

  {
    typedef codi::DerivativeHelper<r6s> DH;

    r6s::TapeType& tape = r6s::getGlobalTape();
    r6s aRev = 2.0;
    // set all first order directions on the primal value
    DH::setDerivativesForward(aRev, 1, 1.0);

    tape.setActive();
    tape.registerInput(aRev);

    r6s cRev = func(aRev);

    tape.registerOutput(cRev);
    // set all first order directions on the adjoint value
    DH::setDerivativesReverse(cRev, 1, 1.0);

    tape.setPassive();
    tape.evaluate();

    cout << "r0s: " << cRev << std::endl;
    cout << "r6s: " << DH::derivative(aRev, 6, 0) << std::endl;
  }

  return 0;
}
