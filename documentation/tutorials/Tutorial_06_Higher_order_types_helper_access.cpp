
#include <codi.hpp>
#include <iostream>

//! [Tutorial 6 - Higher order derivatives]
using t1s = codi::RealForwardGen<double>;
using t2s = codi::RealForwardGen<t1s>;
using t3s = codi::RealForwardGen<t2s>;
using t4s = codi::RealForwardGen<t3s>;
using t5s = codi::RealForwardGen<t4s>;
using t6s = codi::RealForwardGen<t5s>;

using r6s = codi::RealReverseGen<t5s>;

//! [Function]
template<typename T>
T func(const T& x) {
  T t = x * x * x * x * x * x * x;
  return t * 3.0;
}
//! [Function]

int main() {

  {
    using DH = codi::DerivativeAccess<t2s>;

    t2s aFor = 2.0;
    // set all first order directions in order to get the 2. order derivative
    DH::setAllDerivatives(aFor, 1, 1.0);

    t2s cFor = func(aFor);

    std::cout << "t0s:   " << DH::derivative(cFor, 0, 0) << std::endl;
    std::cout << "t1_1s: " << DH::derivative(cFor, 1, 0) << std::endl;
    std::cout << "t1_2s: " << DH::derivative(cFor, 1, 1) << std::endl;
    std::cout << "t2s:   " << DH::derivative(cFor, 2, 0) << std::endl;
  }

  {
    using DH = codi::DerivativeAccess<t6s>;

    t6s aFor = 2.0;

    // set all first order directions in order to get the 6. order derivative
    DH::setAllDerivatives(aFor, 1, 1.0);

    t6s cFor = func(aFor);

    std::cout << "t0s: " << cFor << std::endl;
    std::cout << "t6s: " << DH::derivative(cFor, 6, 0) << std::endl;
  }

  {
    using DH = codi::DerivativeAccess<r6s>;

    r6s::Tape& tape = r6s::getTape();
    r6s aRev = 2.0;
    // set all first order directions on the primal value
    DH::setAllDerivativesForward(aRev, 1, 1.0);

    tape.setActive();
    tape.registerInput(aRev);

    r6s cRev = func(aRev);

    tape.registerOutput(cRev);
    // set all first order directions on the adjoint value
    DH::setAllDerivativesReverse(cRev, 1, 1.0);

    tape.setPassive();
    tape.evaluate();

    std::cout << "r0s: " << cRev << std::endl;
    std::cout << "r6s: " << DH::derivative(aRev, 6, 0) << std::endl;

    tape.reset();
  }

  return 0;
}
//! [Tutorial 6 - Higher order derivatives]
