
#include <codi.hpp>
#include <iostream>

//! [Example 4: Higher order derivatives compile time access]
using t1s = codi::RealForwardGen<double>;
using t2s = codi::RealForwardGen<t1s>;
using t3s = codi::RealForwardGen<t2s>;
using t4s = codi::RealForwardGen<t3s>;
using t5s = codi::RealForwardGen<t4s>;
using t6s = codi::RealForwardGen<t5s>;

using r6s = codi::RealReverseGen<t5s>;

using t1v = codi::RealForwardGen<double, codi::Direction<double, 2>>;
using t2v = codi::RealForwardGen<t1v>;

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
    DH::setAllDerivatives<1>(aFor, 1.0);

    t2s cFor = func(aFor);

    std::cout << "t0s:   " << DH::derivative<0, 0>(cFor) << std::endl;
    std::cout << "t1_1s: " << DH::derivative<1, 0>(cFor) << std::endl;
    std::cout << "t1_2s: " << DH::derivative<1, 1>(cFor) << std::endl;
    std::cout << "t2s:   " << DH::derivative<2, 0>(cFor) << std::endl;
  }

  {
    using DH = codi::DerivativeAccess<t6s>;

    t6s aFor = 2.0;

    // set all first order directions in order to get the 6. order derivative
    DH::setAllDerivatives<1>(aFor, 1.0);

    t6s cFor = func(aFor);

    std::cout << "t0s: " << cFor << std::endl;
    std::cout << "t6s: " << DH::derivative<6, 0>(cFor) << std::endl;
  }

  {
    using DH = codi::DerivativeAccess<r6s>;

    r6s::Tape& tape = r6s::getTape();
    r6s aRev = 2.0;
    // set all first order directions on the primal value
    DH::setAllDerivativesForward<1>(aRev, 1.0);

    tape.setActive();
    tape.registerInput(aRev);

    r6s cRev = func(aRev);

    tape.registerOutput(cRev);
    // set all first order directions on the adjoint value
    DH::setAllDerivativesReverse<1>(cRev, 1.0);

    tape.setPassive();
    tape.evaluate();

    std::cout << "r0s: " << cRev << std::endl;
    std::cout << "r6s: " << DH::derivative<6, 0>(aRev) << std::endl;

    tape.reset();
  }

  {
    using DH = codi::DerivativeAccess<t2v>;

    t2v aFor = 2.0;
    // set all first order directions in order to get the 2. order derivative
    DH::derivative<1, 0>(aFor) = {1.0, 2.0};
    DH::derivative<1, 1>(aFor) = 1.0;

    t2v cFor = func(aFor);

    std::cout << "t0v:   " << DH::derivative<0, 0>(cFor) << std::endl;
    std::cout << "t1_1v: " << DH::derivative<1, 0>(cFor) << std::endl;
    std::cout << "t1_2v: " << DH::derivative<1, 1>(cFor) << std::endl;
    std::cout << "t2v:   " << DH::derivative<2, 0>(cFor) << std::endl;
  }

  return 0;
}
//! [Example 4: Higher order derivatives compile time access]
