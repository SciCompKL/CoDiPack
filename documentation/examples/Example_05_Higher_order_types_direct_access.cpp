
#include <codi.hpp>
#include <iostream>

//! [Example 5: Higher order derivatives direct access]
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
    t2s aFor = 2.0;
    // set all first order directions in order to get the 2. order derivative
    aFor.value().gradient() = 1.0;
    aFor.gradient().value() = 1.0;

    t2s cFor = func(aFor);

    std::cout << "t0s:   " << cFor.value().value() << std::endl;
    std::cout << "t1_1s: " << cFor.value().gradient() << std::endl;
    std::cout << "t1_2s: " << cFor.gradient().value() << std::endl;
    std::cout << "t2s:   " << cFor.gradient().gradient() << std::endl;
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

    std::cout << "t0s: " << cFor << std::endl;
    std::cout << "t6s: " << cFor.gradient().gradient().gradient().gradient().gradient().gradient() << std::endl;
  }

  {
    r6s::Tape& tape = r6s::getTape();
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

    std::cout << "r0s: " << cRev << std::endl;
    std::cout << "r6s: " << aRev.gradient().gradient().gradient().gradient().gradient().gradient() << std::endl;

    tape.reset();
  }

  return 0;
}
//! [Example 5: Higher order derivatives direct access]
