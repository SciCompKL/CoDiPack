Tutorial 7.2: Higher order derivatives without a helper interface {#Tutorial7_2}
============

This is the same tutorial as tutorial [Tutorial 7](@ref Tutorial7) but it will use no helper class to set the derivatives.

The helper class performs some kind of branch selection algorithm on the higher order AD structures.
The user needs to select the branches manually, for which the functionality were available since the first release of CoDiPack.
The code will use the two functions [value](@ref codi::ActiveReal::value) and [gradient](@ref codi::ActiveReal::gradient) of the ActiveReal class in order to select the primal order gradient component of the CoDiPack type.

For a second order type, there are four different possibilites for the selection process:
~~~~{.cpp}
  t2s aFor2 = 2.0;

  aFor2.value().value() = 1.0;       // <-- select the primal value of the type
  aFor2.value().gradient() = 2.0;    // <-- select the first first order derivative
  aFor2.gradient().value() = 2.0;    // <-- select the second first order derivative
  aFor2.gradient().gradient() = 4.0; // <-- select the second order derivative
~~~~

If the "value" branch is selected, then the derivative order does not increase.
On the other hand, if the "gradient" branch is selected, then the derivative order increases by one, therefore in order to select the second order derivative, we have to go two times the gradient path.

For higher order derivatives, the selection becomes quite involved and therefore the [DerivativeHelper](@ref codi::DerivativeHelper) class was written.
Now, for the standard second order example, the code is written as:
~~~~{.cpp}
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
~~~~

It yields the same output as in [Tutorial 7](@ref Tutorial7).

The sixth order derivative example receives now the same changes and we need to set all six sixth order derivatives manually:
~~~~{.cpp}
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
~~~~

The "gradient().gradient().gradient()" etc. call in order to get the sixth order derivative, is quite cumbersome to write and not very intuitive for the reader of the code. It is therefore advisable to use the DerivativeHelper structure whenever possible.

The same changes can be made to the reverse example.
It is now very important, that the path where the gradient is choosen first, is just set before the adjoint evaluation:
~~~~{.cpp}
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
~~~~

Both examples will compute the same values as in the previous tutorial.

The full code of the tutorial is:
~~~~{.cpp}
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
~~~~
