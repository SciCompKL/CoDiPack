Tutorial 7: Higher order derivatives {#Tutorial7}
============

In this tutorial we want to compute higher order derivatives with the help of CoDiPack.
This will be done by using the helper structure [DerivativeHelper](@ref codi::DerivativeHelper), which provides several functions to easily access higher order derivatives.

The function that we will differentiate is implemented as:
~~~~{.cpp}
    template<typename T>
    T func(const T& x) {
      T t = x * x * x * x * x * x * x;
      return t * 3.0;
    }
~~~~

The mathematical representation of the function is
\f[
  y = f(x) = 3*x^7 \eqdot
\f]
As the function is quite simple, all the derivatives up to the 6-th order can be computed by hand and they are:
\f[
  \frac{\d f}{\d x}(x) = 21 * x^6, \quad
  \frac{\d^2 f}{\d^2 x}(x) = 126 * x^5, \quad
  \frac{\d^3 f}{\d^3 x}(x) = 630 * x^4, \quad
  \frac{\d^4 f}{\d^4 x}(x) = 2520 * x^3, \quad
  \frac{\d^5 f}{\d^5 x}(x) = 7560 * x^2, \quad
  \frac{\d^6 f}{\d^6 x}(x) = 15120 * x \eqdot
\f]

Now that we have the gradient formulation of the functions, we can introduce how the same values can be comuted with higher order AD types.
If the genral formulation of real valued function is taken as \f$f: \R \rightarrow \R\f$ with \f$y = f(x)\f$, then the first application of the forward AD mode yields the computation of
\f{align*}{
  y =& f(x) \\
  \dot y^{(1)} =& \frac{\d f}{\d x}(x) \dot x^{(1)} \eqdot
\f}
The superscript \f$\cdot^{(1)}\f$ indicates the first application of the forward AD mode.
If it is now applied a second time we will use the superscript \f$\cdot^{(2)}\f$, to distinguish the first from the second time.
When a dot variable (e.g. \f$\dot x^{(1)}\f$) is differentiated a second time, the superscirpts are merged (e.g. \f$\dot x^{(1,2)}\f$).
The second application of the forward AD mode yields now
\f{align*}{
  y =& f(x) \\
  \dot y^{(1)} =& \frac{\d f}{\d x}(x) \dot x^{(1)} \\
  \dot y^{(2)} =& \frac{\d f}{\d x}(x) \dot x^{(2)} \\
  \dot y^{(1,2)} =& \frac{\d^2 f}{\d^2 x}(x) \dot x^{(1)} \dot x^{(2)} + \frac{\d f}{\d x}(x) \dot x^{(1,2)}\eqdot
\f}

From these equations we learn, that all first order tangent directions \f$\dot x^{(1)}\f$ and \f$\dot x^{(2)}\f$ need to be set in order to get the second order derivative \f$\frac{\d^2 f}{\d^2 x}(x)\f$.
The second order tangent direction \f$\dot x^{(1,2)}\f$ needs to be zero, otherwise the result will contain additional information.

The conclusion from this small example can be extended to an arbitrary derivative order.
When only the highest derivative information is required, then all the first order tangent directions need to be set to one.

The first step for our example is now to define the higher order derivative types.
Because we want to calculate 6-th order derivatives in this tutorial, we generate the types up this order.
The generation can now be done quite simple by using the generic "Gen" types of CoDiPack:

~~~~{.cpp}
typedef codi::RealForwardGen<double> t1s;
typedef codi::RealForwardGen<t1s>    t2s;
typedef codi::RealForwardGen<t2s>    t3s;
typedef codi::RealForwardGen<t3s>    t4s;
typedef codi::RealForwardGen<t4s>    t5s;
typedef codi::RealForwardGen<t5s>    t6s;

typedef codi::RealReverseGen<t5s>    r6s;
~~~~

This will insert the CoDiPack types into each other such that the forward AD mode is applied multiple times, which raises the question how many derivatives are available?
The standard CoDiPack type has two values, the primal and the derivative, if the type is used recursively, then there are \f$2^d\f$ values available where \f$d\f$ is the number of times the recursion is performed and the highest derivative order is then \f$d\f$.
For a sixth order type, this are already 64 values, that contain 6-th, 5-th, 4-th, 3-tr, 2-nd, 1-st order derivatives and the primal value.
In order to know how many derivatives of a certain type are available, the bonimal coefficient can be used:
\f[
  s = \binom{d}{n} \eqdot
\f]
\f$s\f$ is then the number of derivatives, that have the order \f$n\f$.
If we want to set the derivatives of a second order type manually, the following code would be necessary:
~~~~{.cpp}
  t2s aFor = 2.0;

  aFor.value().value() = 1.0;       // <-- select the primal value of the type
  aFor.value().gradient() = 2.0;    // <-- select the first first order derivative
  aFor.gradient().value() = 2.0;    // <-- select the second first order derivative
  aFor.gradient().gradient() = 4.0; // <-- select the second order derivative
~~~~
We have to manually select the correct data, which is still quite easy for second order types, but will get more involved for higher order types as it is shown in [Tutorial 7.2](@ref Tutorial7_2)

The derivative helper [DerivativeHelper](@ref codi::DerivativeHelper) can also be used to select specific derivatives in a more convinient fashion.
The function [derivative](@ref codi::DerivativeHelper::derivative) of the helper structure can be used to select all the derivatives.
The "order" parameter will give the order of the derivative, e.g. 1 for first order derivatives, and the "l" parameter will select the l-th derivative.
"l" can go from 0 to \f$s - 1\f$.

If we want to set all the first order derivatives for the "t2s" type, then we can do this with:
~~~~{.cpp}
    typedef codi::DerivativeHelper<t2s> DH;
    t2s aFor = 2.0;

    DH::derivative(aFor, 1, 0) = 1.0; // set the first first order derivative
    DH::derivative(aFor, 1, 1) = 1.0; // set the second first order derivative
~~~~

This can also be done in a loop or with the other helper function [setDerivatives](@ref codi::DerivativeHelper::setDerivatives), that will also set all derivatives.
If the function is used then the above code will look like:
~~~~{.cpp}
    typedef codi::DerivativeHelper<t2s> DH;
    t2s aFor = 2.0;

    DH::setDerivatives(aFor, 1, 1.0);
~~~~

We are now able to compute the second order derivative by setting all the first order derivatives with the derivative helper.
The code for the second order example is:

~~~~{.cpp}
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
~~~~

This will produce the output:
~~~~
t0s:   384
t1_1s: 1344
t1_2s: 1344
t2s:   4032
~~~~
The results are same as if they would be computed with the equations from the start of this tutorial.

We can now compute the sixth order derivative with the code:
~~~~{.cpp}
  {
    t6s aFor = 2.0;

    // set all first order directions in order to get the 6. order derivative
    typedef codi::DerivativeHelper<t6s> DH;
    DH::setDerivatives(aFor, 1, 1.0);

    t6s cFor = func(aFor);

    cout << "t0s: " << cFor << std::endl;
    cout << "t6s: " << DH::derivative(cFor, 6, 0) << std::endl;
  }
~~~~
This will produce the output:
~~~~
t0s: 384
t6s: 30240
~~~~
The results are same as if they would be computed with the equations from the start of this tutorial.

Higher order derivatves can also be computed with the reverse AD mode.
It is very complex and also not very intuitive to manage several tapes and it does not yield any improvements, therefore it is advisable, that the reverse AD type is only used as the most outer type.
That is, lhe "r6s" types uses the "t5s" forward type as the nested computation type.

The derivative helper has the same functions for the reverse types but more care has to be taken when all first order derivatives are set.
The reverse AD mode of the function \f$f\f$ is described by
\f{align*}{
  y =& f(x) \\
  \bar x =& \frac{\d f}{\d x}(x) \bar y \eqdot
\f}
If we now apply the forward AD mode to these equations, we also get second order derivatives:
\f{align*}{
  y =& f(x) \\
  \dot y =& \frac{\d f}{\d x}(x) \dot x \\
  \bar x =& \frac{\d f}{\d x}(x) \bar y \\
  \dotb x =& \frac{\d^2 f}{\d^2 x}(x) \dot x \bar y + \frac{\d f}{\d x}(x) \dotb y\eqdot
\f}
The fourth equation shows, that we have to set \f$\dot x\f$ and \f$\bar y\f$ in order to get the second order derivative.
The difference is, that \f$\dot x\f$ needs to be set before the function \f$f\f$ is evaluated and \f$\bar y\f$ needs to be set before the reverse mode is evaluated.
For this purpose the derivative helper has two additional functions [setDerivativesForward](@ref codi::DerivativeHelper::setDerivativesForward) and [setDerivativesReverse](@ref codi::DerivativeHelper::setDerivativesReverse), which will set all derivatives of the forward run and the reverse run respectively.

We can now calculate the sixth order derivatives with the reverse mode:
~~~~{.cpp}
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
~~~~
The code uses both methods in order to set the first order derivativs of the appropritate types.
The output of the code is:
~~~~
r0s: 384
r6s: 30240
~~~~
This are the same results as for the forward mode.

With the derivative helper it is now quite easy to set specific derivatives of higher order derivative types.
There are four methods that can be used:
  - [derivative](@ref codi::DerivativeHelper::derivative) for setting single derivatives
  - [setDerivatives](@ref codi::DerivativeHelper::setDerivatives) for setting all derivatives of a specific order
  - [setDerivativesForward](@ref codi::DerivativeHelper::setDerivativesForward) for setting all derivatives of a specific order for the forward run
  - [setDerivativesReverse](@ref codi::DerivativeHelper::setDerivativesReverse) for settingg all derivatives of a specific order for the reverse run

The advantage of these methods is, that the order and position can by defined at runtime.
The restriction with this method is, that all primal values and gradient values need to be of the same type.
Otherwise the compiler will create conversion errors, when the derivative helper is used.

The full code for this tutorial is:
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
~~~~
