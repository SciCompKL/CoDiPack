Tutorial 7.1: Higher order derivatives with the template interface {#Tutorial7_1}
============

This is the same tutorial as [Tutorial 7](@ref Tutorial7) but it will use the template version of the previous methods, therefore just the difference in the usage is explained.

For all of the four methods, that have already been introduced, a different version exists, where the parameters for the order and the index are template parameters.
These methods have the same functionality but can not be used with runtime variables, because all arguments to templates need to be known at compile time.
The advantage of this restriction is, that the compiler can optimize the code, such that no function calls or other computations are performed.
This makes the template methods the optimal choice if they are used in time critical portions of the code.

The template versions have another advantage, because they can also be used with mixed type gradient and primal values, e.g. the gradient can be a vector type.

In [Tutorial 7](@ref Tutorial7) we used the following code to set the first order derivatives:
~~~~{.cpp}
    typedef codi::DerivativeHelper<t2s> DH;
    t2s aFor2 = 2.0;

    DH::derivative(aFor, 1, 0) = 1.0; // set the first first order derivative
    DH::derivative(aFor, 1, 1) = 1.0; // set the second first order derivative
~~~~
If the template version of the function is now used, the code changes to:
~~~~{.cpp}
    typedef codi::DerivativeHelper<t2s> DH;
    t2s aFor2 = 2.0;

    DH::derivative<1,0>(aFor) = 1.0; // set the first first order derivative
    DH::derivative<1,1>(aFor) = 1.0; // set the second first order derivative
~~~~
Instead of giving the parameters to the function call, the parameters are now given as template arguments.

This same change also applies to the helper function [setDerivatives](@ref codi::DerivativeHelper::setDerivatives(Real& value, const Type& derivative)).
If we use this function, then the above code will look like:
~~~~{.cpp}
    typedef codi::DerivativeHelper<t2s> DH;
    t2s aFor2 = 2.0;

    DH::setDerivatives<1>(aFor2, 1.0);
~~~~

The full code for the second order derivatives is now:
~~~~{.cpp}
  {
    typedef codi::DerivativeHelper<t2s> DH;

    t2s aFor2 = 2.0;
    // set all first order directions in order to get the 2. order derivative
    DH::setDerivatives<1>(aFor2, 1.0);

    t2s cFor2 = func(aFor2);

    cout << "t0s:   " << DH::derivative<0, 0>(cFor2) << std::endl;
    cout << "t1_1s: " << DH::derivative<1, 0>(cFor2) << std::endl;
    cout << "t1_2s: " << DH::derivative<1, 1>(cFor2) << std::endl;
    cout << "t2s:   " << DH::derivative<2, 0>(cFor2) << std::endl;
  }
~~~~
This produces the same output as in the previous tutorial.
~~~~
t0s:   384
t1_1s: 1344
t1_2s: 1344
t2s:   4032
~~~~

The changes for the sixth order derivatives with the forward and reverse mode are the same and therefore we will not show them in detail.
They are included in the full code example of this tutorial.

As mentioned at the beginning of this tutorial, it is also possible to handle vector types with the templated functions.
Therefore, we introduce the second order type with a vector adjoint direction:
~~~~{.cpp}
typedef codi::RealForwardGen<double, codi::Direction<double, 2>> t1v;
typedef codi::RealForwardGen<t1v>                                t2v;
~~~~
The derivative function of the helper structure can now be used to specify the first order derivatives:
~~~~{.cpp}
    typedef codi::DerivativeHelper<t2v> DH;

    t2v aFor2 = 2.0;
    DH::derivative<1, 0>(aFor2) = {1.0, 2.0};
    DH::derivative<1, 1>(aFor2) = 1.0;
~~~~
That there are a double value and a vector as a derivative comes from the t1v type, where the primal type is a double.
This type is used in the t2v type as the gradient data and therefore the derivatives can have double values.
Which derivative is a vector and which is a double value can be determined by the derivative order and the selected derivative.
If the derivative order is even, then the first selected derivative will be a double value.
If the derivative order is uneven, then the first selected derivative will be a vector.
The other derivatives of the same order will be alternating between double values and vector values.

The full code for the vector example is then:
~~~~{.cpp}
  {
    typedef codi::DerivativeHelper<t2v> DH;

    t2v aFor2 = 2.0;
    // set all first order directions in order to get the 2. order derivative
    DH::derivative<1, 0>(aFor2) = {1.0, 2.0};
    DH::derivative<1, 1>(aFor2) = 1.0;

    t2v cFor2 = func(aFor2);

    cout << "t0v:   " << DH::derivative<0, 0>(cFor2) << std::endl;
    cout << "t1_1v: " << DH::derivative<1, 0>(cFor2) << std::endl;
    cout << "t1_2v: " << DH::derivative<1, 1>(cFor2) << std::endl;
    cout << "t2v:   " << DH::derivative<2, 0>(cFor2) << std::endl;
  }
~~~~

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

typedef codi::RealForwardGen<double, codi::Direction<double, 2>> t1v;
typedef codi::RealForwardGen<t1v>                                t2v;

template<typename T>
T func(const T& x) {
  T t = x * x * x * x * x * x * x;
  return t * 3.0;
}

int main() {

  {
    typedef codi::DerivativeHelper<t2s> DH;

    t2s aFor2 = 2.0;
    // set all first order directions in order to get the 2. order derivative
    DH::setDerivatives<1>(aFor2, 1.0);

    t2s cFor2 = func(aFor2);

    cout << "t0s:   " << DH::derivative<0, 0>(cFor2) << std::endl;
    cout << "t1_1s: " << DH::derivative<1, 0>(cFor2) << std::endl;
    cout << "t1_2s: " << DH::derivative<1, 1>(cFor2) << std::endl;
    cout << "t2s:   " << DH::derivative<2, 0>(cFor2) << std::endl;
  }

  {
    t6s aFor = 2.0;

    // set all first order directions in order to get the 6. order derivative
    typedef codi::DerivativeHelper<t6s> DH;
    DH::setDerivatives<1>(aFor, 1);

    t6s cFor = func(aFor);

    cout << "t0s: " << cFor << std::endl;
    cout << "t6s: " << DH::derivative<6, 0>(cFor) << std::endl;
  }

  {
    typedef codi::DerivativeHelper<r6s> DH;

    r6s::TapeType& tape = r6s::getGlobalTape();
    r6s aRev = 2.0;
    // set all first order directions on the primal value
    DH::setDerivativesForward<1>(aRev, 1.0);

    tape.setActive();
    tape.registerInput(aRev);

    r6s cRev = func(aRev);

    tape.registerOutput(cRev);
    // set all first order directions on the adjoint value
    DH::setDerivativesReverse<1>(cRev, 1.0);

    tape.setPassive();
    tape.evaluate();

    cout << "r0s: " << cRev << std::endl;
    cout << "r6s: " << DH::derivative<6, 0>(aRev) << std::endl;
  }

  {
    typedef codi::DerivativeHelper<t2v> DH;

    t2v aFor2 = 2.0;
    // set all first order directions in order to get the 2. order derivative
    DH::derivative<1, 0>(aFor2) = {1.0, 2.0};
    DH::derivative<1, 1>(aFor2) = 1.0;

    t2v cFor2 = func(aFor2);

    cout << "t0v:   " << DH::derivative<0, 0>(cFor2) << std::endl;
    cout << "t1_1v: " << DH::derivative<1, 0>(cFor2) << std::endl;
    cout << "t1_2v: " << DH::derivative<1, 1>(cFor2) << std::endl;
    cout << "t2v:   " << DH::derivative<2, 0>(cFor2) << std::endl;
  }

  return 0;
}
~~~~
