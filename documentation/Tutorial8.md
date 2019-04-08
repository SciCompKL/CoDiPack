Tutorial 8: Several tape recordings {#Tutorial8}
============

In this tutorial we want to compute the same function as in [Tutorial 2](@ref Tutorial2),
but now the function will be evaluated at several points.

This will require the same steps as in Tutorial 2 and in addition the tape needs to be reset for the
different evaluations.
If this is not done, the tape will grow with each successive function evaluation.
In addition to the increased memory requirements the evaluation time for the tape will increase with each pass, too.

The recording and evaluation procedure from Tutorial 2 is:
~~~~{.cpp}
      codi::RealReverse x = 4.0;

      codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
      tape.setActive();

      tape.registerInput(x);
      codi::RealReverse y = func(x);
      tape.registerOutput(y);

      tape.setPassive();
      y.setGradient(1.0);
      tape.evaluate();

      std::cout << "f(4.0) = " << y << std::endl;
      std::cout << "df/dx(4.0) = " << x.getGradient() << std::endl;
~~~~

This procedure is now changed such that multiple evaluations are performed:
~~~~{.cpp}
      codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();

      const size_t n = 5;
      double points[n] = {2.0, 2.1, 2.5, 3.0, -1.0};

      for(int i = 0; i < n; i += 1) {
        codi::RealReverse x = points[i];

        ...

      }
~~~~

The change is quite simple but would lead to the above mentioned problem, that the tape size would increase with each loop iteration.
'setActive' does not reset the in CoDiPack since we want to be able to exclude code regions from the taping process.
The reset of the tape has to be done manually with the 'reset' function of the ReverseTapeInterface.

This adds just one call inside of the loop:
~~~~{.cpp}
      codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();

      const size_t n = 5;
      double points[n] = {2.0, 2.1, 2.5, 3.0, -1.0};

      for(int i = 0; i < n; i += 1) {
        tape.reset();

        codi::RealReverse x = points[i];

        ...

      }
~~~~

Now the tape will have always the same size in each loop iteration.

The full code of the example is now:
~~~~{.cpp}
#include <codi.hpp>
#include <iostream>

codi::RealReverse func(const codi::RealReverse& x) {
  return x * x * x;
}

void call(bool reset, bool stats) {
  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();

  const size_t n = 5;
  double points[n] = {2.0, 2.1, 2.5, 3.0, -1.0};

  for(int i = 0; i < n; i += 1) {
    if(reset) {
      tape.reset();
    }
    codi::RealReverse x = points[i];

    tape.setActive();

    tape.registerInput(x);
    codi::RealReverse y = func(x);
    tape.registerOutput(y);

    tape.setPassive();
    y.setGradient(1.0);
    tape.evaluate();

    std::cout << "f(" << x.value() << ") = " << y << std::endl;
    std::cout << "df/dx(" << x.value() << ") = " << x.getGradient() << std::endl;

    if(stas) {
      tape.printStatistics();
    }
  }
}
int main(int nargs, char** args) {

  call(false, true);

  call(true, true);

  return 0;
}
~~~~
