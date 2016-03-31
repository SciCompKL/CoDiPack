Tutorial 5: ReferenceActiveReal Demonstration {#Tutorial5}
============

The [ReferenceActiveReal](@ref codi::ReferenceActiveReal) can help to optimize statements where one variable is used multiple times.
AD could detect if a variable is used multiple times in a statement but this would be a n*log(n) + n operation for each statement where n is the number of arguments.
Most of the times n would be very small but it would still be a big overhead.
With the [ReferenceActiveReal](@ref codi::ReferenceActiveReal) a different approach is taken.
It hides the variable from the tape and accumulates all the computed Jacobians for the variable.
After the statement is evaluated the accumulated Jacobian is forwarded only once to the tape.
Therefore only one Jacobian needs to be stored and not one for each occurrence of the variable in the statement.

The function that demonstrates the effect is implemented as:
~~~~{.cpp}
    void func(const codi::RealReverse& x, codi::RealReverse& y) {
      y = 3.0*x*x*x*x + 5.0*x*x*x - 3.0*x*x + 2.0*x -4.0;
    }
~~~~
and the version with the reference real is:
~~~~{.cpp}
    typedef codi::ReferenceActiveReal<codi::RealReverse> RefReal;

    void funcRef(const codi::RealReverse& x, codi::RealReverse& y) {
      RefReal xRef = x;

      y = 3.0*xRef*xRef*xRef*xRef + 5.0*xRef*xRef*xRef - 3.0*xRef*xRef + 2.0*xRef -4.0;
    }
~~~~
In order to use the reference real you have to declare a variable `xRef` and initialize it with the variable that is used multiple times in the statement.
It is not possible to do this inside the statement like y = RefReal(x) * RefReal(x).
This would yield two Jacobian entries instead of one.

The driver for these two functions is analog to the ones in the previous tutorials.
It computes the derivatives for func and funcRef and prints the results.
It also stores the tape statistics for both runs and prints them at the end.

The main function is then:
~~~~{.cpp}
    int main(int nargs, char** args) {
      codi::RealReverse x = 3.14;
      codi::RealReverse y;

      codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();

      std::cout << "Func with standard codi type." << std::endl;
      tape.setActive();

      tape.registerInput(x);
      func(x, y);
      tape.registerOutput(y);

      tape.setPassive();
      std::cout << "f(3.14) = (" << y << ")" << std::endl;

      y.setGradient(1.0);
      tape.evaluate();

      std::cout << "df/dx = (" << x.getGradient() << ")" << std::endl;

      std::ostringstream run1;
      tape.printStatistics(run1);

      tape.reset();

      std::cout << "Func with reference codi type." << std::endl;
      tape.setActive();

      tape.registerInput(x);
      funcRef(x, y);
      tape.registerOutput(y);

      tape.setPassive();
      std::cout << "f(3.14) = (" << y << ")" << std::endl;

      y.setGradient(1.0);
      tape.evaluate();

      std::cout << "df/dx = (" << x.getGradient() << ")" << std::endl;

      std::ostringstream run2;
      tape.printStatistics(run2);

      tape.reset();

      std::cout << "Statistics for the standard codi type:" << std::endl;
      std::cout << run1.str() << std::endl << std::endl;

      std::cout << "Statistics for the reference codi type:" << std::endl;
      std::cout << run2.str() << std::endl << std::endl;
    }
~~~~

The output of the program is (the statistics output is trimmed to the essential parts):

~~~~
Func with standard codi type.
f(3.14) = (419.132)
df/dx = (502.564)
Func with reference codi type.
f(3.14) = (419.132)
df/dx = (502.564)

Statistics for the standard codi type:
-------------------------------------
Statements
-------------------------------------
  Total Number:              3
-------------------------------------
Jacobi entries
-------------------------------------
  Total Number:             11


Statistics for the reference codi type:
-------------------------------------
Statements
-------------------------------------
  Total Number:              3
-------------------------------------
Jacobi entries
-------------------------------------
  Total Number:              2
~~~~

The number of statements is the same but the number of Jacobian entries changed, because the ReferenceActiveType accumulates the Jacobians in the statements.

The full program is then

~~~~{.cpp}
#include <codi.hpp>
#include <iostream>
#include <sstream>

void func(const codi::RealReverse& x, codi::RealReverse& y) {
  y = 3.0*x*x*x*x + 5.0*x*x*x - 3.0*x*x + 2.0*x -4.0;
}

typedef codi::ReferenceActiveReal<codi::RealReverse> RefReal;

void funcRef(const codi::RealReverse& x, codi::RealReverse& y) {
  RefReal xRef = x;

  y = 3.0*xRef*xRef*xRef*xRef + 5.0*xRef*xRef*xRef - 3.0*xRef*xRef + 2.0*xRef -4.0;
}

int main(int nargs, char** args) {
  codi::RealReverse x = 3.14;
  codi::RealReverse y;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();

  std::cout << "Func with standard codi type." << std::endl;
  tape.setActive();

  tape.registerInput(x);
  func(x, y);
  tape.registerOutput(y);

  tape.setPassive();
  std::cout << "f(3.14) = (" << y << ")" << std::endl;

  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "df/dx = (" << x.getGradient() << ")" << std::endl;

  std::ostringstream run1;
  tape.printStatistics(run1);

  tape.reset();

  std::cout << "Func with reference codi type." << std::endl;
  tape.setActive();

  tape.registerInput(x);
  funcRef(x, y);
  tape.registerOutput(y);

  tape.setPassive();
  std::cout << "f(3.14) = (" << y << ")" << std::endl;

  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "df/dx = (" << x.getGradient() << ")" << std::endl;

  std::ostringstream run2;
  tape.printStatistics(run2);

  tape.reset();

  std::cout << std::endl;
  std::cout << "Statistics for the standard codi type:" << std::endl;
  std::cout << run1.str() << std::endl << std::endl;

  std::cout << "Statistics for the reference codi type:" << std::endl;
  std::cout << run2.str() << std::endl << std::endl;
}
~~~~
