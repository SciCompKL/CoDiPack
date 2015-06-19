# CoDiPack

[CoDiPack](http://www.scicomp.uni-kl.de/software/codi/) (Code Differentiation Package) is a tool for gradient evaluation in computer programs. It supports the features:
  - Forward mode of Algorithmic Differentiation(AD)
  - Reverse mode of Algorithmic Differentiation(AD)
  - Different tape implementations
  - An AdjointMPI interface
  - External functions
  - Higher order derivatives

The design principle for CoDiPack is that it is easy to use.
However, it also gives experienced AD developers the full access to all the data structures.

The [Scientific Computing Group](http://www.scicomp.uni-kl.de) at the TU Kaiserslautern develops CoDiPack and
will enhance and extend CoDiPack in the future.

## Usage

CoDiPack is a header only library.
The only file the user needs to include is `codi.hpp`.
The only other requirement is a c++11 compliant compiler
where one usually needs to specify '--std=c++11' in compiler arguments.
CoDiPack is tested with gcc and the intel compiler.

The file `codi.hpp` defines the datatypes `RealForward`, `RealReverse` and `RealReverseUnchecked`.
For RealForward type implements the forward mode of AD and
the RealReverse type implements the reverse mode of AD.
The third type is also an implementation of the reverse mode of AD but it should only be used by experienced users.
For each type there is also a type with single precession e.g. RealForwardFloat.

## Hello World Example

A very small and simple example for the usage of the RealForward type is the code:

~~~~{.cpp}
    #include <codi.hpp>
    #include <iostream>

    int main(int nargs, char** args) {
      codi::RealForward x = 4.0;
      x.setGradient(1.0);

      codi::RealForward y = x * x;

      std::cout << "f(4.0) = " << y << std::endl;
      std::cout << "df/dx(4.0) = " << y.getGradient() << std::endl;

      return 0;
    }
~~~~

It is compiled with
~~~~{.txt}
  g++  -I<path to codi>/include --std=c++11 -g -o forward forward.cpp
~~~~
for the gcc compiler or with
~~~~{.txt}
  icpc  -I<path to codi>/include --std=c++11 -g -o forward forward.cpp
~~~~
for the intel compiler.

Please visit the [tutorial page](http://www.scicomp.uni-kl.de/codi/db/d3c/tutorialPage.html) for further information.
