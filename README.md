# CoDiPack

[CoDiPack](http://www.scicomp.uni-kl.de/software/codi/) (Code Differentiation Package) is a tool for gradient evaluation in computer programs. It supports the features:
  - Forward mode of Algorithmic Differentiation (AD)
  - Reverse mode of Algorithmic Differentiation (AD)
  - Different tape implementations
  - An AdjointMPI interface
  - External functions
  - Higher order derivatives

The design principle for CoDiPack is that it is easy to use.
However, it also gives experienced AD developers full access to all the data structures.

The [Scientific Computing Group](http://www.scicomp.uni-kl.de) at the TU Kaiserslautern develops CoDiPack and will enhance and extend CoDiPack in the future.
There is a newsletter available at [codi-info@uni-kl.de](https://lists.uni-kl.de/uni-kl/subscribe/codi-info) and if you want to contact us please write a mail to [codi@scicomp.uni-kl.de](mailto:codi@scicomp.uni-kl.de).

[![Build Status](https://travis-ci.org/SciCompKL/CoDiPack.svg?branch=develop)](https://travis-ci.org/SciCompKL/CoDiPack)
[![DOI](https://zenodo.org/badge/37602249.svg)](https://zenodo.org/badge/latestdoi/37602249)

## Usage

CoDiPack is a header only library.
The only file the user needs to include is `codi.hpp`.
The only other requirement is a C++11 compliant compiler
where one usually needs to specify '-std=c++11' in compiler arguments.
CoDiPack is tested with gcc, clang, and the Intel compiler.

The file `codi.hpp` defines several datatypes. The most important ones are:
 - Implementations of the forward mode of AD:
   - codi::RealForward
 - Implementation of the reverse mode of AD:
   - codi::RealReverse (most common type, works everywhere, C-compatible)
   - codi::RealReverseIndex (reduced tape size w.r.t. codi::RealReverse, no C-like memory operations (e.g. memcpy))
   - codi::RealReversePrimal (reduced tape size w.r.t. codi::RealReverseIndex, works everywhere, C-compatible, increased type complexity)
   - codi::RealReversePrimalIndex (reduced tape size w.r.t. codi::RealReversePrimal, no C-like memory operations (e.g. memcpy), increased type complexity)

We recommend to use the codi::RealReverse type when AD is first introduced to an application.
After that there should be no difficulties in replacing the codi::RealReverse type with other types.

For further details please visit our [CoDiPack](http://www.scicomp.uni-kl.de/software/codi/) web page.

### CMake

CMake should be able to find CoDiPack either if `CMAKE_PREFIX_PATH` contains the CoDiPack directory or if the parameter `CoDiPack_DIR` is provided to CMake. If you install CoDiPack into a directory which is in the default search path for CMake then you do not need to specify any additional path.

The path is different if you use the CoDiPack directory or a CMake installation of CoDiPack.

 - CoDiPack directory (e.g. github checkout):
   - `export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:<CoDiPack root>/cmake`
   - `cmake . -DCoDiPack_DIR=<CoDiPack root>/cmake`
 - CoDiPack CMake installation:
   - `export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:<CoDiPack install>/share/CoDiPack/cmake`
   - `cmake . -DCoDiPack_DIR=<CoDiPack install>/share/CoDiPack/cmake`

You also have to add CoDiPack as a target link library in the `CMakeLists.txt` of your project:
~~~~{.cmake}
find_package(CoDiPack CONFIG REQUIRED)
target_link_libraries(<target> CoDiPack)
~~~~

## Miscellaneous information

### Debugging with gdb

The ActiveReal type contains the tape as a static member.
GDB prints the information of these members in its default settings, which makes the output quite verbose.
We recommend to disable the output of the static class members.
This can be done with
~~~~{.txt}
set print static-members off
~~~~

### Intel compiler options

Because CoDiPack relies on inlining of the compiler the performance can drop if it is not done or ignored.
Therefore we recomend to force inlining of CoDiPack with the option
~~~~{.txt}
-DCODI_UseForcedInlines 
~~~~

## Hello World Example

A very small and simple example for the usage of the RealForward type is the following code:

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
  g++  -I<path to codi>/include -std=c++11 -g -o forward forward.cpp
~~~~
for the gcc compiler or with
~~~~{.txt}
  icpc  -I<path to codi>/include -std=c++11 -g -o forward forward.cpp
~~~~
for the Intel compiler.

Please visit the [tutorial page](http://www.scicomp.uni-kl.de/codi/db/d3c/tutorialPage.html) for further information.

## Citation

If you use CoDiPack in one of your applications and write a paper it would be nice if you could cite the paper
[High-Performance Derivative Computations using CoDiPack](https://dl.acm.org/doi/abs/10.1145/3356900).
~~~~{.txt}
@article{SaAlGauTOMS2019,
  title = {High-Performance Derivative Computations using CoDiPack},
  author = {M. Sagebaum, T. Albring, N.R. Gauger},
  url = {https://dl.acm.org/doi/abs/10.1145/3356900},
  year = {2019},
  date = {2019-12-01},
  journal = {ACM Transactions on Mathematical Software (TOMS)},
  volume = {45},
  number = {4},
}
~~~~
