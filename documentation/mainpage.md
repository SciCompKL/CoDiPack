CoDiPack: Fast gradient evaluation in C++ based on Expression Templates.   {#mainpage}
============

CoDiPack (Code Differentiation Package) is a tool for gradient evaluation in computer programs. It supports the features:
  - Forward mode of Algorithmic Differentiation (AD)
  - Reverse mode of Algorithmic Differentiation (AD)
  - Different tape implementations
  - An AdjointMPI interface
  - External functions
  - Higher order derivatives

The design principle for CoDiPack is that it is easy to use.
However, it also gives experienced AD developers full access to all the data structures.

The Scientific Computing Group at the TU Kaiserslautern develops CoDiPack and
will enhance and extend CoDiPack in the future.

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

For the handling of libraries and the memory optimization of the tape there exist several helper structures.
Most of them are introduced in the tutorial section:
 - codi::ExternalFunctionHelper
   - Handle external libraries which can not be handled with AD
   - Optimize large code regions
 - codi::PreaccumulationHelper
   - Reduce memory for code section that have few input and output values but
     are expensive to compute
 - codi::StatementPushHelper
   - Reduce the memory for small code fragments where the derivatives are available from an external source
 - codi::CustomAdjointVectorHelper
   - Evaluate reverse tapes with different vector settings
   - No recompilation of the whole application on a vector dimension change
 - codi::DerivativeAccess
   - More intuitive handling of higher order derivatives

Please visit the \ref TutorialsAndExamples "tutorial page" for further information.

For a full type list of the file 'codi.hpp' please see \ref ActiveTypeList.

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

You can get CoDiPack from https://www.scicomp.uni-kl.de/software/codi.

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
