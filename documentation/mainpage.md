CoDiPack: Fast gradient evaluation in C++ based on Expression Templates.   {#mainpage}
============

CoDiPack (Code Differentiation Package) is a tool for gradient evaluation in computer programs. It supports the features:
  - Forward mode of Algorithmic Differentiation(AD)
  - Reverse mode of Algorithmic Differentiation(AD)
  - Different tape implementations
  - An AdjointMPI interface
  - External functions
  - Higher order derivatives

The design principle for CoDiPack is that it is easy to use.
However, it also gives experienced AD developers the full access to all the data structures.

The Scientific Computing Group at the TU Kaiserslautern develops CoDiPack and
will enhance and extend CoDiPack in the future.

## Usage

CoDiPack is a header only library.
The only file the user needs to include is `codi.hpp`.
The only other requirement is a c++11 compliant compiler
where one usually needs to specify '--std=c++11' in compiler arguments.
CoDiPack is tested with gcc and the intel compiler.

The file `codi.hpp` defines several datatypes. The most important ones are:
 - Implementations of the forward mode of AD:
   - codi::RealForward
 - Implementation of the reverse mode of AD:
   - codi::RealReverse (most common type, works everywhere, c-compatible)
   - codi::RealReverseIndex (reduced tape size w.r.t. codi::RealReverse, no c-like memory operations (e.g. memcpy))
   - codi::RealReversePrimal (reduced tape size w.r.t. codi::RealReverseIndex, works everywhere, c-compatible, increased type complexity)
   - codi::RealReversePrimalIndex (reduced tape size w.r.t. codi::RealReversePrimal, no c-like memory operations (e.g. memcpy), increased type complexity)

We recommend to use the codi::RealReverse type when AD is first introduced to an application.
After that there should be no difficulties in replacing the codi::RealReverse type with other types.

The full type list of the file 'codi.hpp' is:
 - Implementations of the forward mode of AD:
   - codi::RealForward
   - codi::RealForwardFloat
 - Implementation of the reverse mode of AD:
   - codi::RealReverse
   - codi::RealReverseIndex
   - codi::RealReverseUnchecked
   - codi::RealReverseIndexUnchecked
   - codi::RealReverseFloat
   - codi::RealReverseIndexFloat
   - codi::RealReverseUncheckedFloat
   - codi::RealReverseIndexUncheckedFloat
   - codi::RealReversePrimal
   - codi::RealReversePrimalIndex
   - codi::RealReversePrimalUnchecked
   - codi::RealReversePrimalIndexUnchecked
   - codi::RealReversePrimalFloat
   - codi::RealReversePrimalIndexFloat
   - codi::RealReversePrimalUncheckedFloat
   - codi::RealReversePrimalIndexUncheckedFloat
 - Vector versions of the above AD types:
   - codi::RealForwardVec<dim>
   - codi::RealReverseVec<dim>
   - codi::RealReverseIndexVec<dim>
   - codi::RealReversePrimalVec<dim>
   - codi::RealReversePrimalIndexVec<dim>

The reverse types support various use cases. The regular type codi::RealReverse is the most used type and provides
the most common use case. This type can be used in c-like memory operation like memset and memcpy.
The 'Index' variant of the reverse type uses an indexing scheme that reuses freed indices and therefore
reduces the amount of memory that is needed. This type is no longer compatible with c-like memory operations.
The 'Primal' variants implement a different strategy for storing the data.
Instead of storing the partial derivatives for each statement, they store the primal values.
This change reduces the required memory of the 'Primal' types.
The 'Unchecked' variant is also an implementation of the reverse mode of AD but it should only be used by experienced users. This type performs no bounds checking for the memory access.
For each type there is also a type with single precession e.g. codi::RealForwardFloat.
The 'Vec' variant implements the vector mode of the corresponding AD type.
The dimension is fixed and can be defined via the template argument.

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
  g++  -I<path to codi>/include -std=c++11 -g -o forward forward.cpp
~~~~
for the gcc compiler or with
~~~~{.txt}
  icpc  -I<path to codi>/include -std=c++11 -g -o forward forward.cpp
~~~~
for the intel compiler.

Please visit the \ref tutorialPage "tutorial page" for further information.

## Citation

If you use CoDiPack in one of your applications and write a paper it would be nice if you could cite the paper
[High-Performance Derivative Computations using CoDiPack](https://arxiv.org/abs/1709.07229) (submitted to ACM TOMS).
~~~~{.txt}
@article{sagebaum2017high,
  title={{High-Performance Derivative Computations using CoDiPack}},
  author={Sagebaum, Max and Albring, Tim and Gauger, Nicolas R.},
  journal={arXiv preprint arXiv:1709.07229},
  year={2017}
}
~~~~
