Tutorial B1.1: Function objects for the evaluation helper {#TutorialB1_1}
============

This tutorial explains in more details than tutorial [B1](@ref TutorialB1) how function objects can be implemented
such that they can be used together with the codi::EvaluationHelper.

The basic [tutorial](@ref TutorialB1) for the codi::EvaluationHelper computes the angle between two vectors \f$a\f$
and \f$b\f$. The equation is:
\f[
  \alpha = f(a, b) = \arccos\left(\frac{\scalar{a}{b}}{\norm{a} \norm{b}}\right)
\f]

The implementation of the function is:
~~~~{.cpp}
template<typename Real>
void dotWithNorms(Real const* a, Real const* b, size_t n, Real& alpha, Real& aNorm, Real& bNorm) {
  alpha = Real(); // Dot product is accumulated in alpha
  aNorm = Real();
  bNorm = Real();

  for(size_t i = 0; i < n; i += 1) {
    alpha += a[i] * b[i];
    aNorm += a[i] * a[i];
    bNorm += b[i] * b[i];
  }

  aNorm = sqrt(aNorm);
  bNorm = sqrt(bNorm);
  alpha = acos(alpha / (aNorm * bNorm));
}
~~~~

As in the basic tutorial we need to wrap this function now such that the codi::EvaluationHelper can call it from its
routines. As always there are multiple choices available, which are a function pointer, a function object and a lambda
function. But before we present the different ways of implementing the functions objects we need to have a look which
interface is expected from the codi::EvaluationHelper.

In tutorial [B1](@ref TutorialB1) the expected interface is presented as:
~~~~{.cpp}
template<typename Real>
void func(std::vector<Real> const &x, std::vector<Real> &y);
~~~~
which is already generalized for the computation type. The next codi::EvaluationHelper will call the function object
as `func(x,y)` where the type of `x` and `y` are vectors of CoDiPack types. The vector type can be specified by the
user in the more generalized `evalHandle` versions of the codi::EvaluationHelper. The expected function interface can
then be generalized for the vectors:
~~~~{.cpp}
template<template<typename> class VecX, template<typename> class VecY, typename Real>
void func(VecX<Real> const &x, VecY<Real> &y);
~~~~
`VecX` and `VecY` are the types of the vectors the user can provide for the storage of the data. `Real` is the CoDiPack
type used for the evaluation. Since the instatiation of such templates is quite cumbersome and also the autodeduction
of the template parameters might have some problems, the template template paraemters can be left out and be replaced
with regular template parameters:
~~~~{.cpp}
template<typename VecX, typename VecY>
void func(VecX const &x, VecY &y);
~~~~
With this function definition an autodeduction will have no problems. The instatiation of this method would now look like
~~~~{.cpp}
func<std::vector<codi::EvaluationHelper::HessianComputationType>, std::vector<codi::EvaluationHelper::HessianComputationType>
~~~~
for the default values. Since this definition is quite cumbersome we are now presenting different techniques for the
implementation of function objects.


Structures as function objects
------------------------------

Structures can act as function objects if there functional call operator, that is `operator()`, is overloaded. Since
this operator can have template paraemters the user no longer needs to know which kind of vectors or CoDiPack types are
used. In addition the structure can have members which provide additional information for the call operator. The wrapper
object for the function `dotWithNorms` could then be implemented like:
~~~~{.cpp}
struct WrapperDotWithNorms {
  size_t n;

  WrapperDotWithNorms(size_t n) : n(n) {}

  template<typename VecX, typename VecY>
  void operator() (VecX const &x, VecY &y) {
    dotWithNorms(&x[0], &x[this->n], this->n, y[0], y[1], y[2]);
  }
};
~~~~
The member `n` stores now the size of the vectors and does no longer need to be deduced from the input vector `x`.
The instatiation and use of the wrapper object is now quite simple:
~~~~{.cpp}
using EH = codi::EvaluationHelper;

WrapperDotWithNorms wrapDotWithNorms(n);

EH::evalJacobian(wrapDotWithNorms, x, 3, jac);
EH::evalHessian(wrapDotWithNorms, x, 3, hes);
~~~~
The advantage is now that in the `evalJacobian` and `evalHessian` calls the user does not need to specify the used
CoDiPack type. Everything is auto deduced by the compiler during the instantiation of the functions. The structure
function object makes it also simpler to handle arguments, that are not active types.

Lambda functions
----------------

Lambda function have been introduced in C++11 and have been extended in C++14. The restrictions on lambda functions in
C++11 require a similar awarenes of the used vector and CoDiPack types in the wrapper implementation as with the
function definition in tutorial [B1](@ref TutorialB1). The implementation and call would be:
~~~~{.cpp}
using EH = codi::EvaluationHelper;

auto jac = EH::createJacobian(3, xSize);
auto hes = EH::createHessian(3, xSize);

auto lambdaWrapDotWithNorms = [n](std::vector<EH::HessianComputationType> const &x, std::vector<EH::HessianComputationType> &y) {
  dotWithNorms(&x[0], &x[n], n, y[0], y[1], y[2]);
};

EH::evalJacobianAndHessian(lambdaWrapDotWithNorms, x, 3, jac, hes);
~~~~
With the lambda functions in C++11 it is not possible to add a template parameter for the CoDiPack evaluation type.
Therefore the Hessian type is used and the combined evaluation procedure. Otherwise two seperate lambda functions would
need to be specified.

For C++14 the situation is much simpler. Here introduction of generic lambdas provides us with the generalization such
that we can define the function parameters with the auto variable:
~~~~{.cpp}
using EH = codi::EvaluationHelper;

auto lambdaWrapDotWithNorms = [n](auto const&x, auto &y) {
  dotWithNorms(&x[0], &x[n], n, y[0], y[1], y[2]);
};

EH::evalJacobian(lambdaWrapDotWithNorms, x, 3, jac);
EH::evalHessian(lambdaWrapDotWithNorms, x, 3, hes);
~~~~
With this declaration it is now possible to use the lambda function with different CoDiPack types and vector types.

The full code for the tutorial is:
~~~~{.cpp}
#include <codi.hpp>

#include <iostream>

template<typename Real>
void dotWithNorms(Real const* a, Real const* b, size_t n, Real& alpha, Real& aNorm, Real& bNorm) {
  alpha = Real(); // Dot product is accumulated in alpha
  aNorm = Real();
  bNorm = Real();

  for(size_t i = 0; i < n; i += 1) {
    alpha += a[i] * b[i];
    aNorm += a[i] * a[i];
    bNorm += b[i] * b[i];
  }

  aNorm = sqrt(aNorm);
  bNorm = sqrt(bNorm);
  alpha = acos(alpha / (aNorm * bNorm));
}

struct WrapperDotWithNorms {
  size_t n;

  WrapperDotWithNorms(size_t n) : n(n) {}

  template<typename VecX, typename VecY>
  void operator() (VecX const &x, VecY &y) {
    dotWithNorms(&x[0], &x[this->n], this->n, y[0], y[1], y[2]);
  }
};

void printVector(std::string const& name, std::vector<double> const& v, size_t length, size_t offset) {
  std::cout << "Vector " << name << ": {";
  for(size_t i = 0; i < length; i += 1) {
    if(i != 0) {
      std::cout << ", ";
    }
    std::cout << v[offset + i];
  }
  std::cout << "}" << std::endl;
}

template<typename Jac>
void printJacCol(std::string const& text, Jac const &jac, size_t col) {
  std::cout << text <<": {";
  for(size_t j = 0; j < jac.getN(); j += 1) {
    if(j != 0) {
      std::cout << ", ";
    }
    std::cout << jac(col, j);
  }
  std::cout << "}" << std::endl;
}

template<typename Hes>
void printHesForOutput(std::string const& text, Hes const &hes, size_t output) {
  std::cout << text <<": {\n";
  for(size_t j = 0; j < hes.getN(); j += 1) {
    std::cout << "  ";
    for(size_t k = 0; k < hes.getN(); k += 1) {
      if(k != 0) {
        std::cout << ", ";
      }
      std::cout << hes(output, j, k);
    }
    std::cout << "\n";
  }
  std::cout << "}" << std::endl;
}

int main(int nargs, char** args) {

  int mode = 1;
  if(2 <= nargs) {
    mode = std::stoi(args[1]);

    if(mode < 1 || 3 < mode) {
      std::cerr << "Error: Please enter a mode from 1 to 3, it was '" << mode << "'." << std::endl;
      std::cerr << "  Mode 1: Function object" << std::endl;
      std::cerr << "  Mode 2: C++11 lambda" << std::endl;
      std::cerr << "  Mode 1: C++14 generic lambda" << std::endl;

      exit(-1);
    }
  }

  const size_t n = 10;
  const size_t xSize = 2 * n;

  std::vector<double> x(xSize);
  for(size_t i = 0; i < n; i += 1) {
    // vector a
    x[0 + i] = i;
    // vector b
    x[n + i] = pow(-1, i);
  }

  using EH = codi::EvaluationHelper;

  auto jac = EH::createJacobian(3, xSize);
  auto hes = EH::createHessian(3, xSize);

  if(1 == mode) { // Function object

    std::cout << "Using a structure function object." << std::endl;
    WrapperDotWithNorms wrapDotWithNorms(n);

    EH::evalJacobian(wrapDotWithNorms, x, 3, jac);
    EH::evalHessian(wrapDotWithNorms, x, 3, hes);

  } else if(2 == mode) { // C++11 lambda

    std::cout << "Using a C++11 lambda." << std::endl;
    auto lambdaWrapDotWithNorms = [n](std::vector<EH::HessianComputationType> const &x, std::vector<EH::HessianComputationType> &y) {
      dotWithNorms(&x[0], &x[n], n, y[0], y[1], y[2]);
    };

    EH::evalJacobianAndHessian(lambdaWrapDotWithNorms, x, 3, jac, hes);
  } else if(3 == mode) { // C++14 generic lambda
#if 201402L <= __cplusplus

    std::cout << "Using a C++14 generic lambda." << std::endl;
    auto lambdaWrapDotWithNorms = [n](auto const&x, auto &y) {
      dotWithNorms(&x[0], &x[n], n, y[0], y[1], y[2]);
    };

    EH::evalJacobian(lambdaWrapDotWithNorms, x, 3, jac);
    EH::evalHessian(lambdaWrapDotWithNorms, x, 3, hes);
#else
    std::cerr << "Error: Compile with C++14 to use generic lambdas." << std::endl;
    exit(-1);
#endif
  } else {
    std::cerr << "Error: Undefined mode '" << mode << "'." << std::endl;
    exit(-1);
  }

  printVector("a", x, n, 0);
  printVector("b", x, n, n);
  std::cout << std::endl;
  printJacCol("Jacobian with respect to alpha: ", jac, 0);
  printJacCol("Jacobian with respect to aNorm: ", jac, 1);
  printJacCol("Jacobian with respect to bNorm: ", jac, 2);
  std::cout << std::endl;
  printHesForOutput("Hessian with respect to alpha: ", hes, 0);
  printHesForOutput("Hessian with respect to aNorm: ", hes, 1);
  printHesForOutput("Hessian with respect to bNorm: ", hes, 2);

  return 0;
}
~~~~
