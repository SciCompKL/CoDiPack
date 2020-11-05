Tutorial B1: Evaluation helper {#TutorialB1}
============

This tutorial describes the basic use of the codi::EvaluationHelper which can be used to compute Jacobians and Hessians
of an arbitrary function in an automatic fashion. The advantage of the codi::EvaluationHelper is the nearly no knowledge 
about Algorithmic Differentiation or CoDiPack is required. The helper hides all implementational details and provides an
interface where only the function needs to be presented to the codi::EvaluationHelper. With this function it can then
compute the Hessian and Jacobian of that function.

We want to differentiate a function that computes the angle between two vectors \f$a\f$ and \f$b\f$.
Mathematically this is done by computing the dot product of their normalized vectors and taking the arcus cosine function of
the result:
\f[
  \alpha = f(a, b) = \arccos\left(\frac{\scalar{a}{b}}{\norm{a} \norm{b}}\right)
\f]

This function is implemented in such a way that the angle and both vector norms are returned:
~~~~{.cpp}
const size_t n = 10;

void dotWithNorms(double const* a, double const* b, size_t n, double& alpha, double& aNorm, double& bNorm) {
  alpha = double(); // Dot product is accumulated in alpha
  aNorm = double();
  bNorm = double();

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

Since this function has \f$ 2 \cdot n \f$ input values and 3 output values, we would like to compute the Hessian and the
Jacobian with the codi::EvaluationHelper. This class provides several helper functions for the evaluation of Jacobians
and Hessians of function objects. Which kind of function objects are accepted by the codi::EvaluationHelper is
discussed in the sub tutorial [B1.1](@ref TutorialB1_1). Here, we assume the simplest default case where the function declaration is:
~~~~{.cpp}
void func(std::vector<CoDiType> const &x, std::vector<CoDiType> &y);
~~~~
The function uses standard vectors as the array containers and consists of an input vector `x` and an output vector
`y`. The `CoDiType` can be an arbitrary CoDiPack type, but in this first tutorial we will use the default types
provided by the codi::EvaluationHelper, these are codi::EvaluationHelper::JacobianComputationType for the Jacobian
evaluation and codi::EvaluationHelper::HessianComputationType for the Hessian computation. Since we need to be able to
use these two types in the evaluation of `dotWithNorms` a template parameter is added to this function:
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
Now we can instantiate a version of `dotWithNorms` for each of the two evaluation types. In order to provide this
function to the codi::EvaluationHelper we need to bring it into the form `func(x,y)`. We could rewrite `dotWithNorms` a
second time such that this interface is met but in this tutorial we decide to write a wrapper function that will call
`dotWithNorms` and provides the required interface for the codi::EvaluationHelper. The wrapper function is:
~~~~{.cpp}
template<typename Real>
void codiDotWithNormsWrap(std::vector<Real> const &x, std::vector<Real> &y) {
  size_t n = x.size() / 2;
  dotWithNorms(&x[0], &x[n], n, y[0], y[1], y[2]);
}
~~~~
The wrapper assumes that the two vectors `a` and `b` have been packed into the `x` vector. The same assumption is made
for the output values. Here, the first output value is the angle between the vectors. The other two output variables are the
norms of the vectors.

We now have to fix a point at which we want to evaluate the Jacobian or Hessian. For this we create a vector and fill
it with the data for the vectors `a` and `b`:
~~~~{.cpp}
std::vector<double> x(xSize);
for(size_t i = 0; i < n; i += 1) {
  // vector a
  x[0 + i] = i;
  // vector b
  x[n + i] = pow(-1, i);
}
~~~~

Now everything is setup such that the codi::EvaluationHelper can be called. The code is shown now and discussed
afterwards:
~~~~{.cpp}
using EH = codi::EvaluationHelper;

auto jac = EH::createJacobian(3, xSize);
auto hes = EH::createHessian(3, xSize);

EH::evalJacobian(codiDotWithNormsWrap<EH::JacobianComputationType>, x, 3, jac);
EH::evalHessian(codiDotWithNormsWrap<EH::HessianComputationType>, x, 3, hes);
~~~~
The first line creates a shortcut such that `codi::EvaluationHelper` can be accessed via the `EH` abbreviation.
Then the storage for the Jacobian and Hessian is created. The computation of these two objects is then triggered
with calls to `evalJacobian` and `evalHessian`. Since both calls are quite similar, we will only explain the one for
the Jacobian.

The call expects the function object for which the Jacobian shall be computed, the point at which the evaluation shall
be done, the number of output values and the actual storage for the Jacobian. In our case the function object is just
the pointer to the `codiDotWithNormsWrap` wrapper function instantiated with the template parameter for the
[default Jacobian computation type](@ref codi::EvaluationHelper::JacobianComputationType). Because of the template
parameter it is now quite easy to use the function for both evaluations. The vector `x` contains the two packed vectors
`a` and `b`. This packing of the input and output variables is necessary since we needed to fix an interface for the
implementation of all the helpers in CoDiPack. The implementation in CoDiPack provides still a lot of flexibility that
we will show in the sub tutorials. The number of output variables needs to be provided for this function since
the result `y` is not requested with the call. Otherwise, the function call would not know the size for the internal structures.
The final argument needs to be the storage space for the Jacobian. A default value can be generated through the
codi::EvaluationHelper function [createJacobian](@ref codi::EvaluationHelper::createJacobian).

As already stated, the call for the Hessian is nearly the same and is not explained. The EvaluationHelper contains
several other functions which also provide functionality to e.g. compute the Hessian and Jacobian in one call.
The set of all functions is:
[evalPrimal](@ref codi::EvaluationHelper::evalPrimal),
[evalPrimalAndJacobian](@ref codi::EvaluationHelper::evalPrimalAndJacobian),
[evalPrimalAndHessian](@ref codi::EvaluationHelper::evalPrimalAndHessian),
[evalPrimalAndJacobianAndHessian](@ref codi::EvaluationHelper::evalPrimalAndJacobianAndHessian),
[evalJacobian](@ref codi::EvaluationHelper::evalJacobian),
[evalJacobianAndHessian](@ref codi::EvaluationHelper::evalJacobianAndHessian),
[evalHessian](@ref codi::EvaluationHelper::evalHessian).
With these methods the desired combination of information can be acquired without calling multiple functions. We can now
change our example such that the Jacobian and Hessian are computed with one function.:
~~~~{.cpp}
using EH = codi::EvaluationHelper;

auto jac = EH::createJacobian(3, xSize);
auto hes = EH::createHessian(3, xSize);

EH::evalJacobianAndHessian(codiDotWithNormsWrap<EH::HessianComputationType>, x, 3, jac, hes);
~~~~

This concludes the basic introduction to the codi::EvaluationHelper class. The tutorial [B1.1](@ref TutorialB1_1) will
go into more details how function objects can be defined in a more general way and which kind of function objects are
admissible. Tutorial [B1.2](@ref TutorialB1_2) will introduce the creation of handles for the codi::EvaluationHelper which
provide a speedup if a function object needs to evaluated at several different points. The handles also enable the
use of all CoDiPack types.

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

template<typename Real>
void codiDotWithNormsWrap(std::vector<Real> const &x, std::vector<Real> &y) {
  size_t n = x.size() / 2;
  dotWithNorms(&x[0], &x[n], n, y[0], y[1], y[2]);
}

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

    if(mode < 1 || 2 < mode) {
      std::cerr << "Error: Please enter a mode from 1 to 2, it was '" << mode << "'." << std::endl;
      std::cerr << "  Mode 1: evalJacobian and evalHessian call" << std::endl;
      std::cerr << "  Mode 2: evalJacobianAndHessian call" << std::endl;

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

  if(1 == mode) {
    std::cout << "Using evalJacobian and evalHessian." << std::endl;
    EH::evalJacobian(codiDotWithNormsWrap<EH::JacobianComputationType>, x, 3, jac);
    EH::evalHessian(codiDotWithNormsWrap<EH::HessianComputationType>, x, 3, hes);
  } else if(2 == mode) {
    std::cout << "Using evalJacobianAndHessian." << std::endl;
    EH::evalJacobianAndHessian(codiDotWithNormsWrap<EH::HessianComputationType>, x, 3, jac, hes);
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

