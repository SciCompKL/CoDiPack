Tutorial B1.2: Handle creation and advanced uses of the evaluation helper {#TutorialB1_2}
============

This tutorial explains in more detail than tutorial [B1](@ref TutorialB1) how other CoDiPack types than the default ones
in the codi::EvaluationHelper can be used and how the performance can be improved if a function object needs to be
evaluated multiple times.

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

Tutorial [B1.1](@ref TutorialB1_1) describes how function objects can be written to wrap this function. In this tutorial
we want to focus on the handle creation and therefore we take the most general function handle which uses a structure to
create the function handle:
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

The default logic in the codi::EvaluationHelper will create a default handle internally every time an `eval...` function
is called. With the handle creation there is also the creation of the internal data structures which are required to
convert the user provided evaluation points into the CoDiPack data structures. For single calls to an `eval...` function
this is fine but if the function is called millions of times, the overhead will be noticeable. Therefore it is possible
to create the handle up front and then reuse the same data structures for all evaluation points. A second advantage of
the handle creation is that the user can choose the CoDiPack type for the handle. Since the default handle uses the
forward mode it can be very slow for Jacobian or Hessian computations where there are more input values than output
values.

The general functions for the default handle creation are
[createHandleDefault](@ref codi::EvaluationHelper::createHandleDefault),
[createHandleDefaultFixed](@ref codi::EvaluationHelper::createHandleDefaultFixed),
[createHandleDefault2nd](@ref codi::EvaluationHelper::createHandleDefault2nd) and
[createHandleDefaultFixed2nd](@ref codi::EvaluationHelper::createHandleDefaultFixed2nd).
The [createHandleDefault](@ref codi::EvaluationHelper::createHandleDefault) will create a handle with `std::vector` as
storage containers and a first order foward CoDiPack vector type with the vector size 4. The same handle is created with
[createHandleDefaultFixed](@ref codi::EvaluationHelper::createHandleDefaultFixed), but instead `std::array` containers are
used for the storage. This function requires the size of the inputs and outputs as template arguments. The other two
functions [createHandleDefault2nd](@ref codi::EvaluationHelper::createHandleDefault2nd) and
[createHandleDefaultFixed2nd](@ref codi::EvaluationHelper::createHandleDefaultFixed2nd) do the same as there first order
counter parts but they use a second order forward CoDiPack vector type, where both vector dimensions have the size 4.
An example call of all four methods is:
~~~~{.cpp}
using EH = codi::EvaluationHelper;

WrapperDotWithNorms wrapDotWithNorms(n);
const size_t xSize = 2 * n;

auto hDef = EH::createHandleDefault(wrapDotWithNorms, 3, xSize);
auto hDefFixed = EH::createHandleDefaultFixed<3, xSize>(wrapDotWithNorms);
auto hDef2nd = EH::createHandleDefault2nd(wrapDotWithNorms, 3, xSize);
auto hDefFixed2nd = EH::createHandleDefaultFixed2nd<3, xSize>(wrapDotWithNorms);
~~~~
For the `Fixed` versions of the methods, it is necessary that the arguments for the sizes are known at compile time.
Otherwise they behave in the same way.

As we have now the handles available, they can be used in the `evalHandle...` function of the codi::EvaluationHelper.
These are the same functions as in the tutorial [B1](@ref TutorialB1) and have nearly the same argument structure. Instead
of the function object the handle needs to be provided. For all functions that did not perform a primal evaluation (e.g.
`evalHandleJacobian`) it is no longer necessary to provide the number of output arguments. An example use of the handles
is:
~~~~{.cpp}
EH::evalHandleJacobian(hDef, x, jac);
EH::evalHandleJacobian(hDefFixed, x, jac);
EH::evalHandleJacobianAndHessian(hDef2nd, x, jac, hes);
EH::evalHandleJacobianAndHessian(hDefFixed2nd, x, jac, hes);
~~~~

It is now more efficient to evaluate the function object several times:
~~~~{.cpp}
for(size_t j = 0; j < 1000000; j += 1) {
  x[0] = j * 10;
  EH::evalHandleJacobian(hDef, x, jac);

  // perform some output
}
~~~~

The other `eval...` functions from tutorial [B1](@ref TutorialB1) can be called in the same way. Now we will present the
other handle creation functions that provide more flexibility in terms of the used CoDiPack type and for the underlying data
structure.

The next two functions are [createHandle](@ref codi::EvaluationHelper::createHandle) and
[createHandleFixed](@ref codi::EvaluationHelper::createHandleFixed). The first one creates a handle with `std::vector`
as the underlying container type. The second one uses `std::array`. In both functions the used CoDiPack type is given as
a template parameter. All `Real...` type definitions in [codi.hpp](@ref codi.hpp) can be used here. Since the defintion
of higher order types can be quite involved there are two additional definitions, namely
[JacobianComputationType](@ref codi::JacobianComputationType) and [HessianComputationType](@ref codi::HessianComputationType)
in this file in the global `codi` namespace. They represent a Jacobian index management taping approach  with a vector
mode of 4 for the first one. The second one uses a primal value index mangement tape with a vector mode of 4 and the
inner AD type is a forward type with a vector dimension of 4. Some example custom handle creations are now:
~~~~{.cpp}
using EH = codi::EvaluationHelper;

WrapperDotWithNorms wrapDotWithNorms(n);
const size_t xSize = 2 * n;

auto hJac = EH::createHandle<codi::JacobianComputationType>(wrapDotWithNorms, 3, xSize);
auto hPrimalFixed2nd = EH::createHandleFixed<codi::HessianComputationType, 3, xSize>(wrapDotWithNorms);
auto hPrimal = EH::createHandle<codi::RealReversePrimalVec<4>>(wrapDotWithNorms, 3, xSize);
auto hJac2nd = EH::createHandle<codi::RealReverseGen<codi::RealForwardVec<4>, codi::Direction<codi::RealForwardVec<4>, 4>>(wrapDotWithNorms, 3, xSize);
~~~~

The final method for the most flexibility is [createHandleFull](@ref codi::EvaluationHelper::createHandleFull). In
addition to all other methods, the storage types for the input and output vectors can be given as template template
parameters. Since it is no longer possible to deduce a fixed size from the template argumetns, this method requires that
the number of input and output variables are given as arguments. Examples for the use of this method are:
~~~~{.cpp}
using EH = codi::EvaluationHelper;

WrapperDotWithNorms wrapDotWithNorms(n);
const size_t xSize = 2 * n;

auto hFull = EH::createHandleFull<codi::JacobianComputationType, codi::adapters::StdVector, codi::adapters::StdArray<3>::template Type>(wrapDotWithNorms, 3, xSize);
auto hFull2nd = EH::createHandleFull<codi::HessianComputationType, codi::adapters::StdVector, codi::adapters::StdArray<3>::template Type>(wrapDotWithNorms, 3, xSize);
~~~~
`hFull` creates a first order type where the internal storage for the input variables is a `std::vector` and the storage
for the output vectors is a `std::array` of size 3. `hFull2nd` does the same with a second order type.

All versions of `createHandle` and which defaults they use are now summarized:
<table>
<caption id="createHandle_summary">Possible versions of createHandle in the CoDiPack evaluation helper.</caption>
<tr><th> Function <th> AD mode <th> Order <th> Storage container
<tr><td> [createHandleDefault](@ref codi::EvaluationHelper::createHandleDefault) <td> Forward Vec 4 <td> 1st <td> std::vector
<tr><td> [createHandleDefaultFixed](@ref codi::EvaluationHelper::createHandleDefaultFixed) <td> Forward Vec 4 <td> 1st <td> std::array
<tr><td> [createHandleDefault2nd](@ref codi::EvaluationHelper::createHandleDefault2nd) <td> Forward Vec 4 over Forward Vec 4 <td> 2nd <td> std::vector
<tr><td> [createHandleDefaultFixed2nd](@ref codi::EvaluationHelper::createHandleDefaultFixed2nd) <td> Forward Vec 4 over Forward Vec 4 <td> 2nd <td> std::array
<tr><td> [createHandle](@ref codi::EvaluationHelper::createHandle) <td> User defined <td> User defined <td> std::vector
<tr><td> [createHandleFixed](@ref codi::EvaluationHelper::createHandleFixed) <td> User defined <td> User defined <td> std::array
<tr><td> [createHandleFull](@ref codi::EvaluationHelper::createHandleFull) <td> User defined <td> User defined <td> User defined
</table>

The full code of the tutorial is:
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

    if(mode < 1 || 10 < mode) {
      std::cerr << "Error: Please enter a mode from 1 to 10, it was '" << mode << "'." << std::endl;
      std::cerr << "  Mode  1: createHandleDefault" << std::endl;
      std::cerr << "  Mode  2: createHandleDefaultFixed" << std::endl;
      std::cerr << "  Mode  3: createHandleDefault2nd" << std::endl;
      std::cerr << "  Mode  4: createHandleDefaultFixed2nd" << std::endl;
      std::cerr << "  Mode  5: createHandle 1st order Jacobian tape" << std::endl;
      std::cerr << "  Mode  6: createHandleFixed 2nd order primal value tape" << std::endl;
      std::cerr << "  Mode  7: createHandle 1st order primal value tape" << std::endl;
      std::cerr << "  Mode  8: createHandle 2nd order primal value tape" << std::endl;
      std::cerr << "  Mode  9: createHandleFull 1st order Jacobian tape" << std::endl;
      std::cerr << "  Mode 10: createHandleFull 2nd order primal value tape" << std::endl;


      exit(-1);
    }
  }

  constexpr size_t n = 10;
  constexpr size_t xSize = 2 * n;

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
  bool hesEval = false;

  WrapperDotWithNorms wrapDotWithNorms(n);

  switch(mode) {
  case 1: { // createHandleDefault
      auto hDef = EH::createHandleDefault(wrapDotWithNorms, 3, xSize);
      EH::evalHandleJacobian(hDef, x, jac);
      break;
    }
  case 2: { // createHandleDefaultFixed
      auto hDefFixed = EH::createHandleDefaultFixed<3, xSize>(wrapDotWithNorms);
      EH::evalHandleJacobian(hDefFixed, x, jac);
      break;
    }
  case 3: { // createHandleDefault2nd
      auto hDef2nd = EH::createHandleDefault2nd(wrapDotWithNorms, 3, xSize);
      EH::evalHandleJacobianAndHessian(hDef2nd, x, jac, hes);
      hesEval = true;
      break;
    }
  case 4: { // createHandleDefaultFixed2nd
      auto hDefFixed2nd = EH::createHandleDefaultFixed2nd<3, xSize>(wrapDotWithNorms);
      EH::evalHandleJacobianAndHessian(hDefFixed2nd, x, jac, hes);
      hesEval = true;
      break;
    }
  case 5: { // createHandle 1st order Jacobian tape
      auto hJac = EH::createHandle<codi::JacobianComputationType>(wrapDotWithNorms, 3, xSize);
      EH::evalHandleJacobian(hJac, x, jac);
      break;
    }
  case 6: { // createHandleFixed 2nd order primal value tape
      auto hPrimalFixed2nd = EH::createHandleFixed<codi::HessianComputationType, 3, xSize>(wrapDotWithNorms);
      EH::evalHandleJacobianAndHessian(hPrimalFixed2nd, x, jac, hes);
      hesEval = true;
      break;
    }
  case 7: { // createHandle 1st order primal value tape
      auto hPrimal = EH::createHandle<codi::RealReversePrimalVec<4>>(wrapDotWithNorms, 3, xSize);
      EH::evalHandleJacobian(hPrimal, x, jac);
      break;
    }
  case 8: { // createHandle 2nd order Jacobian tape
      auto hJac2nd = EH::createHandle<codi::RealReverseGen<codi::RealForwardVec<4>, codi::Direction<codi::RealForwardVec<4>, 4>>>(wrapDotWithNorms, 3, xSize);
      EH::evalHandleJacobianAndHessian(hJac2nd, x, jac, hes);
      hesEval = true;
      break;
    }
  case 9: { // createHandleFull 1st order Jacobian tape
      auto hFull = EH::createHandleFull<codi::JacobianComputationType, codi::adapters::StdVector, codi::adapters::StdArray<3>::template Type>(wrapDotWithNorms, 3, xSize);
      EH::evalHandleJacobian(hFull, x, jac);
      break;
    }
  case 10: { // createHandleFull 2st order primal value tape
      auto hFull2nd = EH::createHandleFull<codi::HessianComputationType, codi::adapters::StdVector, codi::adapters::StdArray<3>::template Type>(wrapDotWithNorms, 3, xSize);
      EH::evalHandleJacobianAndHessian(hFull2nd, x, jac, hes);
      hesEval = true;
      break;
    }
  default:
    std::cerr << "Error: Undefined mode '" << mode << "'." << std::endl;
    exit(-1);
  }

  printVector("a", x, n, 0);
  printVector("b", x, n, n);
  std::cout << std::endl;
  printJacCol("Jacobian with respect to alpha: ", jac, 0);
  printJacCol("Jacobian with respect to aNorm: ", jac, 1);
  printJacCol("Jacobian with respect to bNorm: ", jac, 2);
  if(hesEval) {
    std::cout << std::endl;
    printHesForOutput("Hessian with respect to alpha: ", hes, 0);
    printHesForOutput("Hessian with respect to aNorm: ", hes, 1);
    printHesForOutput("Hessian with respect to bNorm: ", hes, 2);
  }

  return 0;
}
~~~~









