Tutorial A5: Custom Jacobian and Hessian storage {#TutorialA5}
============

In this tutorial we want to show how custom data structures can be used in the codi::EvaluationHelper for the Jacobian
and Hessian. It will create a small wrapper structure that can be used to map the data of the user to the CoDiPack
interface.

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
For the implementation of the function wrapper please see tutorial [B1.1](@ref TutorialB1_1).

In the B1 tutorials we used the function [createJacobian](@ref codi::EvaluationHelper::createJacobian) and
[createHessian](@ref codi::EvaluationHelper::createHessian) to create the data storage for the derivatives. These functions
return instances of the default CoDiPack implementations codi::Jacobian and codi::Hessian. The interfaces of these two
implementations are already quite small but the algorithms behind the codi::EvaluationHelper require only the functional
call operator, that is `operator ()`, as they are define in [JacobianInterface](@ref codi::JacobianInterface) and
[HessianInterface](@ref codi::HessianInterface). The reference access of both interfaces looks like:
~~~~{.cpp}
// struct JacobianInterface
  T& operator()(const size_t i, const size_t j);

// struct HessianInterface
  T& operator()(const size_t i, const size_t j, const size_t k);
~~~~
For both interfaces the parameter `i` iterates over the output variables and j iterates over the input variables. In the
Hessian version of the operator, the parameter `k` iterates over the input variables from the second order derivation.

In order to demonstrate how a simple wrapper for the data structure can be written. We will create a wrapper for a
Hessian. Our underlying data will be an array which is organized such that all values for one output are stored in a
column major matrix. That is `j` is the fastest running index, `k` the second one and `i` is the slowest running index.
The implementation of the wrapper is then:
~~~~{.cpp}
struct HessianPointer {

  double* data;

  size_t m;
  size_t n;

  HessianPointer(double* data, size_t m, size_t n) : data(data), m(m), n(n) {}

  double operator()(const size_t i, const size_t j, const size_t k) const {
    return data[computeIndex(i,j,k)];
  }

  double& operator()(const size_t i, const size_t j, const size_t k) {
    return data[computeIndex(i,j,k)];
  }

  size_t computeIndex(const size_t i, const size_t j, const size_t k) const {
    return i * n * m + k * n + j;
  }
};
~~~~
We store the data pointer and the sizes as member variables which are initialized in the constructor. The function
`computeIndex` calculates from `i`, `j` and `k` the offset into the data pointer. The two `operator()` functions provide
the implementation of the interface.

This warpper can now be used in the algorithms of the codi::EvaluationHelper. After creating the data manually, we create
the wrapper and provide the wrapper to the eval function call:
~~~~{.cpp}
double* hesData = new double[3 * xSize * xSize];
HessianPointer hes(hesData, 3, xSize);

using EH = codi::EvaluationHelper;

WrapperDotWithNorms wrapDotWithNorms(n);
EH::evalHessian(wrapDotWithNorms, hDef2nd, x, 3, hes);
~~~~

For the codi::JacobianInterface the implementation is nearly the same. With the wrapper a user can directly store the
computed Jacobians and Hessians in its own data structures.


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

void printHesForOutput(std::string const& text, double const *hes, size_t output, size_t m, size_t n) {
  std::cout << text <<": {\n";
  for(size_t j = 0; j < n; j += 1) {
    std::cout << "  ";
    for(size_t k = 0; k < n; k += 1) {
      if(k != 0) {
        std::cout << ", ";
      }
      std::cout << hes[k * n * m + output * n + j];
    }
    std::cout << "\n";
  }
  std::cout << "}" << std::endl;
}

struct HessianPointer {

  double* data;

  size_t m;
  size_t n;

  HessianPointer(double* data, size_t m, size_t n) : data(data), m(m), n(n) {}

  double operator()(const size_t i, const size_t j, const size_t k) const {
    return data[computeIndex(i,j,k)];
  }

  double& operator()(const size_t i, const size_t j, const size_t k) {
    return data[computeIndex(i,j,k)];
  }

  size_t computeIndex(const size_t i, const size_t j, const size_t k) const {
    return k * n * m + i * n + j;
  }
};

int main(int nargs, char** args) {

  constexpr size_t n = 10;
  constexpr size_t xSize = 2 * n;

  std::vector<double> x(xSize);
  for(size_t i = 0; i < n; i += 1) {
    // vector a
    x[0 + i] = i;
    // vector b
    x[n + i] = pow(-1, i);
  }

  double* hesData = new double[3 * xSize * xSize];
  HessianPointer hes(hesData, 3, xSize);

  using EH = codi::EvaluationHelper;

  WrapperDotWithNorms wrapDotWithNorms(n);
  EH::evalHessian(wrapDotWithNorms, hDef2nd, x, 3, hes);

  printVector("a", x, n, 0);
  printVector("b", x, n, n);
  std::cout << std::endl;
  printHesForOutput("Hessian with respect to alpha: ", hesData, 0, 3, xSize);
  printHesForOutput("Hessian with respect to aNorm: ", hesData, 1, 3, xSize);
  printHesForOutput("Hessian with respect to bNorm: ", hesData, 2, 3, xSize);

  delete [] hesData;

  return 0;
}
~~~~
