//! [Example 17 - Evaluation helper]
#include <iostream>

#include <codi.hpp>

#include "outputHelpers.hpp"

//! [Function]
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
//! [Function]

template<typename Real>
void codiDotWithNormsWrap(std::vector<Real> const &x, std::vector<Real> &y) {
  size_t n = x.size() / 2;
  dotWithNorms(&x[0], &x[n], n, y[0], y[1], y[2]);
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

  size_t const n = 10;
  size_t const xSize = 2 * n;

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
//! [Example 17 - Evaluation helper]
