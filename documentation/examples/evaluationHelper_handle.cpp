//! [EvaluationHelper handle]
#include <codi.hpp>

#include "outputHelpers.hpp"

using Real = codi::EvaluationHelper::HessianComputationType;

void func(std::array<Real, 2> const& x, std::array<Real, 4>& y) {
  y[0] = x[0] + x[1];
  y[1] = x[0] - x[1];
  y[2] = x[0] * x[1];
  y[3] = x[0] / x[1];
}

int main(int nargs, char** args) {
  std::vector<double> x = {3.0, 4.0};
  std::vector<double> y(4);

  codi::EvaluationHelper eh;
  auto jac = eh.createJacobian(4,2);
  auto hes = eh.createHessian(4,2);

//! [EvaluationHelper changes]
  auto handle = eh.createHandleDefaultFixed2nd<4, 2>(func);
  eh.evalHandlePrimalAndJacobianAndHessian(handle, x, y, jac, hes);
//! [EvaluationHelper changes]

  std::cout << "Jacobian:" << std::endl;
  std::cout << jac << std::endl;
  printHesForOutput("Hessian with respect to y[0]: ", hes, 0);
  printHesForOutput("Hessian with respect to y[1]: ", hes, 1);
  printHesForOutput("Hessian with respect to y[2]: ", hes, 2);
  printHesForOutput("Hessian with respect to y[3]: ", hes, 3);

  return 0;
}
//! [EvaluationHelper handle]
