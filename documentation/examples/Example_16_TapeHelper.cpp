//! [Example 16 - Tape helper]
#include <iostream>

#include <codi.hpp>

#include "outputHelpers.hpp"

using Real = codi::HessianComputationType;

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

int main(int nargs, char** args) {

  int mode = 1;

  if(2 <= nargs) {
    mode = std::stoi(args[1]);

    if(mode < 1 || 2 < mode) {
      std::cerr << "Error: Please enter a mode from 1 to 2, it was '" << mode << "'." << std::endl;
      std::cerr << "  Mode  1: separate evaluation of Hessian and Jacobian" << std::endl;
      std::cerr << "  Mode  2: combined evaluation of Hessian and Jacobian" << std::endl;

      exit(-1);
    }
  }

//! [Tape recording]
  using TH = codi::TapeHelper<Real>;
  TH th;                                 // Step 1: Create the tape helper

  const size_t n = 10;

  std::vector<Real> a(n);
  std::vector<Real> b(n);
  for(size_t i = 0; i < n; i += 1) {
    a[i] = i;
    b[i] = pow(-1, i);
  }

  th.startRecording();                  // Step 2: Start the recording
  for(size_t i = 0; i < n; i += 1) {
    th.registerInput(a[i]);             // Step 3: Register the inputs
  }
  for(size_t i = 0; i < n; i += 1) {
    th.registerInput(b[i]);
  }

  // Step 4: Perform the computation
  Real alpha, aNorm, bNorm;
  dotWithNorms(a.data(), b.data(), n, alpha, aNorm, bNorm);

  // Step 5: Register the outputs
  th.registerOutput(alpha);
  th.registerOutput(aNorm);
  th.registerOutput(bNorm);

  th.stopRecording();                    // Step 6: Stop the recording
//! [Tape recording]

//! [Hessian evaluation]
  TH::JacobianType& jac = th.createJacobian();
  TH::HessianType& hes = th.createHessian();

  if(1 == mode) {
    th.evalJacobian(jac);
    th.evalHessian(hes);
  } else {
    th.evalHessian(hes, jac);
  }
//! [Hessian evaluation]

  printVector("a", a, n, 0);
  printVector("b", b, n, 0);
  std::cout << std::endl;
  printJacCol("Jacobian with respect to alpha: ", jac, 0);
  printJacCol("Jacobian with respect to aNorm: ", jac, 1);
  printJacCol("Jacobian with respect to bNorm: ", jac, 2);
  std::cout << std::endl;
  printHesForOutput("Hessian with respect to alpha: ", hes, 0);
  printHesForOutput("Hessian with respect to aNorm: ", hes, 1);
  printHesForOutput("Hessian with respect to bNorm: ", hes, 2);


  // Evaluate at different position
  TH::Real* x = th.createPrimalVectorInput();
  TH::Real* y = th.createPrimalVectorOutput();

  for(size_t i = 0; i < n; i += 1) {
    x[0 + i] = i * i;
    x[n + i] = pow(-1, i + 1);
  }

  if(1 == mode) {
    th.evalJacobianAt(x, jac, y);
    // Jacobian evaluation already shifted the point for the evaluation. No second ...At call is required here
    th.evalHessian(hes);
  } else {
    th.evalHessianAt(x, hes, y, jac);
  }

  std::cout << std::endl;
  std::cout << "Reevaluation at new location:" << std::endl;
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

  // Perform a regular AD reverse mode interpretation of the tape.
//! [Reverse evaluation]
  TH::Gradient* x_b = th.createGradientVectorInput();
  TH::Gradient* y_b = th.createGradientVectorOutput();

  y_b[0] = {1.0, 0.0, 0.0, 0.0};
  y_b[1] = {0.0, 1.0, 0.0, 0.0};
  y_b[2] = {0.0, 0.0, 1.0, 0.0};

  th.evalReverse(y_b, x_b);
//! [Reverse evaluation]

  std::cout << "Reverse evaluation for alpha_b:" << std::endl;
  printVectorDim("a_b", x_b, n, 0, 0);
  printVectorDim("b_b", x_b, n, n, 0);
  std::cout << std::endl;
  std::cout << "Reverse evaluation for aNorm_b:" << std::endl;
  printVectorDim("a_b", x_b, n, 0, 1);
  printVectorDim("b_b", x_b, n, n, 1);
  std::cout << std::endl;
  std::cout << "Reverse evaluation for bNorm_b:" << std::endl;
  printVectorDim("a_b", x_b, n, 0, 2);
  printVectorDim("b_b", x_b, n, n, 2);

  // Clean up vectors
  th.deleteGradientVector(x_b);
  th.deleteGradientVector(y_b);

  th.deletePrimalVector(x);
  th.deletePrimalVector(y);

  th.deleteJacobian(jac);
  th.deleteHessian(hes);

  return 0;
}
//! [Example 16 - Tape helper]
