/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2020 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * CoDiPack is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU
 * General Public License along with CoDiPack.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *     Max Sagebaum
 *     Tim Albring
 *     Johannes Blühdorn
 */
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

template<typename Vec>
void printVector(std::string const& name, Vec const& v, size_t length, size_t offset) {
  std::cout << "Vector " << name << ": {";
  for(size_t i = 0; i < length; i += 1) {
    if(i != 0) {
      std::cout << ", ";
    }
    std::cout << v[offset + i];
  }
  std::cout << "}" << std::endl;
}

template<typename Vec>
void printVectorDim(std::string const& name, Vec const& v, size_t length, size_t offset, size_t dim) {
  std::cout << "Vector " << name << ": {";
  for(size_t i = 0; i < length; i += 1) {
    if(i != 0) {
      std::cout << ", ";
    }
    std::cout << v[offset + i][dim];
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
      std::cerr << "  Mode  1: separate evaluation of Hessian and Jacobian" << std::endl;
      std::cerr << "  Mode  2: combined evaluation of Hessian and Jacobian" << std::endl;

      exit(-1);
    }
  }

  using TH = codi::TapeHelper<codi::HessianComputationType>;
  TH th;

  const size_t n = 10;

  std::vector<codi::HessianComputationType> a(n);
  std::vector<codi::HessianComputationType> b(n);
  for(size_t i = 0; i < n; i += 1) {
    a[i] = i;
    b[i] = pow(-1, i);
  }

  th.startRecording();
  for(size_t i = 0; i < n; i += 1) {
    th.registerInput(a[i]);
  }
  for(size_t i = 0; i < n; i += 1) {
    th.registerInput(b[i]);
  }

  codi::HessianComputationType alpha, aNorm, bNorm;
  dotWithNorms(a.data(), b.data(), n, alpha, aNorm, bNorm);

  th.registerOutput(alpha);
  th.registerOutput(aNorm);
  th.registerOutput(bNorm);

  th.stopRecording();

  TH::JacobianType& jac = th.createJacobian();
  TH::HessianType& hes = th.createHessian();

  if(1 == mode) {
    th.evalJacobian(jac);
    th.evalHessian(hes);
  } else {
    th.evalHessian(hes, jac);
  }

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

  // Evaluate gradient
  TH::GradientValue* x_b = th.createGradientVectorInput();
  TH::GradientValue* y_b = th.createGradientVectorOutput();

  y_b[0] = {1.0, 0.0, 0.0, 0.0};
  y_b[1] = {0.0, 1.0, 0.0, 0.0};
  y_b[2] = {0.0, 0.0, 1.0, 0.0};

  th.evalReverse(y_b, x_b);

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
