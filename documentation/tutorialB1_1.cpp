/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
