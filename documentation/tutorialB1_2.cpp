/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
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
