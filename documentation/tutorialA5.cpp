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
  EH::evalHessian(wrapDotWithNorms, x, 3, hes);

  printVector("a", x, n, 0);
  printVector("b", x, n, n);
  std::cout << std::endl;
  printHesForOutput("Hessian with respect to alpha: ", hesData, 0, 3, xSize);
  printHesForOutput("Hessian with respect to aNorm: ", hesData, 1, 3, xSize);
  printHesForOutput("Hessian with respect to bNorm: ", hesData, 2, 3, xSize);

  delete [] hesData;

  return 0;
}
