/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#include <codi.hpp>

#include "../include/tapeReadWriteBase.hpp"

void openFile(FILE*& fileHandle, std::string const& name, std::string const& mode) {
  fileHandle = fopen(name.c_str(), mode.c_str());
  if (nullptr == fileHandle) {
    CODI_EXCEPTION("Could not open file %s", name.c_str());
  }
}

template<typename Tape, typename Id, typename Grad>
void evalTape(Tape& tape, std::vector<Id> const& x_id, std::vector<Id> const& y_id, std::vector<Grad> const& seed,
              std::vector<Grad>& grad) {
  for (size_t i = 0; i < x_id.size(); i += 1) {
    tape.gradient(y_id[i]) = seed[i];
  }

  tape.evaluate();

  for (size_t i = 0; i < x_id.size(); i += 1) {
    grad[i] = tape.gradient(x_id[i]);
  }
}

template<typename Grad>
void writeGrad(std::vector<Grad> const& grad, FILE*& resultHandle) {
  for (size_t i = 0; i < grad.size(); i += 1) {
    fprintf(resultHandle, "%zu: %0.12e\n", i, grad[i]);
  }
}
template<typename Real>
void runTest(std::string const& name, FILE*& resultHandle) {
  fprintf(resultHandle, "Running %s:\n", name.c_str());

  using Id = typename Real::Identifier;
  using Tape = typename Real::Tape;

  size_t const n = 5;

  std::vector<Real> x(n);
  std::vector<Real> y(n);

  std::vector<Id> x_id(n);
  std::vector<Id> y_id(n);

  std::vector<typename Real::Real> x_grad(n);
  std::vector<typename Real::Real> y_grad(n);

  Tape& tape = Real::getTape();
  tape.setActive();

  func(tape, x, y, x_id, y_id);

  tape.setPassive();

  for (size_t i = 0; i < n; i += 1) {
    y_grad[i] = 1.0;
  }

  evalTape(tape, x_id, y_id, y_grad, x_grad);

  writeGrad(x_grad, resultHandle);

  tape.writeTape(codi::createWriter<Real>(name + "Text.txt", x_id, y_id, codi::FileType::Text));
  tape.writeTape(codi::createWriter<Real>(name + "Binary.dat", x_id, y_id, codi::FileType::Binary));
  tape.writeTape(codi::createWriter<Real>(name + "Graph.dot", x_id, y_id, codi::FileType::Graph));
  if (codi::TapeTraits::isPrimalValueTape<Tape>) {
    tape.writeTape(codi::createWriter<Real>(name + "Math.txt", x_id, y_id, codi::FileType::Math));
  }

  tape.resetHard();
}

int main(int nargs, char** args) {
  FILE* resultHandle;
  openFile(resultHandle, "tapeWrite.out", "w");

  runTest<codi::RealReverse>("jacobian_linear", resultHandle);
  runTest<codi::RealReverseIndexGen<double, double, codi::ReuseIndexManager<int>>>("jacobian_reuse", resultHandle);
  runTest<codi::RealReverseIndex>("jacobian_multiuse", resultHandle);

  runTest<codi::RealReversePrimal>("primal_linear", resultHandle);
  runTest<codi::RealReversePrimalIndexGen<double, double, codi::ReuseIndexManager<int>>>("primal_reuse", resultHandle);
  runTest<codi::RealReversePrimalIndex>("primal_multiuse", resultHandle);
}
