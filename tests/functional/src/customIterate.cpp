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
#include "../include/multLowLevelFunction.hpp"

template<typename Tape>
struct IdStats {
    using Real = typename Tape::Real;
    using Identifier = typename Tape::Identifier;

    using EvalHandle = typename Tape::EvalHandle;

    Tape& tape;
    std::vector<int> idUse;
    Real* primals;

    IdStats(Tape& tape) : tape(tape), idUse() {
      idUse.resize(tape.getIndexManager().getLargestCreatedIndex() + 1);
    }

    void count(Identifier id) {
      idUse[id] +=1;
    }

    void handleStatement(Identifier& lhsIndex, codi::Config::ArgumentSize const& size, Real const* jacobians,
                         Identifier const* rhsIdentifiers) {
      count(lhsIndex);
      for(codi::Config::ArgumentSize i = 0; i < size; i += 1) {
        count(rhsIdentifiers[i]);
      }
    }

    void handleStatement(EvalHandle const& evalHandle, codi::Config::ArgumentSize const& nPassiveValues, size_t& linearAdjointPosition, char* stmtData) {
      using StatementEvaluator = typename Tape::StatementEvaluator;

      codi::WriteInfo writeInfo;
      StatementEvaluator::template call<codi::StatementCall::WriteInformation, Tape>(
          evalHandle, writeInfo, primals, nPassiveValues, stmtData);

      StatementEvaluator::template call<codi::StatementCall::IterateInputs, Tape>(evalHandle, linearAdjointPosition, reinterpret_cast<void (*)(int*, void*)>(countId), this, nPassiveValues, stmtData);
      StatementEvaluator::template call<codi::StatementCall::IterateOutputs, Tape>(evalHandle, linearAdjointPosition, reinterpret_cast<void (*)(int*, void*)>(countId), this, nPassiveValues, stmtData);

      if(Tape::LinearIndexHandling) {
        linearAdjointPosition += writeInfo.numberOfOutputArguments;
      }
    }

    //void handleStatement(EvalHandle const& evalHandle, Config::ArgumentSize const& nPassiveValues, char* stmtData);
    void handleLowLevelFunction(codi::LowLevelFunctionEntry<Tape, Real, Identifier> const& func, codi::ByteDataView& llfData) {
      func.template call<codi::LowLevelFunctionEntryCallKind::IterateOutputs>(&tape, llfData, reinterpret_cast<void (*)(int*, void*)>(countId), this);
      llfData.reset();
      func.template call<codi::LowLevelFunctionEntryCallKind::IterateInputs>(&tape, llfData, reinterpret_cast<void (*)(int*, void*)>(countId), this);
    }

    static void countId(Identifier* id, IdStats* data) {
      data->count(*id);
    }

    void print(std::ofstream& out) {
      for(size_t i = 0; i < idUse.size(); i += 1) {
        if(idUse[i] != 0) {
          out << i << " " << idUse[i] << "\n";
        }
      }
      out.flush();
    }

    void eval() {
      if constexpr (codi::TapeTraits::isPrimalValueTape<Tape>) {
        primals = tape.getPrimalVector();
      }
      tape.iterateForward(*this);
    }
};

template<typename Real>
void runTest(std::ofstream& out, std::string const& name) {
  out << "Running: " << name << std::endl;
  using Id = typename Real::Identifier;
  using Tape = typename Real::Tape;

  size_t const n = 5;

  std::vector<Real> x(n);
  std::vector<Real> y(n);

  std::vector<Id> x_id(n);
  std::vector<Id> y_id(n);

  Tape& tape = Real::getTape();
  tape.setActive();

  func(tape, x, y, x_id, y_id);

  Real w;
  MultLowLevelFunction<Real>::evalAndStore(x[0], x[1], w);
  y[0] += w;

  tape.setPassive();

  IdStats<Tape> stats(tape);
  stats.eval();

  stats.print(out);

  tape.resetHard();
}

int main(int nargs, char** args) {
  std::ofstream out("run.out");

  runTest<codi::RealReverse>(out, "jacobian_linear");
  runTest<codi::RealReverseIndex>(out, "jacobian_multiuse");

  runTest<codi::RealReversePrimal>(out, "primal_linear");
  runTest<codi::RealReversePrimalIndex>(out, "primal_multiuse");
}
