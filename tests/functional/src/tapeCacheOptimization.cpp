/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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

#include "../../../include/codi/tools/identifierCacheOptimizer.hpp"

#include "../include/tapeReadWriteBase.hpp"
#include "../include/inputLowLevelFunction.hpp"
#include "../include/multLowLevelFunction.hpp"
#include "../include/outputLowLevelFunction.hpp"

bool constexpr isCacheOptimization = true;

template<typename T_Real>
struct Test {
  using Real = CODI_DD(T_Real, codi::RealReverseIndex);
  using Id = typename Real::Identifier;
  using IdVec = std::vector<Id>;

  static std::array<typename Real::Real, 100> buffer;

  #define REAL_IN(name, value) \
    Real name{value}; \
    Real::getTape().registerInput(name); \
    in.push_back(name.getIdentifier());

  #define REAL_VEC_IN(name, value) \
    name = value; \
    Real::getTape().registerInput(name); \
    in.push_back(name.getIdentifier());

  #define REAL_OUT(name, value) \
    Real name{value}; \
    Real::getTape().registerOutput(name); \
    out.push_back(name.getIdentifier());

  #define REAL_VEC_OUT(name, value) \
    name = value; \
    Real::getTape().registerOutput(name); \
    out.push_back(name.getIdentifier());


  static void test2Ops(IdVec& in, IdVec& out) {
    REAL_IN(x, 2.0);

    Real t = x * x;
    Real t2 = x * t;

    REAL_OUT(y, t2);
  }

  static void testOverwrite(IdVec& in, IdVec& out) {
    REAL_IN(x, 2.0);

    Real t = x * x;
    t = x * t;

    REAL_OUT(y, t);
  }

  static void testLLF(IdVec& in, IdVec& out) {
    REAL_IN(x, 2.0);

    Real t;
    MultLowLevelFunction<Real>::evalAndStore(x, x, t);
    t = x * t;

    REAL_OUT(y, t);
  }

  static void testLLFVec2(IdVec& in, IdVec& out) {
    std::array<Real, 2> x;
    REAL_VEC_IN(x[0], 2.0);
    REAL_VEC_IN(x[1], 5.0);

    std::array<Real, 2> t;
    MultLowLevelFunction<Real>::evalAndStore(x.data(), x.data(), t.data(), 2);
    for(int i = 0; i < 2; i += 1) {
      t[i] = x[i] * t[i];
    }

    std::array<Real, 2> y;
    REAL_VEC_OUT(y[0], t[0]);
    REAL_VEC_OUT(y[1], t[1]);
  }

  static void testLLFInputOutput(IdVec& in, IdVec& out) {

    std::array<Real, 2> x;
    REAL_VEC_IN(x[0], 2.0);
    REAL_VEC_IN(x[1], 5.0);

    InputLowLevelFunction<Real>::evalAndStore(x.data(), 2, buffer.data()); // Data "send" here

    std::array<Real, 2> t;
    for(int i = 0; i < 2; i += 1) {
      t[i] = x[i] * x[i];
    }

    std::array<Real, 2> t2;
    OutputLowLevelFunction<Real>::evalAndStore(t2.data(), 2, buffer.data()); // Data "received" here

    for(int i = 0; i < 2; i += 1) {
      t2[i] = t2[i] * t[i];
    }

    std::array<Real, 2> y;
    REAL_VEC_OUT(y[0], t2[0]);
    REAL_VEC_OUT(y[1], t2[1]);
  }

  static void testUnusedInput(IdVec& in, IdVec& out) {
    REAL_IN(x, 2.0);
    REAL_IN(unused, 1000.0);

    Real t = x * x;
    Real t2 = x * t;

    REAL_OUT(y, t2);
  }

  static void testUnusedOutput(IdVec& in, IdVec& out) {
    REAL_IN(x, 2.0);

    Real t = x * x;
    Real t2 = x * t;

    REAL_OUT(y, t2);
    REAL_OUT(unused, 1000.0);
  }

  static void testUnusedIntermediate(IdVec& in, IdVec& out) {
    REAL_IN(x, 2.0);

    Real t = x * x;
    Real t2 = x * t;

    Real unused = t2 * t;

    REAL_OUT(y, t2);
  }

  static void testDuplicatedInput(IdVec& in, IdVec& out) {
    REAL_IN(x, 2.0);
    in.push_back(x.getIdentifier());

    Real t = x * x;
    Real t2 = x * t;

    Real unused = t2 * t;

    REAL_OUT(y, t2);
  }

  static void testPassiveInput(IdVec& in, IdVec& out) {
    REAL_IN(x, 2.0);
    in.push_back(0);

    Real t = x * x;
    Real t2 = x * t;

    Real unused = t2 * t;

    REAL_OUT(y, t2);
  }

  static void testLLFPassiveOutput(IdVec& in, IdVec& out) {

    std::array<Real, 3> x;
    REAL_VEC_IN(x[0], 2.0);
    x[1] =  5.0;
    x[2] =  25.0;

    InputLowLevelFunction<Real>::evalAndStore(x.data(), 3, buffer.data()); // Data "send" here

    std::array<Real, 3> t;
    for(int i = 0; i < 3; i += 1) {
      t[i] = x[i] * x[i];
    }

    std::array<Real, 3> t2;
    OutputLowLevelFunction<Real>::evalAndStore(t2.data(), 3, buffer.data()); // Data "received" here

    for(int i = 0; i < 3; i += 1) {
      t2[i] = t2[i] * t[i];
    }

    std::array<Real, 3> y;
    REAL_VEC_OUT(y[0], t2[0]);
    REAL_VEC_OUT(y[1], t2[1]);
    REAL_VEC_OUT(y[2], t2[2]);
  }


  template<typename Test>
  void runTest(std::ofstream& out, std::string const& name, Test&& test) {
    out << "Running: " << name << std::endl;
    using Tape = typename Real::Tape;

    std::vector<Id> xId(0);
    std::vector<Id> yId(0);
    auto iterX = [&xId](auto&& func) {
      for(Id& cur : xId) {
        func(cur);
      }
    };

    auto iterY = [&yId](auto&& func) {
      for(Id& cur : yId) {
        func(cur);
      }
    };

    Tape& tape = Real::getTape();
    tape.setActive();

    test(xId, yId);

    if(isCacheOptimization) {
      codi::IdentifierCacheOptimizerHotCold<Tape> co{tape};
      co.eval(iterX, iterY);
    }

    for(size_t curY = 0; curY < yId.size(); curY += 1) {
      tape.gradient(yId[curY]) = 1.0;
      tape.evaluate();

      for(size_t curX = 0; curX < xId.size(); curX += 1) {
        out << "d y_" << curY << "/ d x_" << curX << " = " << tape.gradient(xId[curX]) << "\n";
        tape.gradient(xId[curX]) = 0.0;
      }
    }

    tape.writeTape(codi::createWriter<Real>(name + ".txt", xId, yId, codi::FileType::Text));

    tape.resetHard();
  }

  void runAllTests(std::ofstream& out) {
    runTest(out, "test2Ops", test2Ops);
    runTest(out, "testOverwrite", testOverwrite);
    runTest(out, "testLLF", testLLF);
    runTest(out, "testLLFVec2", testLLFVec2);
    runTest(out, "testLLInputOutput", testLLFInputOutput);
    runTest(out, "testUnusedInput", testUnusedInput);
    runTest(out, "testUnusedOutput", testUnusedOutput);
    runTest(out, "testUnusedIntermediate", testUnusedIntermediate);
    runTest(out, "testDuplicatedInput", testDuplicatedInput);
    runTest(out, "testPassiveInput", testPassiveInput);
    runTest(out, "testLLFPassiveOutput", testLLFPassiveOutput);
  }
};

template<typename Real>
std::array<typename Real::Real, 100> Test<Real>::buffer = {};

int main(int nargs, char** args) {
  std::ofstream out("run.out");

  out << "RealReverseIndex:" << std::endl;
  Test<codi::RealReverseIndex> t1;
  t1.runAllTests(out);

  out << "RealReversePrimalIndex:" << std::endl;
  Test<codi::RealReverseIndex> t2;
  t2.runAllTests(out);

  return 0;
}
