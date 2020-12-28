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

#include <toolDefines.h>

#include <iostream>
#include <vector>

#include "../output.hpp"

struct ReverseDriverBase {

    virtual ~ReverseDriverBase() {}

    virtual NUMBER::GradientValue& getGradient(NUMBER& number)  = 0;

    virtual void evaluate() = 0;

    virtual void doPreEvaluate() {}

    virtual void doLoopCleanup() {}

    void run() {

      using GT = codi::GradientValueTraits<Gradient>;
      constexpr size_t gradDim = GT::getVectorSize();

      int evalPoints = getEvalPointsCount();
      int inputs = getInputCount();
      int outputs = getOutputCount();
      NUMBER* x = new NUMBER[inputs];
      NUMBER* y = new NUMBER[outputs];

      codi::Jacobian<std::vector<double>> jac(outputs, inputs);

      NUMBER::TapeType& tape = NUMBER::getGlobalTape();
      tape.resize(10000, 10000);
      tape.setExternalFunctionChunkSize(1000);
    #if PRIMAL
      tape.setConstantDataSize(10000);
    #endif

      for(int curPoint = 0; curPoint < evalPoints; ++curPoint) {
        std::cout << "Point " << curPoint << " : {";

        for(int i = 0; i < inputs; ++i) {
          if(i != 0) {
            std::cout << ", ";
          }
          double val = getEvalPoint(curPoint, i);
          std::cout << val;

          x[i] = (NUMBER)(val);
        }
        std::cout << "}\n";

        for(int i = 0; i < outputs; ++i) {
          y[i] = 0.0;
        }

        int runs = outputs / gradDim;
        if(outputs % gradDim != 0) {
          runs += 1;
        }
        for(int curOut = 0; curOut < runs; ++curOut) {
          size_t curSize = gradDim;
          if((curOut + 1) * gradDim  > (size_t)outputs) {
            curSize = outputs % gradDim;
          }

          tape.setActive();

          for(int i = 0; i < inputs; ++i) {
            tape.registerInput(x[i]);
          }

          func(x, y);

          for(int i = 0; i < outputs; ++i) {
            tape.registerOutput(y[i]);
          }

          tape.setPassive();

          doPreEvaluate();

          for(size_t curDim = 0; curDim < curSize; ++curDim) {
            if(y[curOut * gradDim + curDim].isActive()) {
              GT::at(getGradient(y[curOut * gradDim + curDim]), curDim) = 1.0;
            }
          }

          evaluate();

          for(size_t curDim = 0; curDim < curSize; ++curDim) {
            for(int curIn = 0; curIn < inputs; ++curIn) {
#if SECOND_ORDER
              jac(curOut * gradDim + curDim, curIn) = getGradient(x[curIn]).getValue();
#else
              jac(curOut * gradDim + curDim, curIn) = GT::at(getGradient(x[curIn]), curDim);
#endif
            }
          }

          tape.reset();

          doLoopCleanup();
        }

        writeOutputJacobian(jac);
      }
    }
};
