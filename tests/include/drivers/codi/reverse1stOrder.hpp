#pragma once

#include <stdio.h>

#include <codi.hpp>
#include <codi/tools/data/jacobian.hpp>

#include "../driverInterface.hpp"
#include "../../tests/listsTests.hpp"
#include "../../output.hpp"


struct CoDiReverse1stOrder : public DriverInterface<CODI_TYPE> {
  public:

    using Number = CODI_TYPE;

    std::string getName() {
      return CODI_TO_STRING(CODI_TYPE_NAME);
    }

    DriverOrder getOrder() {
      return DriverOrder::Deriv1st;
    }

    TestVector<Number> getTests() {
      TestVector<Number> tests;

      listTestAll(tests);

      return tests;
    }

    void runTest(TestInterface<Number>* test, FILE* out) {

      //using GT = codi::GradientValueTraits<Gradient>;
      constexpr size_t gradDim = 1; //GT::getVectorSize();

      int evalPoints = test->getEvalPointsCount();
      int inputs = test->getInputCount();
      int outputs = test->getOutputCount();
      Number* x = new Number[inputs];
      Number* y = new Number[outputs];

      codi::Jacobian<double> jac(outputs, inputs);

      Number::Tape& tape = Number::getGlobalTape();

      for(int curPoint = 0; curPoint < evalPoints; ++curPoint) {
        fprintf(out, "Point %d : {", curPoint);

        for(int i = 0; i < inputs; ++i) {
          if(i != 0) {
            fprintf(out, ", ");
          }
          double val = test->getEvalPoint(curPoint, i);
          fprintf(out, "%0.0f", val);

          x[i] = (Number)(val);
        }
        fprintf(out, "}\n");

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

          test->func(x, y);

          for(int i = 0; i < outputs; ++i) {
            tape.registerOutput(y[i]);
          }

          tape.setPassive();

          for(size_t curDim = 0; curDim < curSize; ++curDim) {
            if(tape.isIdentifierActive(y[curOut * gradDim + curDim].getIdentifier())) {
              y[curOut * gradDim + curDim].gradient() = 1.0;
            }
          }

          tape.evaluate();

          for(size_t curDim = 0; curDim < curSize; ++curDim) {
            for(int curIn = 0; curIn < inputs; ++curIn) {
              //jac(curOut * gradDim + curDim, curIn) = GT::at(getGradient(x[curIn]), curDim);
              jac(curOut * gradDim + curDim, curIn) = x[curIn].getGradient();
            }
          }

          tape.reset();
        }

        writeOutputJacobian(out, jac);
      }

      delete [] x;
      delete [] y;

    }
};
