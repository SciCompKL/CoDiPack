#pragma once

#include <codi.hpp>
#include <codi/tools/data/jacobian.hpp>

#include "../driverInterface.hpp"
#include "../../tests/allTests.hpp"
#include "../../output.hpp"

#include DRIVER_TESTS_INC

template<typename Number>
void createTests(TestVector<Number>& tests) {}

template<typename Number, typename Test, typename... Args>
void createTests(TestVector<Number>& tests) {
  tests.push_back(TestInfo<Number>(new Test(), Test::template func<Number>));

  createTests<Number, Args...>(tests);
}


struct CoDiForward1stOrder : public DriverInterface<CODI_TYPE> {
  public:

    using Number = CODI_TYPE;

    std::string getName() {
      return CODI_TO_STRING(CODI_TYPE_NAME);
    }

    DriverOrder getOrder() {
      return DriverOrder::Deriv1st;
    }

    TestVector<Number> getTestInfos() {
      TestVector<Number> testInfos;

      createTests<Number, DRIVER_TESTS>(testInfos);

      return testInfos;
    }

    void runTest(TestInfo<Number>& info, FILE* out) {

      TestInterface* test = info.test;

      //using GT = codi::GradientValueTraits<Gradient>;
      constexpr size_t gradDim = 1; // GT::getVectorSize();

      int evalPoints = test->getEvalPointsCount();
      int inputs = test->getInputCount();
      int outputs = test->getOutputCount();
      Number* x = new Number[inputs];
      Number* y = new Number[outputs];

      codi::Jacobian<double> jac(outputs, inputs);

      using Gradient = typename Number::Gradient;

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

        int runs = inputs / gradDim;
        if(inputs % gradDim != 0) {
          runs += 1;
        }
        for(int curIn = 0; curIn < runs; ++curIn) {
          size_t curSize = gradDim;
          if((curIn + 1) * gradDim  > (size_t)inputs) {
            curSize = inputs % gradDim;
          }
          for(size_t curDim = 0; curDim < curSize; ++curDim) {
            //GT::at(x[curIn * gradDim + curDim].gradient(), curDim) = 1.0;
            x[curIn * gradDim + curDim].gradient() = 1.0;
          }

          for(int i = 0; i < outputs; ++i) {
            y[i].setGradient(Gradient());
          }

          info.func(x, y);

          for(size_t curDim = 0; curDim < curSize; ++curDim) {
            for(int curOut = 0; curOut < outputs; ++curOut) {
              //jac(curOut, curIn * gradDim + curDim) = GT::at(y[curOut].getGradient(), curDim);
              jac(curOut, curIn * gradDim + curDim) = y[curOut].getGradient();
            }
          }

          for(size_t curDim = 0; curDim < curSize; ++curDim) {
            x[curIn * gradDim + curDim].setGradient(Gradient());
          }
        }

        writeOutputJacobian(out, jac);
      }

      delete [] x;
      delete [] y;

    }
};
