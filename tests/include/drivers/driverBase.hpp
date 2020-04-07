#pragma once

#include <codi/tools/data/jacobian.hpp>

#include "../output.hpp"
#include "../tests/allTests.hpp"
#include "driverInterface.hpp"

template<typename _Number>
struct DriverBase : public DriverInterface<_Number> {
  public:

    using Number = DECLARE_DEFAULT(_Number, double);

    virtual void createAllTests(TestVector<Number>& tests) = 0;
    virtual void evaluateJacobian(TestInfo<Number>& info, Number* x, size_t inputs, Number* y, size_t outputs, codi::Jacobian<double>& jac) = 0;

  private:

    std::string name;

  public:

    DriverBase(std::string const& name) :
      name(name) {}

    std::string getName() {
      return name;
    }

    TestVector<Number> getTestInfos() {
      TestVector<Number> testInfos;

      createAllTests(testInfos);

      return testInfos;
    }

    void runTest(TestInfo<Number>& info, FILE* out) {

      TestInterface* test = info.test;

      int evalPoints = test->getEvalPointsCount();
      int inputs = test->getInputCount();
      int outputs = test->getOutputCount();
      Number* x = new Number[inputs];
      Number* y = new Number[outputs];

      codi::Jacobian<double> jac(outputs, inputs);

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

        evaluateJacobian(info, x, inputs, y, outputs, jac);

        writeOutputJacobian(out, jac);
      }

      delete [] x;
      delete [] y;

    }

  protected:
    template<typename Number>
    void createTests(TestVector<Number>& tests) {
      (void)tests;
    }

    template<typename Number, typename Test, typename... Args>
    void createTests(TestVector<Number>& tests) {
      tests.push_back(TestInfo<Number>(new Test(), Test::template func<Number>));

      createTests<Number, Args...>(tests);
    }

};
