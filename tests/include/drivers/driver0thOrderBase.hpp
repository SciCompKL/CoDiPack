#pragma once

#include <vector>

#include "../output.hpp"
#include "driverBase.hpp"

template<typename _Number>
struct Driver0thOrderBase : public DriverBase<_Number> {
  public:

    using Number = CODI_DECLARE_DEFAULT(_Number, double);
    using Base = DriverBase<Number>;

    virtual void evaluatePrimal(TestInfo<Number>& info, Number* x, size_t inputs, Number* y, size_t outputs, std::vector<double>& primals) = 0;

  public:

    Driver0thOrderBase(std::string const& name) : Base(name) {}

    void runTest(TestInfo<Number>& info, FILE* out) {

      TestInterface* test = info.test;

      int evalPoints = test->getEvalPointsCount();
      int inputs = test->getInputCount();
      int outputs = test->getOutputCount();
      Number* x = new Number[inputs];
      Number* y = new Number[outputs];

      std::vector<double> primals(outputs);

      for(int curPoint = 0; curPoint < evalPoints; ++curPoint) {

        Base::prepare(x, y, curPoint, test, out);

        evaluatePrimal(info, x, inputs, y, outputs, primals);

        writeOutputPrimal(out, primals);
      }

      delete [] x;
      delete [] y;

    }
};
