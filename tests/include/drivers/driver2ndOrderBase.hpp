#pragma once

#include <codi/tools/data/hessian.hpp>

#include "../output.hpp"
#include "driverBase.hpp"

template<typename _Number>
struct Driver2ndOrderBase : public DriverBase<_Number> {
  public:

    using Number = CODI_DECLARE_DEFAULT(_Number, double);
    using Base = DriverBase<Number>;

    virtual void evaluateHessian(TestInfo<Number>& info, Number* x, size_t inputs, Number* y, size_t outputs, codi::Hessian<double>& hes) = 0;

  public:

    Driver2ndOrderBase(std::string const& name) : Base(name) {}

    void runTest(TestInfo<Number>& info, FILE* out) {

      TestInterface* test = info.test;

      int evalPoints = test->getEvalPointsCount();
      int inputs = test->getInputCount();
      int outputs = test->getOutputCount();
      Number* x = new Number[inputs];
      Number* y = new Number[outputs];

      codi::Hessian<double> hes(outputs, inputs);

      for (int curPoint = 0; curPoint < evalPoints; ++curPoint) {

        Base::prepare(x, y, curPoint, test, out);

        evaluateHessian(info, x, inputs, y, outputs, hes);

        writeOutputHessian(out, hes);
      }

      delete [] x;
      delete [] y;
    }
};
