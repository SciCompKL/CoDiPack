#pragma once

#include <codi/tools/data/jacobian.hpp>

#include "../output.hpp"
#include "driverBase.hpp"

template<typename _Number>
struct Driver1stOrderBase : public DriverBase<_Number> {
  public:

    using Number = CODI_DECLARE_DEFAULT(_Number, double);
    using Base = DriverBase<Number>;

    virtual void evaluateJacobian(TestInfo<Number>& info, Number* x, size_t inputs, Number* y, size_t outputs, codi::Jacobian<double>& jac) = 0;

  public:

    Driver1stOrderBase(std::string const& name) : Base(name) {}

    void runTest(TestInfo<Number>& info, FILE* out) {

      TestInterface* test = info.test;

      int evalPoints = test->getEvalPointsCount();
      int inputs = test->getInputCount();
      int outputs = test->getOutputCount();
      Number* x = new Number[inputs];
      Number* y = new Number[outputs];

      codi::Jacobian<double> jac(outputs, inputs);

      for(int curPoint = 0; curPoint < evalPoints; ++curPoint) {

        Base::prepare(x, y, curPoint, test, out);

        evaluateJacobian(info, x, inputs, y, outputs, jac);

        writeOutputJacobian(out, jac);
      }

      delete [] x;
      delete [] y;

    }
};
