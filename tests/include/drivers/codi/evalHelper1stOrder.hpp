#pragma once

#include "../driver1stOrderBase.hpp"

#include DRIVER_TESTS_INC

struct CoDiEvalHelper1stOrder : public Driver1stOrderBase<CODI_TYPE> {
  public:

    using Number = CODI_TYPE;

    using Base = Driver1stOrderBase;

    using Gradient = Number::Gradient;

    CoDiEvalHelper1stOrder() : Base(CODI_TO_STRING(CODI_TYPE_NAME)) {}

    void createAllTests(TestVector<Number>& tests) {
      createTests<Number, DRIVER_TESTS>(tests);
    }

    void evaluateJacobian(TestInfo<Number>& info, Number* x, size_t inputs, Number* /*y*/, size_t outputs, codi::Jacobian<double>& jac) {

      std::vector<double> xVec(inputs);

      for (size_t i = 0; i < inputs; ++i) {
        xVec[i] = codi::RealTraits::getPassiveValue(x[i]);
      }

      auto evalFunc = [&](std::vector<Number>& x, std::vector<Number>& y) {
        info.func(x.data(), y.data());
      };

      auto handle = codi::EvaluationHelper::template createHandle<Number>(evalFunc, outputs, inputs);

      codi::EvaluationHelper::evalHandleJacobian(handle, xVec, jac);

      // evaluate a second time to force at least one tape reset.
      codi::EvaluationHelper::evalHandleJacobian(handle, xVec, jac);
    }
};
