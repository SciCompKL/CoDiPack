#pragma once

#include "../driver2ndOrderBase.hpp"

#include DRIVER_TESTS_INC

struct CoDiEvalHelper2ndOrder : public Driver2ndOrderBase<CODI_TYPE> {
  public:

    using Number = CODI_TYPE;

    using Base = Driver2ndOrderBase;

    using Gradient = Number::Gradient;

    CoDiEvalHelper2ndOrder() : Base(CODI_TO_STRING(CODI_TYPE_NAME)) {}

    void createAllTests(TestVector<Number>& tests) {
      createTests<Number, DRIVER_TESTS>(tests);
    }

    void evaluateHessian(TestInfo<Number>& info, Number* x, size_t inputs, Number* /*y*/, size_t outputs, codi::Hessian<double>& hes) {

      std::vector<double> xVec(inputs);

      for (size_t i = 0; i < inputs; ++i) {
        xVec[i] = codi::RealTraits::getPassiveValue(x[i]);
      }

      auto evalFunc = [&](std::vector<Number>& x, std::vector<Number>& y) {
        info.func(x.data(), y.data());
      };

      auto handle = codi::EvaluationHelper::template createHandle<Number>(evalFunc, outputs, inputs);

      codi::EvaluationHelper::evalHandleHessian(handle, xVec, hes);

      // evaluate a second time to force at least one tape reset.
      codi::EvaluationHelper::evalHandleHessian(handle, xVec, hes);
    }
};
