#pragma once

#include "../driver0thOrderBase.hpp"

#include DRIVER_TESTS_INC

struct CoDiEvalHelper0thOrder : public Driver0thOrderBase<CODI_TYPE> {
  public:

    using Number = CODI_TYPE;

    using Base = Driver0thOrderBase;

    using Gradient = Number::Gradient;

    CoDiEvalHelper0thOrder() : Base(CODI_TO_STRING(CODI_TYPE_NAME)) {}

    void createAllTests(TestVector<Number>& tests) {
      createTests<Number, DRIVER_TESTS>(tests);
    }

    void evaluatePrimal(TestInfo<Number>& info, Number* x, size_t inputs, Number* /*y*/, size_t outputs, std::vector<double>& primals) {

      std::vector<double> xVec(inputs);

      for (size_t i = 0; i < inputs; ++i) {
        xVec[i] = codi::RealTraits::getPassiveValue(x[i]);
      }

      auto evalFunc = [&](std::vector<Number>& x, std::vector<Number>& y) {
        info.func(x.data(), y.data());
      };

      auto handle = codi::EvaluationHelper::template createHandle<Number>(evalFunc, outputs, inputs);

      codi::EvaluationHelper::evalHandlePrimal(handle, xVec, primals);
    }
};
