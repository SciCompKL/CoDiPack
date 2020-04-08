#pragma once

#include <codi.hpp>

#include "../driverPrimalBase.hpp"

#include DRIVER_TESTS_INC

struct CoDiPrimal : public DriverPrimalBase<CODI_TYPE> {
  public:

    using Number = CODI_TYPE;

    using Base = DriverPrimalBase<Number>;

    CoDiPrimal() : Base(CODI_TO_STRING(CODI_TYPE_NAME)) {}

    void createAllTests(TestVector<Number>& tests) {
      createTests<Number, DRIVER_TESTS>(tests);
    }

    void evaluatePrimal(TestInfo<Number>& info, Number* x, size_t inputs, Number* y, size_t outputs, std::vector<double>& primals) {
      codi::CODI_UNUSED(inputs);

      info.func(x, y);

      for(size_t curOut = 0; curOut < outputs; curOut += 1) {
        primals[curOut] = codi::getPassiveValue(y[curOut]);
      }
    }
};
