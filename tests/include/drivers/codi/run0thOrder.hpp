#pragma once

#include <codi.hpp>

#include "../driver0thOrderBase.hpp"

#include DRIVER_TESTS_INC

struct CoDi0thOrder : public Driver0thOrderBase<CODI_TYPE> {
  public:

    using Number = CODI_TYPE;

    using Base = Driver0thOrderBase<Number>;

    CoDi0thOrder() : Base(CODI_TO_STRING(CODI_TYPE_NAME)) {}

    void createAllTests(TestVector<Number>& tests) {
      createTests<Number, DRIVER_TESTS>(tests);
    }

    void evaluatePrimal(TestInfo<Number>& info, Number* x, size_t inputs, Number* y, size_t outputs,
                        std::vector<double>& primals) {
      codi::CODI_UNUSED(inputs);

      info.func(x, y);

      for (size_t curOut = 0; curOut < outputs; curOut += 1) {
        primals[curOut] = codi::RealTraits::getPassiveValue(y[curOut]);
      }
    }
};
