#pragma once

#include "../../../testInterface.hpp"

#include "multiplyExternalFunctionHelper.hpp"

struct TestExternalFunctionHelper : public TestInterface {
  public:
    NAME("ExternalFunctionHelper")
    IN(2)
    OUT(1)
    POINTS(1) = {{2.0, 3.0}};

    static int constexpr ITER = 5;

    template<typename Number>
    static void func(Number* x, Number* y) {
      Number w[ITER];

      w[0] = x[0];
      for(int i = 1; i < ITER; ++i) {
        w[i] = MultiplyExternalFunctionHelper<Number>::create(x[1], w[i-1], false);
      }

      y[0] = w[ITER - 1]*w[ITER - 1];
    }
};
