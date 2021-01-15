#include "../../testInterface.hpp"

struct TestPreaccumulationForward : public TestInterface {
  public:
    NAME("PreaccumulationForward")
    IN(2)
    OUT(4)
    POINTS(1) = {{1.0, 0.5}};

    template<typename Number>
    static void evalFunc(Number* x, Number* y) {
      y[0] = x[0];
      y[1] = x[1];
      for (int i = 0; i < 1000; ++i) {
        Number xTemp = y[0];
        Number yTemp = y[1];

        Number xSqr = xTemp * xTemp;
        Number ySqr = yTemp * yTemp;

        y[0] = xSqr - ySqr - 0.65;
        y[1] = 2.0 * yTemp * xTemp;
      }

      y[2] = x[0] * x[0];
      y[3] = x[1] * x[1];
    }

    template<typename Number>
    static void func(Number* x, Number* y) {

      codi::PreaccumulationHelper<Number> ph;

      ph.start(x[0], x[1]);

      evalFunc(x, y);

      ph.finish(false, y[0], y[1], y[2], y[3]);
    }
};
