#include "../../testInterface.hpp"

struct TestPreaccumulationForwardInvalidAdjoint : public TestInterface {
  public:
    NAME("PreaccumulationForwardInvalidAdjoint")
    IN(2)
    OUT(4)
    POINTS(1) = {{1.0, 0.5}};

    template<typename Number>
    static void evalFunc(Number* x, Number* y) {
      Number temp1 = x[0] * x[1];
      Number temp2 = x[0] / x[1];
      Number temp3 = x[0] + x[1];
      Number temp4 = x[0] - x[1];
      Number temp5 = temp1 * temp3;
      Number temp6 = temp2 * temp4;


      y[0] = temp5 * temp5;
      y[1] = temp6 * temp6;
      y[2] = temp5 * temp5;
      y[3] = temp6 * temp6;
    }

    template<typename Number>
    static void func(Number* x, Number* y) {

      codi::PreaccumulationHelper<Number> ph;

      ph.start(x[0], x[1]);

      evalFunc(x, y);

      ph.finish(false, y[0], y[1], y[2], y[3]);

      Number temp1 = y[0] + y[1];
      Number temp2 = y[2] + y[3];

      y[0] = temp1 + temp2;
      y[1] = temp1 - temp2;
      y[2] = temp1 * temp2;
      y[3] = temp1 / temp2;
    }
};
