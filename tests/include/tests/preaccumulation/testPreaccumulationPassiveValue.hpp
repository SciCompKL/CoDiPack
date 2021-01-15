#include "../../testInterface.hpp"

struct TestPreaccumulationPassiveValue : public TestInterface {
  public:
    NAME("PreaccumulationPassiveValue")
    IN(2)
    OUT(2)
    POINTS(1) = {{1.0, 0.5}};

    template<typename Number>
    static void evalFunc(Number* x, Number* y) {
      y[0] = codi::RealTraits::getPassiveValue(x[0]); // kill x dependency
      y[1] = x[1];

      Number two = 2.0;
      Number zeroSix = 0.65;
      for (int i = 0; i < 5; ++i) {
        Number xTemp = y[0];
        Number yTemp = y[1];

        y[0] = xTemp * xTemp - yTemp * yTemp - zeroSix;
        y[1] = two * yTemp * xTemp;
      }
    }

    template<typename Number>
    static void func(Number* x, Number* y) {

      codi::PreaccumulationHelper<Number> ph;

      ph.start(x[0], x[1]);

      evalFunc(x, y);

      ph.finish(false, y[0], y[1]);
    }
};
