#include "../../../testInterface.hpp"

struct TestPreaccumulationLargeStatement : public TestInterface {
  public:
    NAME("PreaccumulationLargeStatement")
    IN(2)
    OUT(2)
    POINTS(1) = {{1.0, 0.5}};

    template<typename Number>
    static void evalFunc(Number* x, Number* y, size_t size) {
      y[0] = x[0];
      y[1] = x[0];
      for (size_t i = 1; i < size; ++i) {
        y[0] += x[i];
        y[1] = max(y[1], x[i]);
      }
    }

    template<typename Number>
    static void func(Number* x, Number* y) {
      codi::PreaccumulationHelper<Number> ph;

      size_t const size = 256 * 3;
      Number intermediate[size];

      for (size_t i = 0; i < size; ++i) {
        intermediate[i] = x[0] * (double)i + x[1];
      }

      ph.start();
      for (size_t i = 0; i < size; ++i) {
        ph.addInput(intermediate[i]);
      }

      evalFunc(intermediate, y, size);

      ph.addOutput(y[0]);
      ph.addOutput(y[1]);
      ph.finish(false);
    }
};
