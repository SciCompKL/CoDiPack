#include "../../testInterface.hpp"

struct TestOutput : public TestInterface {
  public:
    NAME("Output")
    IN(2)
    OUT(2)
    POINTS(1) = {{1.0, 2.0}};

    template<typename Number>
    static void func(Number* x, Number* y) {
      Number t1 = x[0] * x[1];
      y[0] = t1;
      y[1] = t1;
    }
};
