#include "../../testInterface.hpp"
#include "multiplyExternalFunction.hpp"

struct TestExtFunctionCall : public TestInterface {
  public:
    NAME("ExtFunctionCall")
    IN(2)
    OUT(1)
    POINTS(1) = {{2.0, 3.0}};

    template<typename Number>
    static void func(Number* x, Number* y) {
      Number w = MultiplyExternalFunction<Number>::create(x[0], x[1]);
      y[0] = w * w;
    }
};
