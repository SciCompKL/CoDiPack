#include "../../testInterface.hpp"

struct TestCopy : public TestInterface {
  public:
    NAME("Copy")
    IN(1)
    OUT(1)
    POINTS(1) = {{1.0}};

    template<typename Number>
    static void func(Number* x, Number* y) {
      y[0] = x[0];
    }
};
