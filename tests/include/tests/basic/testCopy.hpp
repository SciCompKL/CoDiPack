#include "../../testInterface.hpp"

template<typename Number>
struct TestCopy : public TestInterface<Number> {

    NAME("Copy")
    IN(1)
    OUT(1)
    POINTS(1) = {{1.0}};

    void func(Number* x, Number* y) {
      y[0] = x[0];
    }
};
