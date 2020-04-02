#include "../../testInterface.hpp"

template<typename Number>
struct TestExpr : public TestInterface<Number> {
    NAME("Expr")
    IN(2)
    OUT(1)
    POINTS(1) = {{1.0, 2.0}};

    void func(Number* x, Number* y) {
      y[0] = x[0] + x[1];
    }
};
