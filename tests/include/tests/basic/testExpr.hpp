#include "../../testInterface.hpp"

struct TestExpr : public TestInterface {
    NAME("Expr")
    IN(2)
    OUT(1)
    POINTS(1) = {{1.0, 2.0}};

    template<typename Number>
    static void func(Number* x, Number* y) {
      y[0] = x[0] + x[1];
    }
};
