#include "../../testInterface.hpp"

struct TestTwoArgumentExpr2 : public TestInterface {
  public:
    NAME("TwoArgumentExpr2")
    IN(2)
    OUT(12)
    POINTS(18) =  // clang-format off
    {
      {-10.0,   -10},
      {-10.0,    -5},
      {-10.0,     5},
      {-10.0,    10},
      { -5.0,   -10},
      { -5.0,    -5},
      { -5.0,     5},
      { -5.0,    10},
      {  0.0,     5},
      {  0.0,    10},
      {  5.0,   -10},
      {  5.0,    -5},
      {  5.0,     5},
      {  5.0,    10},
      { 10.0,   -10},
      { 10.0,    -5},
      { 10.0,     5},
      { 10.0,    10}
    };  // clang-format on

    template<typename Number>
    static void func(Number* x, Number* y) {
      y[0] = x[0] / x[1];            // R x (R \ {0})
      y[1] = 5.00 / x[1];            // R x (R \ {0})
      y[2] = x[0] / 5.00;            // R x (R \ {0})
      y[3] = atan2(x[0], x[1]);      // R x R \ {0, 0}
      y[4] = atan2(5.00, x[1]);      // R x R \ {0, 0}
      y[5] = atan2(x[0], 5.00);      // R x R \ {0, 0}
      y[6] = remainder(x[0], x[1]);  // R x (R \ {0})
      y[7] = remainder(5.0, x[1]);   // R x (R \ {0})
      y[8] = remainder(x[0], 5.0);   // R x (R \ {0})
      y[9] = hypot(x[0], x[1]);      // R x R \ {0, 0})
      y[10] = hypot(5.0, x[1]);      // R x R \ {0, 0})
      y[11] = hypot(x[0], 5.0);      // R x R \ {0, 0})
    }
};
