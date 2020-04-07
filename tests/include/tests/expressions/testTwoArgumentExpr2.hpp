#include "../../testInterface.hpp"

struct TestTwoArgumentExpr2 : public TestInterface {
  public:
    NAME("TwoArgumentExpr2")
    IN(2)
    OUT(12)
    POINTS(18) =
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
    };
    
    template<typename Number>
    static void func(Number* x, Number* y) {
      y[0] = x[0] / x[1];  // R x (R \ {0})
      y[1] = 5.00 / x[1];  // R x (R \ {0})
      y[2] = x[0] / 5.00;  // R x (R \ {0})
      y[3] =   pow(x[0], x[1]);  // R x R
      y[4] =   pow(5.00, x[1]);  // R x R
      y[5] =   pow(x[0], 5.00);  // R x R
      y[6] = atan2(x[0], x[1]);  // R x R \ {0, 0}
      y[7] = atan2(5.00, x[1]);  // R x R \ {0, 0}
      y[8] = atan2(x[0], 5.00);  // R x R \ {0, 0}
      y[9] = copysign(x[0], x[1]); // R x R
      y[10]= copysign(5.00, x[1]); // R x R
      y[11]= copysign(x[0], 5.00); // R x R
    }
};
