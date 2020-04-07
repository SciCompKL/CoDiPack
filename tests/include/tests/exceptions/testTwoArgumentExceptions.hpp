#include "../../testInterface.hpp"

struct TestTwoArgumentExceptions : public TestInterface {
  public:
    NAME("TwoArgumentExceptions")
    IN(2)
    OUT(9)
    POINTS(4) =
    {
      {0.0,   0.0},
      {1.0,   0.0},
      {0.0,  -5.0},
      {-1.0,  0.5}
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
    }
};
