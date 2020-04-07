#include "../../testInterface.hpp"

struct TestIncrementOperators : public TestInterface {
  public:
    NAME("IncrementOperators")
    IN(1)
    OUT(8)
    POINTS(3) =
    {
      {-1.0},
      { 0.0},
      { 1.0}
    };
    
    template<typename Number>
    static void func(Number* x, Number* y) {
      y[0] = x[0];
      y[1] = ++y[0];
    
      y[2] = x[0];
      y[3] = y[2]++;
    
      y[4] = x[0];
      y[5] = --y[4];
    
      y[6] = x[0];
      y[7] = y[6]--;
    }
};
