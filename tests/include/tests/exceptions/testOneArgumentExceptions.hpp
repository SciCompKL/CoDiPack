#define _USE_MATH_DEFINES
#include <cmath>

#include "../../testInterface.hpp"

struct TestOneArgumentExceptions : public TestInterface {
  public:
    NAME("OneArgumentExceptions")
    IN(1)
    OUT(7)
    POINTS(6) =
    {
      {0.5 * M_PI},
      { 0.000},
      { 1.000},
      {-1.000},
      { 2.000},
      {-2.000},
    };
    
    template<typename Number>
    static void func(Number* x, Number* y) {
      y[0] =   tan(x[0]);  // R \ {(1/2 + i) * PI | for all i in Z}
      y[1] =   log(x[0]);  // (0, inf)
      y[2] = log10(x[0]);  // (0, inf)
      y[3] =  sqrt(x[0]);  // [0 , inf)
      y[4] = atanh(x[0]);  // (-1, 1)
      y[5] =  asin(x[0]);  // [-1, 1]
      y[6] =  acos(x[0]);  // [-1, 1]
    }
};
