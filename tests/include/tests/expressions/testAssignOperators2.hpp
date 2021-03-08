#include "../../testInterface.hpp"

struct TestAssignOperators2 : public TestInterface {
  public:
    NAME("AssignOperators2")
    IN(2)
    OUT(3)
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
      y[0] = x[0];
      y[0] /= x[1];  // R x (R \ {0})
      y[1] = 5.00;
      y[1] /= x[1];  // R x (R \ {0})
      y[2] = x[0];
      y[2] /= 5.00;  // R x (R \ {0})
    }
};
