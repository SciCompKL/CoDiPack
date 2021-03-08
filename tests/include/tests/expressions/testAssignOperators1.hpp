#include "../../testInterface.hpp"

struct TestAssignOperators1 : public TestInterface {
  public:
    NAME("AssignOperators1")
    IN(2)
    OUT(9)
    POINTS(25) =  // clang-format off
    {
      {-10.0,   -10},
      {-10.0,    -5},
      {-10.0,     0},
      {-10.0,     5},
      {-10.0,    10},
      { -5.0,   -10},
      { -5.0,    -5},
      { -5.0,     0},
      { -5.0,     5},
      { -5.0,    10},
      {  0.0,   -10},
      {  0.0,    -5},
      {  0.0,     0},
      {  0.0,     5},
      {  0.0,    10},
      {  5.0,   -10},
      {  5.0,    -5},
      {  5.0,     0},
      {  5.0,     5},
      {  5.0,    10},
      { 10.0,   -10},
      { 10.0,    -5},
      { 10.0,     0},
      { 10.0,     5},
      { 10.0,    10}
    };  // clang-format on

    template<typename Number>
    static void func(Number* x, Number* y) {
      y[0] = x[0];
      y[0] += x[1];  // R x R
      y[1] = 5.00;
      y[1] += x[1];  // R x R
      y[2] = x[0];
      y[2] += 5.00;  // R x R
      y[3] = x[0];
      y[3] -= x[1];  // R x R
      y[4] = 5.00;
      y[4] -= x[1];  // R x R
      y[5] = x[0];
      y[5] -= 5.00;  // R x R
      y[6] = x[0];
      y[6] *= x[1];  // R x R
      y[7] = 5.00;
      y[7] *= x[1];  // R x R
      y[8] = x[0];
      y[8] *= 5.00;  // R x R
    }
};
