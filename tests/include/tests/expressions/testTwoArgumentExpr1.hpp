#include "../../testInterface.hpp"

struct TestTwoArgumentExpr1 : public TestInterface {
  public:
    NAME("TwoArgumentExpr1")
    IN(2)
    OUT(24)
    POINTS(25) =
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
    };
    
    template<typename Number>
    static void func(Number* x, Number* y) {
      y[0] = x[0] + x[1];  // R x R
      y[1] = 5.00 + x[1];  // R x R
      y[2] = x[0] + 5.00;  // R x R
      y[3] = x[0] - x[1];  // R x R
      y[4] = 5.00 - x[1];  // R x R
      y[5] = x[0] - 5.00;  // R x R
      y[6] = x[0] * x[1];  // R x R
      y[7] = 5.00 * x[1];  // R x R
      y[8] = x[0] * 5.00;  // R x R
      y[9]  = min(x[0], x[1]);  // R x R
      y[10] = min(5.00, x[1]);  // R x R
      y[11] = min(x[0], 5.00);  // R x R
      y[12] = max(x[0], x[1]);  // R x R
      y[13] = max(5.00, x[1]);  // R x R
      y[14] = max(x[0], 5.00);  // R x R
      y[15] = fmin(x[0], x[1]); // R x R the results for the points with the same values are reversed here because the other one uses the standard template.
      y[16] = fmin(5.00, x[1]); // R x R
      y[17] = fmin(x[0], 5.00); // R x R
      y[18] = fmax(x[0], x[1]); // R x R the results for the points with the same values are reversed here because the other one uses the standard template.
      y[19] = fmax(5.00, x[1]); // R x R
      y[20] = fmax(x[0], 5.00); // R x R
      y[21] = remainder(x[0], x[1]); // R x R
      y[22] = remainder(5.0, x[1]); // R x R
      y[23] = remainder(x[0], 5.0); // R x R
    }
};
