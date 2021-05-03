#include "../../testInterface.hpp"

struct TestOneArgumentExpr2 : public TestInterface {
  public:
    NAME("OneArgumentExpr2")
    IN(1)
    OUT(4)
    POINTS(20) =  // clang-format off
    {
      {  0.5000},
      {  1.0000},
      {  1.5000},
      {  2.0000},
      {  2.5000},
      {  3.0000},
      {  3.5000},
      {  4.0000},
      {  4.5000},
      {  5.0000},
      {  5.5000},
      {  6.0000},
      {  6.5000},
      {  7.0000},
      {  7.5000},
      {  8.0000},
      {  8.5000},
      {  9.0000},
      {  9.5000},
      { 10.0000}
    };  // clang-format on

    template<typename Number>
    static void func(Number* x, Number* y) {
      y[0] = log(x[0]);     // (0, inf)
      y[1] = log10(x[0]);   // (0, inf)
      y[2] = sqrt(x[0]);    // [0, inf)
      y[3] = tgamma(x[0]);  // R currently only defined for positive arguments
    }
};