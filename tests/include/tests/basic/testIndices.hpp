#include "../../testInterface.hpp"

struct TestIndices : public TestInterface {
  public:
    NAME("Indices")
    IN(2)
    OUT(1)
    POINTS(1) = {{1.0, 2.0}};
    
    template<typename Number>
    static void func(Number* x, Number* y) {
      y[0] = x[0] + x[1];
      Number c;
      {
        Number b = x[0] * x[0];
        c = b + x[0];
      }
      Number d = 2.72;
    
#if REVERSE_TAPE
      Number::getGlobalTape().registerInput(d);
#endif
    
      y[0] = d * c;
    }
};
