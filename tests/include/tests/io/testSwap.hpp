#include "../../testInterface.hpp"

struct TestSwap : public TestInterface {
  public:
    NAME("Swap")
    IN(1)
    OUT(1)
    POINTS(1) = {{1.0}};

    template<typename Number>
    static void func(Number* x, Number* y) {
      y[0] = x[0];

#if REVERSE_TAPE
      typename Number::Tape other;
      Number::getTape().swap(other);
      Number::getTape().swap(other);
#endif
    }
};
