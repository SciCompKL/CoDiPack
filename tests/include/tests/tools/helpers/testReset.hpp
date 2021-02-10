#include "../../../testInterface.hpp"

struct TestReset : public TestInterface {
  public:
    NAME("Reset")
    IN(2)
    OUT(1)
    POINTS(1) = {{1.0, 0.5}};

    template<typename Number>
    static void func(Number* x, Number* y) {

      #if REVERSE_TAPE
        typename Number::Tape& tape = Number::getGlobalTape();
      #endif

      Number a = x[0] * x[1];
      Number b = x[0] / sin(x[1]);
      Number c = b * a;

      #if REVERSE_TAPE
        typename Number::Tape::Position pos = tape.getPosition();
      #endif

      b = a * x[0];

      #if REVERSE_TAPE
        tape.resetTo(pos);
      #endif

      y[0] = c * a;
    }
};
