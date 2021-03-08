#include "../../testInterface.hpp"

#include <type_traits>
#include <vector>

#include "../../../../include/codi.hpp"

struct TestNumericLimits : public TestInterface {
  public:
    NAME("NumericLimits")
    IN(1)
    OUT(1)
    POINTS(1) = {{1.0}};

    template<typename Number>
    static void func(Number* x, Number* y) {
      y[0] = x[0] + (double)std::numeric_limits<Number>::min_exponent;
    }
};
