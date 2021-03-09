#include <type_traits>
#include <vector>

#include "../../../../include/codi.hpp"
#include "../../testInterface.hpp"

struct TestReferenceActiveType : public TestInterface {
  public:
    NAME("ReferenceActiveType")
    IN(1)
    OUT(1)
    POINTS(1) = {// clang-format off
      {0.5}
    };  // clang-format on

#if !defined(DOUBLE)
    template<typename Number>
    using RefReal = codi::ReferenceActiveType<Number>;
#else
    template<typename Number>
    using RefReal = Number;
#endif

    template<typename Number>
    static void func(Number* x, Number* y) {
      RefReal<Number> xRef = x[0];
      y[0] = 3.0 * xRef * xRef * xRef * xRef + 5.0 * xRef * xRef * xRef - 3.0 * xRef * xRef + 2.0 * xRef - 4.0;
    }
};
