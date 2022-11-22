#pragma once

#include "reverseDriver.hpp"

template<typename Number>
struct ReverseDriverPrimal : public ReverseDriver<Number> {

    using Tape = typename Number::Tape;

    void evaluate(Tape& tape) {
      tape.evaluate();

#if (TEST_NAME != TestPreacc)
      tape.evaluatePrimal();
#endif
    }
};
