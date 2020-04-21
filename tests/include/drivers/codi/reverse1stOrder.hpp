#pragma once

#include "reverse1stOrderBase.hpp"

struct CoDiReverse1stOrder : public CoDiReverse1stOrderBase {
  public:

    using Number = CODI_TYPE;

    using Base = CoDiReverse1stOrderBase;

    using Gradient = Number::Gradient;

    Gradient& accessGradient(Number &value) {
      return value.gradient();
    }

    void cleanup() {}

    void evaluate() {
      Number::getGlobalTape().evaluate();
    }

    void prepare() {}
};
