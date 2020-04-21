#pragma once

#include "reverse1stOrderBase.hpp"

struct CoDiReverse1stOrderVectorHelper : public CoDiReverse1stOrderBase {
  public:

    using Number = CODI_TYPE;

    using Base = CoDiReverse1stOrderBase;

    using Gradient = Number::Gradient;


    codi::CustomGradientVectorHelper<Number, Gradient> helper;

    CoDiReverse1stOrderVectorHelper() : Base(), helper() {}

    Gradient& accessGradient(Number &value) {
      return helper.gradient(value.getIdentifier());
    }

    void cleanup() {
      helper.clearAdjoints();
    }

    void evaluate() {
      helper.evaluate();
    }

    void prepare() {}
};
