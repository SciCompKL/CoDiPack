#pragma once

#include <codi.hpp>
#include <codi/tools/data/jacobian.hpp>

#include "../driver1stOrderBase.hpp"

#include DRIVER_TESTS_INC

struct CoDiForward1stOrder : public Driver1stOrderBase<CODI_TYPE> {
  public:

    using Number = CODI_TYPE;

    using Base = Driver1stOrderBase<Number>;

    CoDiForward1stOrder() : Base(CODI_TO_STRING(CODI_TYPE_NAME)) {}

    void createAllTests(TestVector<Number>& tests) {
      createTests<Number, DRIVER_TESTS>(tests);
    }

    void evaluateJacobian(TestInfo<Number>& info, Number* x, size_t inputs, Number* y, size_t outputs, codi::Jacobian<double>& jac) {

      using Gradient = typename Number::Gradient;
      size_t constexpr gradDim = codi::GradientTraits::dim<Gradient>();

      size_t runs = inputs / gradDim;
      if(inputs % gradDim != 0) {
        runs += 1;
      }
      for(size_t curIn = 0; curIn < runs; ++curIn) {
        size_t curSize = gradDim;
        if((curIn + 1) * gradDim  > (size_t)inputs) {
          curSize = inputs % gradDim;
        }
        for(size_t curDim = 0; curDim < curSize; ++curDim) {
          codi::GradientTraits::at(x[curIn * gradDim + curDim].gradient(), curDim) = 1.0;
        }

        for(size_t i = 0; i < outputs; ++i) {
          y[i].setGradient(Gradient());
        }

        info.func(x, y);

        for(size_t curDim = 0; curDim < curSize; ++curDim) {
          for(size_t curOut = 0; curOut < outputs; ++curOut) {
#ifdef SECOND_ORDER
            jac(curOut, curIn * gradDim + curDim) = codi::GradientTraits::at(y[curOut].getGradient(), curDim).value();
#else
            jac(curOut, curIn * gradDim + curDim) = codi::GradientTraits::at(y[curOut].getGradient(), curDim);
#endif
          }
        }

        for(size_t curDim = 0; curDim < curSize; ++curDim) {
          x[curIn * gradDim + curDim].setGradient(Gradient());
        }
      }
    }
};
