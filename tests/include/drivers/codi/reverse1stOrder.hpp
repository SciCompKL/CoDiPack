#pragma once

#include <codi.hpp>
#include <codi/tools/data/jacobian.hpp>

#include "../driverBase.hpp"

#include DRIVER_TESTS_INC

struct CoDiReverse1stOrder : public DriverBase<CODI_TYPE> {
  public:

    using Number = CODI_TYPE;

    using Base = DriverBase<Number>;

    CoDiReverse1stOrder() : Base(CODI_TO_STRING(CODI_TYPE_NAME)) {}

    void createAllTests(TestVector<Number>& tests) {
      createTests<Number, DRIVER_TESTS>(tests);
    }

    void evaluateJacobian(TestInfo<Number>& info, Number* x, size_t inputs, Number* y, size_t outputs, codi::Jacobian<double>& jac) {

      //using GT = codi::GradientValueTraits<Gradient>;
      constexpr size_t gradDim = 1; //GT::getVectorSize();

      Number::Tape& tape = Number::getGlobalTape();

      size_t runs = outputs / gradDim;
      if(outputs % gradDim != 0) {
        runs += 1;
      }
      for(size_t curOut = 0; curOut < runs; ++curOut) {
        size_t curSize = gradDim;
        if((curOut + 1) * gradDim  > (size_t)outputs) {
          curSize = outputs % gradDim;
        }

        tape.setActive();

        for(size_t i = 0; i < inputs; ++i) {
          tape.registerInput(x[i]);
        }

        info.func(x, y);

        for(size_t i = 0; i < outputs; ++i) {
          tape.registerOutput(y[i]);
        }

        for(size_t curDim = 0; curDim < curSize; ++curDim) {
          if(tape.isIdentifierActive(y[curOut * gradDim + curDim].getIdentifier())) {
            y[curOut * gradDim + curDim].gradient() = 1.0;
          }
        }

        tape.evaluate();

        for(size_t curDim = 0; curDim < curSize; ++curDim) {
          for(size_t curIn = 0; curIn < inputs; ++curIn) {
            //jac(curOut * gradDim + curDim, curIn) = GT::at(getGradient(x[curIn]), curDim);
            jac(curOut * gradDim + curDim, curIn) = x[curIn].getGradient();
          }
        }

        tape.reset();
      }
    }
};
