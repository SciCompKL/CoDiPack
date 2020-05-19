#pragma once

#include <codi.hpp>
#include <codi/tools/data/jacobian.hpp>

#include "../driver1stOrderBase.hpp"

#include DRIVER_TESTS_INC

struct CoDiReverse1stOrderBase : public Driver1stOrderBase<CODI_TYPE> {
  public:

    using Number = CODI_TYPE;

    using Base = Driver1stOrderBase<Number>;

    using Gradient = Number::Gradient;

    CoDiReverse1stOrderBase() : Base(CODI_TO_STRING(CODI_TYPE_NAME)) {}

    virtual Gradient& accessGradient(Number& value) = 0;
    virtual void cleanup() = 0;
    virtual void evaluate() = 0;
    virtual void prepare() = 0;

    void createAllTests(TestVector<Number>& tests) {
      createTests<Number, DRIVER_TESTS>(tests);
    }

    void evaluateJacobian(TestInfo<Number>& info, Number* x, size_t inputs, Number* y, size_t outputs, codi::Jacobian<double>& jac) {

      //using GT = codi::GradientValueTraits<Gradient>;
      constexpr size_t gradDim = 1; //GT::getVectorSize();

      Number::Tape& tape = Number::getGlobalTape();

      // Set sizes for Jacobian tapes
      if(tape.hasOption(codi::ConfigurationOption::JacobianSize)) {
        tape.setOption(codi::ConfigurationOption::JacobianSize, 1000);
      }
      if(tape.hasOption(codi::ConfigurationOption::StatementSize)) {
        tape.setOption(codi::ConfigurationOption::StatementSize, 1000);
      }
      if(tape.hasOption(codi::ConfigurationOption::ExternalFunctionsSize)) {
        tape.setOption(codi::ConfigurationOption::ExternalFunctionsSize, 1000);
      }

      size_t runs = outputs / gradDim;
      if(outputs % gradDim != 0) {
        runs += 1;
      }
      for(size_t curOut = 0; curOut < runs; ++curOut) {
        size_t curSize = gradDim;
        if((curOut + 1) * gradDim  > (size_t)outputs) {
          curSize = outputs % gradDim;
        }

        prepare();

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
            accessGradient(y[curOut * gradDim + curDim]) = 1.0;
          }
        }

        evaluate();

        for(size_t curDim = 0; curDim < curSize; ++curDim) {
          for(size_t curIn = 0; curIn < inputs; ++curIn) {
            //jac(curOut * gradDim + curDim, curIn) = GT::at(getGradient(x[curIn]), curDim);
#ifdef SECOND_ORDER
            jac(curOut * gradDim + curDim, curIn) = accessGradient(x[curIn]).value();
#else
            jac(curOut * gradDim + curDim, curIn) = accessGradient(x[curIn]);
#endif
          }
        }

        tape.reset();

        cleanup();
      }
    }
};
