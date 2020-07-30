#pragma once

#include "../../../../include/codi.hpp"

template<typename _Number, typename = void>
struct MultiplyExternalFunction {
  public:
    using Number = CODI_DECLARE_DEFAULT(_Number, codi::ActiveType<CODI_ANY>);

    static Number create(Number const& x1, Number const& x2) {
      return x1 * x2;
    }
};


template<typename _Number>
struct MultiplyExternalFunction<_Number, codi::enableIfReverseTape<typename _Number::Tape>> {
  public:
    using Number = CODI_DECLARE_DEFAULT(_Number, codi::ActiveType<CODI_ANY>);

    using Tape = CODI_DECLARE_DEFAULT(typename Number::Tape, CODI_TEMPLATE(codi::FullTapeInterface<double, double, int, codi::EmptyPosition>));
    using Real = typename Number::Real;
    using Identifier = typename Number::Identifier;

    using VAI = codi::VectorAccessInterface<Real, Identifier>;

    static void extFuncReverse(void* t, void* d, void* i) {
      codi::CODI_UNUSED(t);

      VAI* vai = (VAI*)i;

      codi::ExternalFunctionData *data = static_cast<codi::ExternalFunctionData*>(d);

      Real x1_v, x2_v;
      Identifier x1_i, x2_i, w_i;
      data->getData(x1_v);
      data->getData(x1_i);
      data->getData(x2_v);
      data->getData(x2_i);
      data->getData(w_i);

      size_t dim = vai->getVectorSize();

      for(size_t i = 0; i < dim; ++i) {

        Real w_b = vai->getAdjoint(w_i, i);
        vai->resetAdjoint(w_i, i);

        vai->updateAdjoint(x1_i, i, x2_v*w_b);
        vai->updateAdjoint(x2_i, i, x1_v*w_b);
      }
    }

    static void extFuncPrimal(void* t, void* d, void* i) {
      codi::CODI_UNUSED(t);

      VAI* vai = (VAI*)i;

      codi::ExternalFunctionData *data = static_cast<codi::ExternalFunctionData*>(d);

      Identifier x1_i, x2_i, w_i;
      Real& x1_v = data->getDataRef<Real>();
      data->getData(x1_i);
      Real& x2_v = data->getDataRef<Real>();
      data->getData(x2_i);
      data->getData(w_i);

      x1_v = vai->getPrimal(x1_i); // Data is overwritten here
      x2_v = vai->getPrimal(x2_i); // Data is overwritten here

      Real z = x1_v * x2_v;
      vai->setPrimal( w_i, z);
    }

    static void extFuncForward(void* t, void* d, void* i) {
      codi::CODI_UNUSED(t);

      VAI* vai = (VAI*)i;

      codi::ExternalFunctionData *data = static_cast<codi::ExternalFunctionData*>(d);

      Identifier x1_i, x2_i, w_i;
      Real& x1_v = data->getDataRef<Real>();
      data->getData(x1_i);
      Real& x2_v = data->getDataRef<Real>();
      data->getData(x2_i);
      data->getData(w_i);

      if(vai->hasPrimals()) {
        x1_v = vai->getPrimal(x1_i); // Data is overwritten here
        x2_v = vai->getPrimal(x2_i); // Data is overwritten here
      }

      size_t dim = vai->getVectorSize();

      for(size_t i = 0; i < dim; ++i) {

        Real x1_d = vai->getAdjoint(x1_i, i);
        Real x2_d = vai->getAdjoint(x2_i, i);

        Real w_d = x1_d * x2_v + x1_v * x2_d;
        vai->resetAdjoint(w_i, i);
        vai->updateAdjoint(w_i, i, w_d);
      }

      Real z = x1_v * x2_v;
      vai->setPrimal( w_i, z);
    }

    static void delFunc(void* tape, void* d){
      codi::CODI_UNUSED(tape);

      codi::ExternalFunctionData *data = static_cast<codi::ExternalFunctionData*>(d);
      delete data;
    }

    static Number create(Number const& x1, Number const& x2) {
      Tape& tape = Number::getGlobalTape();
      codi::ExternalFunctionData *data = new codi::ExternalFunctionData;
      Number w;
      tape.setPassive();
      w = x1 * x2;
      tape.setActive();

      tape.registerExternalFunctionOutput(w);
      data->addData(x1.getValue());
      data->addData(x1.getIdentifier());
      data->addData(x2.getValue());
      data->addData(x2.getIdentifier());
      data->addData(w.getIdentifier());
      tape.pushExternalFunction(codi::ExternalFunction(extFuncReverse, extFuncForward, extFuncPrimal, data, delFunc));

      return w;
    }
};
