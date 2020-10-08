#pragma once

#include <vector>

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../tapes/interfaces/fullTapeInterface.hpp"
#include "../../tapes/aux/vectorAccessInterface.hpp"
#include "../../traits/tapeTraits.hpp"
#include "../data/externalFunctionUserData.hpp"


/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Type>
  struct ExternalFunctionHelper {
    public:

      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

      using Real = typename Type::Real;
      using Identifier = typename Type::Identifier;
      using Tape = CODI_DECLARE_DEFAULT(typename Type::Tape, CODI_TEMPLATE(FullTapeInterface<double, double, int, CODI_ANY>));

      using ReverseFunc = void (*)(Real const* x, Real* x_b, size_t m, Real const* y, Real const* y_b, size_t n, ExternalFunctionUserData* d);
      using ForwardFunc = void (*)(Real const* x, Real const* x_d, size_t m, Real* y, Real* y_d, size_t n, ExternalFunctionUserData* d);
      using PrimalFunc = void (*)(Real const* x, size_t m, Real* y, size_t n, ExternalFunctionUserData* d);

    private:
      struct EvalData {
        public:

          std::vector<Identifier> inputIndices;
          std::vector<Identifier> outputIndices;

          std::vector<Real> inputValues;
          std::vector<Real> outputValues;
          std::vector<Real> oldPrimals;

          ReverseFunc reverseFunc;
          ForwardFunc forwardFunc;
          PrimalFunc primalFunc;

          ExternalFunctionUserData userData;

          EvalData() :
            inputIndices(0),
            outputIndices(0),
            inputValues(0),
            outputValues(0),
            oldPrimals(0),
            reverseFunc(nullptr),
            forwardFunc(nullptr),
            primalFunc(nullptr) {}

          static void delFunc(Tape* t, void* d) {
            CODI_UNUSED(t);

            EvalData* data = (EvalData*)d;

            delete data;
          }

          static void evalForwFuncStatic(Tape* t, void* d, VectorAccessInterface<Real, Identifier>* ra) {
            EvalData* data = (EvalData*)d;

            if(nullptr != data->forwardFunc) {
              data->evalForwFunc(t, ra);
            } else {
              CODI_EXCEPTION("Calling forward evaluation in external function helper without a forward function pointer.");
            }
          }

          void evalForwFunc(Tape* t, VectorAccessInterface<Real, Identifier>* ra) {
            CODI_UNUSED(t);

            Real* x_d = new Real[inputIndices.size()];
            Real* y_d = new Real[outputIndices.size()];

            if(isPrimalValueTape<Tape>::value) {
              for(size_t i = 0; i < inputIndices.size(); ++i) {
                inputValues[i] = ra->getPrimal(inputIndices[i]);
              }
            }

            for(size_t dim = 0; dim < ra->getVectorSize(); ++dim) {

              for(size_t i = 0; i < inputIndices.size(); ++i) {
                x_d[i] = ra->getAdjoint(inputIndices[i], dim);
              }

              forwardFunc(inputValues.data(), x_d, inputIndices.size(), outputValues.data(), y_d, outputIndices.size(), &userData);

              for(size_t i = 0; i < outputIndices.size(); ++i) {
                ra->resetAdjoint(outputIndices[i], dim);
                ra->updateAdjoint(outputIndices[i], dim, y_d[i]);
              }
            }

            if(isPrimalValueTape<Tape>::value) {
              for(size_t i = 0; i < outputIndices.size(); ++i) {
                ra->setPrimal(outputIndices[i], outputValues[i]);
              }
            }

            delete [] x_d;
            delete [] y_d;
          }

          static void evalPrimFuncStatic(Tape* t, void* d, VectorAccessInterface<Real, Identifier>* ra) {
            EvalData* data = (EvalData*)d;

            if(nullptr != data->primalFunc) {
              data->evalPrimFunc(t, ra);
            } else {
              CODI_EXCEPTION("Calling primal evaluation in external function helper without a primal function pointer.");
            }
          }

          void evalPrimFunc(Tape* t, VectorAccessInterface<Real, Identifier>* ra) {
            CODI_UNUSED(t);

            if(isPrimalValueTape<Tape>::value) {
              for(size_t i = 0; i < inputIndices.size(); ++i) {
                inputValues[i] = ra->getPrimal(inputIndices[i]);
              }
            }

            primalFunc(inputValues.data(), inputIndices.size(), outputValues.data(), outputIndices.size(), &userData);

            if(isPrimalValueTape<Tape>::value) {
              for(size_t i = 0; i < outputIndices.size(); ++i) {
                ra->setPrimal(outputIndices[i], outputValues[i]);
              }
            }
          }

          static void evalRevFuncStatic(Tape* t, void* d, VectorAccessInterface<Real, Identifier>* ra) {
            EvalData* data = (EvalData*)d;

            if(nullptr != data->reverseFunc) {
              data->evalRevFunc(t, ra);
            } else {
              CODI_EXCEPTION("Calling reverse evaluation in external function helper without a reverse function pointer.");
            }
          }

          void evalRevFunc(Tape* t, VectorAccessInterface<Real, Identifier>* ra) {
            CODI_UNUSED(t);

            Real* x_b = new Real[inputIndices.size()];
            Real* y_b = new Real[outputIndices.size()];

            for(size_t dim = 0; dim < ra->getVectorSize(); ++dim) {

              for(size_t i = 0; i < outputIndices.size(); ++i) {
                y_b[i] = ra->getAdjoint(outputIndices[i], dim);
                ra->resetAdjoint(outputIndices[i], dim);
              }

              reverseFunc(inputValues.data(), x_b, inputIndices.size(), outputValues.data(), y_b, outputIndices.size(), &userData);

              for(size_t i = 0; i < inputIndices.size(); ++i) {
                ra->updateAdjoint(inputIndices[i], dim, x_b[i]);
              }
            }

            if(Tape::RequiresPrimalRestore) {
              for(size_t i = 0; i < outputIndices.size(); ++i) {
                ra->setPrimal(outputIndices[i], oldPrimals[i]);
              }
            }

            delete [] x_b;
            delete [] y_b;
          }
      };

    public:

      std::vector<Type*> outputValues;

      bool storeInputPrimals;
      bool storeOutputPrimals;
      bool primalFuncUsesADType;

      EvalData* data;

      ExternalFunctionHelper(bool primalFuncUsesADType = false) :
        outputValues(),
        storeInputPrimals(true),
        storeOutputPrimals(true),
        primalFuncUsesADType(primalFuncUsesADType),
        data(nullptr) {
        data = new EvalData();
      }

      ~ExternalFunctionHelper() {
        delete data;
      }

      void disableInputPrimalStore() {
        storeInputPrimals = false;
      }

      void disableOutputPrimalStore() {
        storeOutputPrimals = false;
      }

      void addInput(Type const& input) {
        if(Type::getGlobalTape().isActive()) {
          data->inputIndices.push_back(input.getIdentifier());
        }

        // ignore the setting at this place and the active check
        // We might need the values for the evaluation.
        if (!primalFuncUsesADType || storeInputPrimals){
          data->inputValues.push_back(input.getValue());
        }
      }

    private:

      void addOutputToData(Type& output) {
        Real oldPrimal = Type::getGlobalTape().registerExternalFunctionOutput(output);

        data->outputIndices.push_back(output.getIdentifier());
        if(storeOutputPrimals) {
          data->outputValues.push_back(output.getValue());
        }
        if(Tape::RequiresPrimalRestore) {
          data->oldPrimals.push_back(oldPrimal);
        }
      }

    public:

      void addOutput(Type& output) {
        if(Type::getGlobalTape().isActive()) {
          outputValues.push_back(&output);
        }
      }

      template<typename Data>
      void addUserData(Data const& data) {
        this->data->userData.addData(data);
      }

      ExternalFunctionUserData& getExternalFunctionUserData() {
        return this->data->userData;
      }

      template<typename FuncObj, typename ... Args>
      void callPrimalFuncWithADType(FuncObj& func, Args&& ... args) {

        bool isTapeActive = Type::getGlobalTape().isActive();

        if(isTapeActive) {
          Type::getGlobalTape().setPassive();
        }

        func(std::forward<Args>(args)...);

        if(isTapeActive) {
          Type::getGlobalTape().setActive();

          for(size_t i = 0; i < outputValues.size(); ++i) {
            addOutputToData(*outputValues[i]);
          }
        }
      }

      void callPrimalFunc(PrimalFunc func) {
        if (!primalFuncUsesADType){
          // First set the function here in the external function data so that it can be used for primal evaluation of the tape.
          data->primalFunc = func;

          Real* y = new Real[outputValues.size()];

          func(data->inputValues.data(), data->inputValues.size(), y, outputValues.size(), &data->userData);

          // ok now set the primal values on the output values and add them to the data for the reverse evaluation
          for(size_t i = 0; i < outputValues.size(); ++i) {
            outputValues[i]->setValue(y[i]);

            if(Type::getGlobalTape().isActive()) {
              addOutputToData(*outputValues[i]);
            }
          }

          delete [] y;

        } else {
          CODI_EXCEPTION("callPrimalFunc() not available if external function helper is initialized with passive function mode enabled. Use callPrimalFuncWithADType() instead.");
        }
      }

      void addToTape(ReverseFunc reverseFunc, ForwardFunc forwardFunc = nullptr, PrimalFunc primalFunc = nullptr) {
        if(Type::getGlobalTape().isActive()) {

          data->reverseFunc = reverseFunc;
          data->forwardFunc = forwardFunc;

          if(nullptr != primalFunc) {
            // Only overwrite if the user provides one. Otherwise it is set in the callPrimalFunc method
            data->primalFunc = primalFunc;
          }

          // clear now the primal values if they are not required
          if(!storeInputPrimals) {
            data->inputValues.clear();
          }

          Type::getGlobalTape().pushExternalFunction(ExternalFunction<Tape>::create(
              EvalData::evalRevFuncStatic,
              data,
              EvalData::delFunc,
              EvalData::evalForwFuncStatic,
              EvalData::evalPrimFuncStatic));

          data = nullptr;
        } else {
          // clear the assembled data
          delete data;
        }

        // create a new data object for the next call
        data = new EvalData();
        outputValues.clear();
      }
  };
}
