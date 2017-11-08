/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2017 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * CoDiPack is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU
 * General Public License along with CoDiPack.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#pragma once

#include "dataStore.hpp"
#include "../adjointInterface.hpp"

#include <vector>

namespace codi {

  template<typename CoDiType>
  class ExternalFunctionData {
    public:
      typedef typename CoDiType::Real Real;
      typedef typename CoDiType::PassiveReal PassiveReal;
      typedef typename CoDiType::GradientData GradientData;
      typedef typename CoDiType::GradientValue GradientValue;

      typedef typename CoDiType::TapeType Tape;

      typedef void (*ReverseFunc)(const Real* x, Real* x_b, size_t m, const Real* y, const Real* y_b, size_t n, DataStore* d);

      std::vector<GradientData> inputIndices;
      std::vector<GradientData> outputIndices;

      std::vector<Real> inputValues;
      std::vector<Real> outputValues;
      std::vector<Real> oldPrimals;

      ReverseFunc revFunc;

      DataStore userData;

      static void delFunc(void* t, void* d) {
        ExternalFunctionData<CoDiType>* data = (ExternalFunctionData<CoDiType>*)d;

        delete data;
      }

      static void evalRevFuncStatic(void* t, void* d, void* ra) {
        ExternalFunctionData<CoDiType>* data = (ExternalFunctionData<CoDiType>*)d;

        data->evalRevFunc((Tape*)t, (AdjointInterface<Real>*)ra);
      }

      void evalRevFunc(Tape* t, AdjointInterface<Real>* ra) {
        Real* x_b = new Real[inputIndices.size()];
        Real* y_b = new Real[outputIndices.size()];

        for(size_t dim = 0; dim < ra->getVectorSize(); ++dim) {

          for(size_t i = 0; i < outputIndices.size(); ++i) {
            y_b[i] = ra->getAdjoint(outputIndices[i], dim);
          }

          revFunc(inputValues.data(), x_b, inputIndices.size(), outputValues.data(), y_b, outputIndices.size(), &userData);

          for(size_t i = 0; i < inputIndices.size(); ++i) {
            ra->updateAdjoint(inputIndices[i], dim, x_b[i]);
          }
        }

        if(Tape::RequiresPrimalReset) {
          for(size_t i = 0; i < outputIndices.size(); ++i) {
            ra->resetPrimal(outputIndices[i], oldPrimals[i]);
          }
        }
      }
  };

  template<typename CoDiType>
  class ExternalFunctionHelper {
    public:
      typedef typename CoDiType::Real Real;
      typedef typename CoDiType::GradientData GradientData;
      typedef typename CoDiType::GradientValue GradientValue;

      typedef typename CoDiType::TapeType Tape;

      typedef void (*PrimalFunc)(const Real* x, size_t m, Real* y, size_t n, DataStore* d);
      typedef typename ExternalFunctionData<CoDiType>::ReverseFunc ReverseFunc;

      std::vector<CoDiType*> outputValues;

      bool storeInputPrimals;
      bool storeOutputPrimals;
      bool isPassiveExtFunc;
      bool isTapeActive;

      ExternalFunctionData<CoDiType> *data;

      ExternalFunctionHelper() :
        outputValues(),
        storeInputPrimals(true),
        storeOutputPrimals(true),
        isPassiveExtFunc(false),
        isTapeActive(CoDiType::getGlobalTape().isActive()),
        data(nullptr) {
        data = new ExternalFunctionData<CoDiType>();
      }

      ~ExternalFunctionHelper() {
        if(!isTapeActive) {
          delete data;
        }
      }

      void disableInputPrimalStore() {
        storeInputPrimals = false;
      }

      void disableOutputPrimalStore() {
        storeOutputPrimals = false;
      }

      void addInput(const CoDiType& input) {
        if(isTapeActive) {
          data->inputIndices.push_back(input.getGradientData());
        }

        // ignore the setting at this place and the active check
        // We might need the values for the evaluation.
        data->inputValues.push_back(input.getValue());
      }

    private:
      void addOutputToData(CoDiType& output) {
        Real oldPrimal = CoDiType::getGlobalTape().registerExtFunctionOutput(output);

        data->outputIndices.push_back(output.getGradientData());
        if(storeOutputPrimals) {
          data->outputValues.push_back(output.getValue());
        }
        if(Tape::RequiresPrimalReset) {
          data->oldPrimals.push_back(oldPrimal);
        }
      }

    public:
      void addOutput(CoDiType& output) {
        if(isTapeActive) {
          if(isPassiveExtFunc) {
            addOutputToData(output);
          } else {
            outputValues.push_back(&output);
          }
        }
      }

      template<typename Data>
      void addUserData(const Data& data) {
        if(isTapeActive) {
          this->data->userData.addData(data);
        }
      }

      template<typename FuncObj, typename ... Args>
      void callPassiveFunc(FuncObj& func, Args&& ... args) {
        isPassiveExtFunc = true;

        if(isTapeActive) {
          CoDiType::getGlobalTape().setPassive();
        }

        func(std::forward<Args>(args)...);

        if(isTapeActive) {
          CoDiType::getGlobalTape().setActive();
        }
      }

      void callPrimalFunc(PrimalFunc func) {
        Real* y = new Real[outputValues.size()];

        func(data->inputValues.data(), data->inputValues.size(), y, outputValues.size(), &data->userData);

        // ok now set the primal values on the output values and add them to the data for the reverse evaluation
        for(size_t i = 0; i < outputValues.size(); ++i) {
          outputValues[i]->setValue(y[i]);

          addOutputToData(*outputValues[i]);
        }

        delete [] y;
      }

      void addToTape(ReverseFunc func) {
        if(isTapeActive) {

          data->revFunc = func;

          // clear now the primal values if they are not required
          if(!storeInputPrimals) {
            data->inputValues.clear();
          }

          CoDiType::getGlobalTape().pushExternalFunctionHandle(ExternalFunctionData<CoDiType>::evalRevFuncStatic, data, ExternalFunctionData<CoDiType>::delFunc);
        }
      }
  };
}
