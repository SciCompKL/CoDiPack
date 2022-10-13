/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <omp.h>
#include <vector>

#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../misc/macros.hpp"
#include "../../tapes/misc/vectorAccessInterface.hpp"
#include "../../tapes/interfaces/fullTapeInterface.hpp"
#include "../../traits/tapeTraits.hpp"
#include "../data/externalFunctionUserData.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Helper class for the implementation of an external function in CoDiPack.
   *
   * See ExternalFunctionHelper for the basic concepts.
   *
   * This is a modified implementation for use in an OpenMP parallel environment. Is is intended to be used in a shared
   * fashion where multiple threads access one external function helper concurrently. For details, see the documentation
   * of the member functions.
   *
   * @tparam T_Type  The CoDiPack type that is used outside of the external function.
   */
  template<typename T_Type>
  struct OpenMPExternalFunctionHelper {
    public:

      /// See OpenMPExternalFunctionHelper.
      using Type = CODI_DD(T_Type, CODI_T(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

      using Real = typename Type::Real;              ///< See LhsExpressionInterface.
      using Identifier = typename Type::Identifier;  ///< See LhsExpressionInterface.

      /// See LhsExpressionInterface.
      using Tape = CODI_DD(typename Type::Tape, CODI_T(FullTapeInterface<double, double, int, CODI_ANY>));

      /// \copydoc codi::ExternalFunctionHelper::ReverseFunc <br><br>
      /// Note that the parameters are shared by all threads executing the ReverseFunc.
      using ReverseFunc = void (*)(Real const* x, Real* x_b, size_t m, Real const* y, Real const* y_b, size_t n,
                                   ExternalFunctionUserData* d);

      /// \copydoc codi::ExternalFunctionHelper::ForwardFunc <br><br>
      /// Note that the parameters are shared by all threads executing the ForwardFunc.
      using ForwardFunc = void (*)(Real const* x, Real const* x_d, size_t m, Real* y, Real* y_d, size_t n,
                                   ExternalFunctionUserData* d);

      /// \copydoc codi::ExternalFunctionHelper::PrimalFunc <br><br>
      /// Note that the parameters are shared by all threads executing the PrimalFunc.
      using PrimalFunc = void (*)(Real const* x, size_t m, Real* y, size_t n, ExternalFunctionUserData* d);

    private:
      struct EvalData {
        public:

          std::vector<Identifier> inputIndices;
          std::vector<Identifier> outputIndices;

          std::vector<Real> inputValues;
          std::vector<Real> outputValues;
          std::vector<Real> oldPrimals;

          Real* x_d;
          Real* y_d;
          Real* x_b;
          Real* y_b;

          ReverseFunc reverseFunc;
          ForwardFunc forwardFunc;
          PrimalFunc primalFunc;

          ExternalFunctionUserData userData;

          EvalData()
              : inputIndices(0),
                outputIndices(0),
                inputValues(0),
                outputValues(0),
                oldPrimals(0),
                x_d(nullptr),
                y_d(nullptr),
                x_b(nullptr),
                y_b(nullptr),
                reverseFunc(nullptr),
                forwardFunc(nullptr),
                primalFunc(nullptr) {}

          static void delFunc(Tape* t, void* d) {
            CODI_UNUSED(t);

            EvalData* data = (EvalData*)d;

            delete data;
          }

          /// Must be called by all threads of the current team.
          static void evalForwFuncStatic(Tape* t, void* d, VectorAccessInterface<Real, Identifier>* ra) {
            EvalData* data = (EvalData*)d;

            if (nullptr != data->forwardFunc) {
              data->evalForwFunc(t, ra);
            } else {
              CODI_EXCEPTION(
                  "Calling forward evaluation in external function helper without a forward function pointer.");
            }
          }

          /// Compared to the standard implementation, synchronization is introduced between the evaluation steps.
          /// Must be called by all threads of the current team.
          void evalForwFunc(Tape* t, VectorAccessInterface<Real, Identifier>* ra) {
            CODI_UNUSED(t);

            #pragma omp master
            {
              x_d = new Real[inputIndices.size()];
              y_d = new Real[outputIndices.size()];

              if (TapeTraits::IsPrimalValueTape<Tape>::value) {
                for (size_t i = 0; i < inputIndices.size(); ++i) {
                  inputValues[i] = ra->getPrimal(inputIndices[i]);
                }
              }
            }

            #pragma omp barrier

            for (size_t dim = 0; dim < ra->getVectorSize(); ++dim) {
              #pragma omp master
              {
                for (size_t i = 0; i < inputIndices.size(); ++i) {
                  x_d[i] = ra->getAdjoint(inputIndices[i], dim);
                }
              }

              #pragma omp barrier

              forwardFunc(inputValues.data(), x_d, inputIndices.size(), outputValues.data(), y_d, outputIndices.size(),
                          &userData);

              #pragma omp barrier

              #pragma omp master
              {
                for (size_t i = 0; i < outputIndices.size(); ++i) {
                  ra->resetAdjoint(outputIndices[i], dim);
                  ra->updateAdjoint(outputIndices[i], dim, y_d[i]);
                }
              }

              #pragma omp barrier
            }

            #pragma omp master
            {
              if (TapeTraits::IsPrimalValueTape<Tape>::value) {
                for (size_t i = 0; i < outputIndices.size(); ++i) {
                  ra->setPrimal(outputIndices[i], outputValues[i]);
                }
              }

              delete[] x_d;
              delete[] y_d;
            }

            #pragma omp barrier
          }

          /// Must be called by all threads of the current team.
          static void evalPrimFuncStatic(Tape* t, void* d, VectorAccessInterface<Real, Identifier>* ra) {
            EvalData* data = (EvalData*)d;

            if (nullptr != data->primalFunc) {
              data->evalPrimFunc(t, ra);
            } else {
              CODI_EXCEPTION(
                  "Calling primal evaluation in external function helper without a primal function pointer.");
            }
          }

          /// Compared to the standard implementation, synchronization is introduced between the evaluation steps.
          /// Must be called by all threads of the current team.
          void evalPrimFunc(Tape* t, VectorAccessInterface<Real, Identifier>* ra) {
            CODI_UNUSED(t);

            #pragma omp master
            {
              if (TapeTraits::IsPrimalValueTape<Tape>::value) {
                for (size_t i = 0; i < inputIndices.size(); ++i) {
                  inputValues[i] = ra->getPrimal(inputIndices[i]);
                }
              }
            }

            #pragma omp barrier

            primalFunc(inputValues.data(), inputIndices.size(), outputValues.data(), outputIndices.size(), &userData);

            #pragma omp barrier

            #pragma omp master
            {
              if (TapeTraits::IsPrimalValueTape<Tape>::value) {
                for (size_t i = 0; i < outputIndices.size(); ++i) {
                  ra->setPrimal(outputIndices[i], outputValues[i]);
                }
              }
            }

            #pragma omp barrier
          }

          /// Must be called by all threads of the current team.
          static void evalRevFuncStatic(Tape* t, void* d, VectorAccessInterface<Real, Identifier>* ra) {
            EvalData* data = (EvalData*)d;

            if (nullptr != data->reverseFunc) {
              data->evalRevFunc(t, ra);
            } else {
              CODI_EXCEPTION(
                  "Calling reverse evaluation in external function helper without a reverse function pointer.");
            }
          }

          /// Compared to the standard implementation, synchronization is introduced between the evaluation steps.
          /// Must be called by all threads of the current team.
          void evalRevFunc(Tape* t, VectorAccessInterface<Real, Identifier>* ra) {
            CODI_UNUSED(t);

            #pragma omp master
            {
              x_b = new Real[inputIndices.size()];
              y_b = new Real[outputIndices.size()];
            }

            #pragma omp barrier

            for (size_t dim = 0; dim < ra->getVectorSize(); ++dim) {
              #pragma omp master
              {
                for (size_t i = 0; i < outputIndices.size(); ++i) {
                  y_b[i] = ra->getAdjoint(outputIndices[i], dim);
                  ra->resetAdjoint(outputIndices[i], dim);
                }
              }

              #pragma omp barrier

              reverseFunc(inputValues.data(), x_b, inputIndices.size(), outputValues.data(), y_b, outputIndices.size(),
                          &userData);

              #pragma omp barrier

              #pragma omp master
              {
                for (size_t i = 0; i < inputIndices.size(); ++i) {
                  ra->updateAdjoint(inputIndices[i], dim, x_b[i]);
                }
              }

              #pragma omp barrier
            }

            #pragma omp master
            {
              if (Tape::RequiresPrimalRestore) {
                for (size_t i = 0; i < outputIndices.size(); ++i) {
                  ra->setPrimal(outputIndices[i], oldPrimals[i]);
                }
              }

              delete[] x_b;
              delete[] y_b;
            }

            #pragma omp barrier
          }
      };

    protected:

      std::vector<Type*> outputValues;  ///< \copydoc ExternalFunctionHelper::outputValues

      bool storeInputPrimals;     ///< \copydoc ExternalFunctionHelper::storeInputPrimals
      bool storeOutputPrimals;    ///< \copydoc ExternalFunctionHelper::storeOutputPrimals
      bool primalFuncUsesADType;  ///< \copydoc ExternalFunctionHelper::primalFuncUsesADType

      EvalData* data;  ///< \copydoc ExternalFunctionHelper::data

      Real* y;  ///< Shared array of output variables.

    public:

      /// Constructor
      OpenMPExternalFunctionHelper(bool primalFuncUsesADType = false)
          : outputValues(),
            storeInputPrimals(true),
            storeOutputPrimals(true),
            primalFuncUsesADType(primalFuncUsesADType),
            data(nullptr),
            y(nullptr) {
        data = new EvalData();
      }

      /// Destructor
      ~OpenMPExternalFunctionHelper() {
        delete data;
      }

      /// \copydoc ExternalFunctionHelper::disableInputPrimalStore <br><br>
      /// Not thread-safe.
      void disableInputPrimalStore() {
        storeInputPrimals = false;
      }

      /// \copydoc ExternalFunctionHelper::disableOutputPrimalStore <br><br>
      /// Not thread-safe.
      void disableOutputPrimalStore() {
        storeOutputPrimals = false;
      }

      /// \copydoc ExternalFunctionHelper::addInput <br><br>
      /// Not thread-safe.
      void addInput(Type const& input) {
        if (Type::getTape().isActive()) {
          data->inputIndices.push_back(input.getIdentifier());
        }

        // Ignore the setting at this place and the active check,
        // we might need the values for the evaluation.
        if (!primalFuncUsesADType || storeInputPrimals) {
          data->inputValues.push_back(input.getValue());
        }
      }

    private:

      void addOutputToData(Type& output) {
        Real oldPrimal = Type::getTape().registerExternalFunctionOutput(output);

        data->outputIndices.push_back(output.getIdentifier());
        if (storeOutputPrimals) {
          data->outputValues.push_back(output.getValue());
        }
        if (Tape::RequiresPrimalRestore) {
          data->oldPrimals.push_back(oldPrimal);
        }
      }

    public:

      /// \copydoc ExternalFunctionHelper::addOutput <br><br>
      /// Not thread-safe.
      void addOutput(Type& output) {
        if (Type::getTape().isActive()) {
          outputValues.push_back(&output);
        }
      }

      /// \copydoc ExternalFunctionHelper::addUserData <br><br>
      /// Not thread-safe.
      template<typename Data>
      void addUserData(Data const& data) {
        this->data->userData.addData(data);
      }

      /// \copydoc ExternalFunctionHelper::getExternalFunctionUserData <br><br>
      ExternalFunctionUserData& getExternalFunctionUserData() {
        return this->data->userData;
      }

      /// \copydoc ExternalFunctionHelper::callPrimalFuncWithADType <br><br>
      /// Must be called by all threads of the current team.
      template<typename FuncObj, typename... Args>
      void callPrimalFuncWithADType(FuncObj& func, Args&&... args) {
        bool isTapeActive = Type::getTape().isActive();

        if (isTapeActive) {
          Type::getTape().setPassive();
        }

        func(std::forward<Args>(args)...);

        #pragma omp barrier

        if (isTapeActive) {
          Type::getTape().setActive();

          #pragma omp master
          {
            for (size_t i = 0; i < outputValues.size(); ++i) {
              addOutputToData(*outputValues[i]);
            }
          }
        }

        #pragma omp barrier
      }

      /// \copydoc ExternalFunctionHelper::callPrimalFunc <br><br>
      /// Must be called by all threads of the current team.
      void callPrimalFunc(PrimalFunc func) {
        if (!primalFuncUsesADType) {
          // Store the primal function in the external function data so that it can be used for primal evaluations of
          // the tape.

          #pragma omp master
          {
            data->primalFunc = func;

            y = new Real[outputValues.size()];
          }

          #pragma omp barrier

          func(data->inputValues.data(), data->inputValues.size(), y, outputValues.size(), &data->userData);

          #pragma omp barrier

          #pragma omp master
          {
            // Set the primal values on the output values and add them to the data for the reverse evaluation.
            for (size_t i = 0; i < outputValues.size(); ++i) {
              outputValues[i]->setValue(y[i]);

              if (Type::getTape().isActive()) {
                addOutputToData(*outputValues[i]);
              }
            }

            delete[] y;
          }

          #pragma omp barrier

        } else {
          CODI_EXCEPTION(
              "callPrimalFunc() not available if external function helper is initialized with passive function mode "
              "enabled. Use callPrimalFuncWithADType() instead.");
        }
      }

      /// \copydoc ExternalFunctionHelper::addToTape <br><br>
      /// Must be called by all threads of the current team.
      void addToTape(ReverseFunc reverseFunc, ForwardFunc forwardFunc = nullptr, PrimalFunc primalFunc = nullptr) {
        if (Type::getTape().isActive()) {

          #pragma omp master
          {
            data->reverseFunc = reverseFunc;
            data->forwardFunc = forwardFunc;

            if (nullptr != primalFunc) {
              // Only overwrite the primal function if the user provides one, otherwise it is set in the callPrimalFunc
              // method.
              data->primalFunc = primalFunc;
            }

            // Clear the primal values if they are not required.
            if (!storeInputPrimals) {
              data->inputValues.clear();
            }
          }

          #pragma omp barrier

          // make sure that the delete handle is pushed onto only one tape in a parallel context
          if (omp_in_parallel() && omp_get_thread_num() == 0) {
            Type::getTape().pushExternalFunction(
                ExternalFunction<Tape>::create(EvalData::evalRevFuncStatic, data, EvalData::delFunc,
                                               EvalData::evalForwFuncStatic, EvalData::evalPrimFuncStatic));
          } else {
            Type::getTape().pushExternalFunction(
                ExternalFunction<Tape>::create(EvalData::evalRevFuncStatic, data, nullptr,
                                               EvalData::evalForwFuncStatic, EvalData::evalPrimFuncStatic));
          }

          #pragma omp barrier

          #pragma omp master
          {
            data = nullptr;
          }
        } else {
          #pragma omp master
          {
            // Clear the assembled data.
            delete data;
          }
        }

        #pragma omp master
        {
          // Create a new data object for the next call.
          data = new EvalData();
          outputValues.clear();
        }

        #pragma omp barrier
      }
  };
}
