/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <vector>

#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../misc/macros.hpp"
#include "../../tapes/interfaces/fullTapeInterface.hpp"
#include "../../tapes/misc/vectorAccessInterface.hpp"
#include "../../traits/tapeTraits.hpp"
#include "../data/externalFunctionUserData.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Helper class for the implementation of an external function in CoDiPack.
   *
   * The class helps the user to handle parts where the CoDiPack types cannot be applied or where a more efficient
   * gradient computation is available.
   *
   * The procedure of pushing an external function with the helper is as follows.
   * 1. All function inputs and outputs are specified.
   * 2. The primal function is called. There are two modes for this, which are explained below.
   * 3. The manual reverse implementation is provided to the helper, which embeds it into the tape and prepares itself
   *    for the next external function push.
   *
   * The first mode of operation assumes that the primal function has an implementation without a CoDiPack type.
   * To use this mode, invoke the primal function via callPrimalFunc. An example is:
   * \snippet examples/externalFunctionHelper.cpp Mode 1: Implemented primal function
   *
   * The second mode of operation assumes that the primal function is evaluated with the CoDiPack type. To use it, the
   * helper's constructor has to be called with the option true and the primal call is performed via
   * callPrimalFuncWithADType. An example is:
   * \snippet examples/externalFunctionHelper.cpp Mode 2: Passive primal function
   *
   * The function implementations must follow the definitions of ExternalFunctionHelper::ReverseFunc,
   * ExternalFunctionHelper::ForwardFunc and ExternalFunctionHelper::PrimalFunc, with an exception for the latter if the
   * second mode is used. The implementations from the examples above are:
   * \snippet examples/externalFunctionHelper.cpp Function implementations
   *
   * The ExternalFunctionHelper works with all tapes. It is also able to handle situations where the tape is currently
   * not recording. All necessary operations are performed in such a case but no external function is recorded.
   * If enablePrimalValueTapeSpecialization is active, primal values are recovered from the tape instead of being stored 
   * in the ExternalFunctionHelper instance, reducing its memory footprint.
   *
   * The storing of primal inputs and outputs can be disabled. Outputs can be discarded if they are recomputed in the
   * derivative computation or if the derivative does not depend on them. Inputs can be discarded if the derivative does
   * not depend on them.
   *
   * @tparam T_Type  The CoDiPack type that is used outside of the external function.
   */
  template<typename T_Type>
  struct ExternalFunctionHelper {
    public:

      /// See ExternalFunctionHelper.
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);

      using Real = typename Type::Real;              ///< See LhsExpressionInterface.
      using Identifier = typename Type::Identifier;  ///< See LhsExpressionInterface.

      using Tape = CODI_DD(typename Type::Tape, CODI_DEFAULT_TAPE);  ///< See LhsExpressionInterface.

      /// Function interface for the reverse AD call of an external function.
      using ReverseFunc = void (*)(Real const* x, Real* x_b, size_t m, Real const* y, Real const* y_b, size_t n,
                                   ExternalFunctionUserData* d);

      /// Function interface for the forward AD call of an external function.
      using ForwardFunc = void (*)(Real const* x, Real const* x_d, size_t m, Real* y, Real* y_d, size_t n,
                                   ExternalFunctionUserData* d);

      /// Function interface for the primal call of an external function.
      using PrimalFunc = void (*)(Real const* x, size_t m, Real* y, size_t n, ExternalFunctionUserData* d);

    private:

      static constexpr bool IsPrimalValueTape = TapeTraits::IsPrimalValueTape<Tape>::value;

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

          bool provideInputValues;
          bool provideOutputValues;
          bool getPrimalsFromPrimalValueVector;

          EvalData()
              : inputIndices(0),
                outputIndices(0),
                inputValues(0),
                outputValues(0),
                oldPrimals(0),
                reverseFunc(nullptr),
                forwardFunc(nullptr),
                primalFunc(nullptr),
                provideInputValues(true),
                provideOutputValues(true),
                getPrimalsFromPrimalValueVector(false) {}

          static void delFunc(Tape* t, void* d) {
            CODI_UNUSED(t);

            EvalData* data = (EvalData*)d;

            delete data;
          }

          static void evalForwFuncStatic(Tape* t, void* d, VectorAccessInterface<Real, Identifier>* ra) {
            EvalData* data = (EvalData*)d;

            if (nullptr != data->forwardFunc) {
              data->evalForwFunc(t, ra);
            } else {
              CODI_EXCEPTION(
                  "Calling forward evaluation in external function helper without a forward function pointer.");
            }
          }

          CODI_INLINE void evalForwFunc(Tape* t, VectorAccessInterface<Real, Identifier>* ra) {
            CODI_UNUSED(t);

            std::vector<Real> x_d(inputIndices.size());
            std::vector<Real> y_d(outputIndices.size());

            initRun(ra);

            for (size_t dim = 0; dim < ra->getVectorSize(); ++dim) {
              for (size_t i = 0; i < inputIndices.size(); ++i) {
                x_d[i] = ra->getAdjoint(inputIndices[i], dim);
              }

              forwardFunc(inputValues.data(), x_d.data(), inputIndices.size(), outputValues.data(), y_d.data(),
                          outputIndices.size(), &userData);

              for (size_t i = 0; i < outputIndices.size(); ++i) {
                ra->resetAdjoint(outputIndices[i], dim);
                ra->updateAdjoint(outputIndices[i], dim, y_d[i]);
              }
            }

            finalizeRun(ra);
          }

          static void evalPrimFuncStatic(Tape* t, void* d, VectorAccessInterface<Real, Identifier>* ra) {
            EvalData* data = (EvalData*)d;

            if (nullptr != data->primalFunc) {
              data->evalPrimFunc(t, ra);
            } else {
              CODI_EXCEPTION(
                  "Calling primal evaluation in external function helper without a primal function pointer.");
            }
          }

          CODI_INLINE void evalPrimFunc(Tape* t, VectorAccessInterface<Real, Identifier>* ra) {
            CODI_UNUSED(t);

            initRun(ra);

            primalFunc(inputValues.data(), inputIndices.size(), outputValues.data(), outputIndices.size(), &userData);

            finalizeRun(ra);
          }

          static void evalRevFuncStatic(Tape* t, void* d, VectorAccessInterface<Real, Identifier>* ra) {
            EvalData* data = (EvalData*)d;

            if (nullptr != data->reverseFunc) {
              data->evalRevFunc(t, ra);
            } else {
              CODI_EXCEPTION(
                  "Calling reverse evaluation in external function helper without a reverse function pointer.");
            }
          }

          CODI_INLINE void evalRevFunc(Tape* t, VectorAccessInterface<Real, Identifier>* ra) {
            CODI_UNUSED(t);

            std::vector<Real> x_b(inputIndices.size());
            std::vector<Real> y_b(outputIndices.size());

            initRun(ra, true);

            for (size_t dim = 0; dim < ra->getVectorSize(); ++dim) {
              for (size_t i = 0; i < outputIndices.size(); ++i) {
                y_b[i] = ra->getAdjoint(outputIndices[i], dim);
                ra->resetAdjoint(outputIndices[i], dim);
              }

              reverseFunc(inputValues.data(), x_b.data(), inputIndices.size(), outputValues.data(), y_b.data(),
                          outputIndices.size(), &userData);

              for (size_t i = 0; i < inputIndices.size(); ++i) {
                ra->updateAdjoint(inputIndices[i], dim, x_b[i]);
              }
            }

            finalizeRun(ra, true);
          }

        private:

          CODI_INLINE void initRun(VectorAccessInterface<Real, Identifier>* ra, bool isReverse = false) {
            if (getPrimalsFromPrimalValueVector && provideOutputValues) {
              outputValues.resize(outputIndices.size());

              if (isReverse) {  // Provide result values for reverse evaluations.
                for (size_t i = 0; i < outputIndices.size(); ++i) {
                  outputValues[i] = ra->getPrimal(outputIndices[i]);
                }
              }
            }

            // Restore the old primals for reverse evaluations, before the inputs are read.
            if (isReverse && Tape::RequiresPrimalRestore) {
              for (size_t i = 0; i < outputIndices.size(); ++i) {
                ra->setPrimal(outputIndices[i], oldPrimals[i]);
              }
            }

            if (getPrimalsFromPrimalValueVector && provideInputValues) {
              inputValues.resize(inputIndices.size());

              for (size_t i = 0; i < inputIndices.size(); ++i) {
                inputValues[i] = ra->getPrimal(inputIndices[i]);
              }
            }
          }

          CODI_INLINE void finalizeRun(VectorAccessInterface<Real, Identifier>* ra, bool isReverse = false) {
            if (getPrimalsFromPrimalValueVector && !isReverse) {
              for (size_t i = 0; i < outputIndices.size(); ++i) {
                if (Tape::RequiresPrimalRestore) {
                  oldPrimals[i] = ra->getPrimal(outputIndices[i]);
                }
                ra->setPrimal(outputIndices[i], outputValues[i]);
              }
            }

            if (getPrimalsFromPrimalValueVector && provideInputValues) {
              inputValues.clear();
              inputValues.shrink_to_fit();
            }
            if (getPrimalsFromPrimalValueVector && provideOutputValues) {
              outputValues.clear();
              outputValues.shrink_to_fit();
            }
          }
      };

    protected:

      std::vector<Type*> outputValues;  ///< References to output values.

      bool storeInputPrimals;              ///< If input primals are stored. Can be disabled by the user.
      bool storeOutputPrimals;             ///< If output primals are stored. Can be disabled by the user.
      bool storeInputOutputForPrimalEval;  ///< If a primal call with a self-implemented function will be done.

      EvalData* data;  ///< External function data.

    public:

      /// Constructor
      ExternalFunctionHelper(bool primalFuncUsesADType = false)
          : outputValues(),
            storeInputPrimals(true),
            storeOutputPrimals(true),
            storeInputOutputForPrimalEval(!primalFuncUsesADType),
            data(nullptr) {
        data = new EvalData();
      }

      /// Destructor
      ~ExternalFunctionHelper() {
        delete data;
      }

      /// Load primal values from the primal value vector of the tape. (Has no effect on Jacobian tapes.)
      void enablePrimalValueTapeSpecialization() {
        if (IsPrimalValueTape) {
          storeInputPrimals = false;
          storeOutputPrimals = false;
          data->getPrimalsFromPrimalValueVector = true;
        }
      }

      /// Do not store primal input values. In function calls, pointers to primal inputs will be null.
      void disableInputPrimalStore() {
        storeInputPrimals = false;
        data->provideInputValues = false;
      }

      /// Do not store primal output values. In function calls, pointers to primal outputs will be null.
      void disableOutputPrimalStore() {
        storeOutputPrimals = false;
        data->provideOutputValues = false;
      }

      /// Add an input value.
      CODI_INLINE void addInput(Type const& input) {
        if (Type::getTape().isActive()) {
          Identifier identifier = input.getIdentifier();
          if (!Type::getTape().isIdentifierActive(identifier)) {
            // Register input values for primal value tapes when they are restored from the tape, otherwise the primal
            // values can not be restored. For a lot of inactive inputs, this can inflate the number of identifiers
            // quite a lot. This is especially true for reuse index tapes.
            if (data->getPrimalsFromPrimalValueVector) {
              Type temp = input;
              Type::getTape().registerInput(temp);
              identifier = temp.getIdentifier();
            }
          }

          data->inputIndices.push_back(identifier);
        }

        // Ignore the setting at this place and the active check,
        // we might need the values for the evaluation.
        if (storeInputOutputForPrimalEval || storeInputPrimals) {
          data->inputValues.push_back(input.getValue());
        }
      }

    private:

      CODI_INLINE void addOutputToData(Type& output) {
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

      /// Add an output value.
      CODI_INLINE void addOutput(Type& output) {
        if (Type::getTape().isActive() || storeInputOutputForPrimalEval) {
          outputValues.push_back(&output);
        }
      }

      /// Add user data. See ExternalFunctionUserData for details.
      template<typename Data>
      CODI_INLINE void addUserData(Data const& data) {
        this->data->userData.addData(data);
      }

      /// Get a reference to the full user data created for this external function.
      /// See ExternalFunctionUserData for details.
      ExternalFunctionUserData& getExternalFunctionUserData() {
        return this->data->userData;
      }

      /// This is intended for primal functions that are implemented with the AD type. It is ensured that no data is
      /// recorded on the tape. All output values are registered as outputs of this external function.
      template<typename FuncObj, typename... Args>
      CODI_INLINE void callPrimalFuncWithADType(FuncObj& func, Args&&... args) {
        bool isTapeActive = Type::getTape().isActive();

        if (isTapeActive) {
          Type::getTape().setPassive();
        }

        func(std::forward<Args>(args)...);

        if (isTapeActive) {
          Type::getTape().setActive();

          for (size_t i = 0; i < outputValues.size(); ++i) {
            addOutputToData(*outputValues[i]);
          }
        }
      }

      /// Call the primal function with the values extracted from the inputs. The output values are set on the specified
      /// outputs and registered as outputs of this external functions.
      CODI_INLINE void callPrimalFunc(PrimalFunc func) {
        if (storeInputOutputForPrimalEval) {
          // Store the primal function in the external function data so that it can be used for primal evaluations of
          // the tape.
          data->primalFunc = func;

          std::vector<Real> y(outputValues.size());

          func(data->inputValues.data(), data->inputValues.size(), y.data(), outputValues.size(), &data->userData);

          // Set the primal values on the output values and add them to the data for the reverse evaluation.
          for (size_t i = 0; i < outputValues.size(); ++i) {
            outputValues[i]->setValue(y[i]);

            if (Type::getTape().isActive()) {
              addOutputToData(*outputValues[i]);
            }
          }
        } else {
          CODI_EXCEPTION(
              "callPrimalFunc() not available if external function helper is initialized with passive function mode "
              "enabled. Use callPrimalFuncWithADType() instead.");
        }
      }

      /// Add the external function to the tape.
      CODI_INLINE void addToTape(ReverseFunc reverseFunc, ForwardFunc forwardFunc = nullptr,
                                 PrimalFunc primalFunc = nullptr) {
        if (Type::getTape().isActive()) {
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
            data->inputValues.shrink_to_fit();
          }

          Type::getTape().pushExternalFunction(
              ExternalFunction<Tape>::create(EvalData::evalRevFuncStatic, data, EvalData::delFunc,
                                             EvalData::evalForwFuncStatic, EvalData::evalPrimFuncStatic));

          data = nullptr;
        } else {
          // Clear the assembled data.
          delete data;
        }

        // Create a new data object for the next call.
        data = new EvalData();
        outputValues.clear();
      }
  };
}
