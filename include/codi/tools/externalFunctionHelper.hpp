/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief The data object for the external function helper.
   *
   * It stores all the data for the external function and provides the static functions for the registration on the
   * tapes.
   *
   * @tparam CoDiType  The CoDiPack type that is used in the application.
   */
  template<typename CoDiType>
  class ExternalFunctionData {
    public:

      typedef typename CoDiType::Real Real; /**< The floating point calculation type in the CoDiPack types. */
      typedef typename CoDiType::GradientData GradientData; /**< The type for the identification of gradients. */
      typedef typename CoDiType::GradientValue GradientValue; /**< The type for the gradient computation */

      /** The type of the tape implementation. */
      typedef typename CoDiType::TapeType Tape;

      /**
       * @brief The function for the reverse evaluation.
       *
       * If the primal function evaluates y = f(x), then this function evaluates
       * \f[\bar x = \frac{d f}{d x}(x) \bar y\f]
       *
       * If the user disabled the storing of the primal values. Then the corresponding vectors are null pointers.
       */
      typedef void (*ReverseFunc)(const Real* x, Real* x_b, size_t m, const Real* y, const Real* y_b, size_t n, DataStore* d);

      std::vector<GradientData> inputIndices; /**< The storage for the identifiers of the input values. */
      std::vector<GradientData> outputIndices; /**< The storage for the identifiers of the output values. */

      std::vector<Real> inputValues; /**< The storage for the primal values of the input values. */
      std::vector<Real> outputValues; /**< The storage for the primal values of the output values. */
      std::vector<Real> oldPrimals; /**< The old value in a primal value tape, that are overwritten by the output values. */

      ReverseFunc revFunc; /**< The reverse function provided by the user. */

      DataStore userData; /**< The data manager for the user data. */

      /**
       * @brief The delete function that is registered on the tape.
       *
       * It calls delete on the data object.
       *
       * @param[in] t  unused
       * @param[in] d  An instance of this class.
       */
      static void delFunc(void* t, void* d) {
        CODI_UNUSED(t);

        ExternalFunctionData<CoDiType>* data = (ExternalFunctionData<CoDiType>*)d;

        delete data;
      }

      /**
       * @brief Reverse evaluation function that is registered on the tape.
       *
       * The method casts the data object to an instance of this class and calls the evalRevFunc.
       *
       * @param[in,out]  t  The tape which evaluates this function.
       * @param[in,out]  d  An instance of this class.
       * @param[in,out] ra  The helper structure for the access to the adjoint and primal vector.
       */
      static void evalRevFuncStatic(void* t, void* d, void* ra) {
        ExternalFunctionData<CoDiType>* data = (ExternalFunctionData<CoDiType>*)d;

        data->evalRevFunc((Tape*)t, (AdjointInterface<Real, GradientData>*)ra);
      }

      /**
       * @brief The reverse evaluation function.
       *
       * This function retrieves the adjoint values of the output. Afterwards the user defined evaluation function is
       * called. The adjoint of the input values is then used for the update of the tape adjoints.
       *
       * If the adjoint interface specifies a vector mode the function is evaluated multiple times.
       *
       * @param[in,out]  t  The tape which evaluates this function.
       * @param[in,out] ra  The helper structure for the access to the adjoint and primal vector.
       */
      void evalRevFunc(Tape* t, AdjointInterface<Real, GradientData>* ra) {
        CODI_UNUSED(t);

        Real* x_b = new Real[inputIndices.size()];
        Real* y_b = new Real[outputIndices.size()];

        for(size_t dim = 0; dim < ra->getVectorSize(); ++dim) {

          for(size_t i = 0; i < outputIndices.size(); ++i) {
            y_b[i] = ra->getAdjoint(outputIndices[i], dim);
            ra->resetAdjoint(outputIndices[i], dim);
          }

          revFunc(inputValues.data(), x_b, inputIndices.size(), outputValues.data(), y_b, outputIndices.size(), &userData);

          for(size_t i = 0; i < inputIndices.size(); ++i) {
            ra->updateAdjoint(inputIndices[i], dim, x_b[i]);
          }
        }

        if(Tape::RequiresPrimalReset) {
          for(size_t i = 0; i < outputIndices.size(); ++i) {
            ra->setPrimal(outputIndices[i], oldPrimals[i]);
          }
        }

        delete [] x_b;
        delete [] y_b;
      }
  };

  /**
   * @brief Helper class for the simple implementation of an external function in CoDiPack.
   *
   * The class helps the user the handle parts where the CoDiPack types can not be applied or where a more efficient
   * gradient computation is available.
   *
   * Lets assume that the function which needs to be handled is called 'func' and the mathematical model for func is
   * \f[ y = f(x) \eqdot \f]
   * In order to add an external function for the 'func' the user needs to implement a function that computes the
   * adjoint model of the function which is:
   * \f[ \bar{x} = \frac{df}{dx}^T(x)\cdot \bar{y}. \f]
   *
   * Assume the user has implemented this adjoint model in the function 'func_adj'. This function needs to have the
   * layout defined by ExternalFunctionData::ReverseFunc which is
   *
   * revFunc(x, x_b, m, y, y_b, n, d)
   *
   * where 'm' is the size of 'x' and 'x_b' and 'n' is the size of 'y' and 'y_b'.
   * 'x_b' corresponds to the the \f$ \bar x \f$ value. The same is true for 'y_b'.
   * 'd' is a DataStore object that is used to store used defined data.
   *
   * The first operation mode for the external function helper is the one where the code can not be handled with AD,
   * here the user needs to define a function 'func_prim' which evaluates the function with in the unmodified way (The is
   * with no AD tool). The evaluation schedule is:
   *
   * \code{.cpp}
   * ExternalFunctionHelper<CoDiType> eh;
   *
   * for each xVar in x
   *   eh.addInput(xVar)
   *
   * for each yVar in y
   *   eh.addOutput(yVar)
   *
   * // if user data is required
   * for each data in d
   *   eh.addUserData(data)
   *
   * eh.callPrimalFunc(func_prim);
   * eh.addToTape(func_adj);
   * \endcode
   *
   * This will ensure, that the adjoint data of the external function is correctly handled. Before any function is called
   * on the helper structure the user can decide to disable the storing of the primal values for 'x' or 'y'. The
   * respective methods are ExternalFunctionHelper::disableInputPrimalStore and ExternalFunctionHelper::disableOutputPrimalStore.
   *
   * The second mode of operation can be used if the code section for the external function is already differentiated
   * with a CoDiPack type but the there is a more efficient way for the computation of the derivatives available.
   * The evaluation schedule for the second use case is:
   *
   * \code{.cpp}
   * ExternalFunctionHelper<CoDiType> eh(true);
   *
   * for each xVar in x
   *   eh.addInput(xVar)
   *
   * eh.callPassiveFunc(func, <args for function>);
   *
   * for each yVar in y
   *   eh.addOutput(yVar)
   *
   * // if user data is required
   * for each data in d
   *   eh.addUserData(data)
   *
   * eh.addToTape(func_adj);
   * \endcode
   *
   * It is important, that here the function with the CoDiPack type is called. The helper will ensure that no tape is
   * stored for the call to the function.
   *
   * The ExternalFunctionHelper works with all tapes. It is also able to handle situations where the tape is currently
   * not recording. All necessary operations are performed in such a case but no external function is recorded.
   *
   * @tparam CoDiType  This needs to be one of the CoDiPack types defined through an ActiveReal
   */
  template<typename CoDiType>
  class ExternalFunctionHelper {
    public:
      typedef typename CoDiType::Real Real; /**< The floating point calculation type in the CoDiPack types. */
      typedef typename CoDiType::GradientData GradientData; /**< The type for the identification of gradients. */
      typedef typename CoDiType::GradientValue GradientValue; /**< The type for the gradient computation */

      /** The type of the tape implementation. */
      typedef typename CoDiType::TapeType Tape;

      /**
       * @brief The function for the primal evaluation.
       *
       * This function needs to be provided if use case one is used from the class documentation.
       */
      typedef void (*PrimalFunc)(const Real* x, size_t m, Real* y, size_t n, DataStore* d);

      /** @brief Forward definition of the reverse evaluation function. */
      typedef typename ExternalFunctionData<CoDiType>::ReverseFunc ReverseFunc;

      /**
       * @brief Pointer array to the output values.
       *
       * In use case one the pointers to the output values are stored, since they need to be modified after the primal
       * function call.
       */
      std::vector<CoDiType*> outputValues;

      bool storeInputPrimals; /**< If false the storing of the primal input values is omitted. */
      bool storeOutputPrimals;  /**< If false the storing of the primal output values is omitted. */
      bool isPassiveExtFunc;  /**< Flag if a passive function was called by the user. */
      bool isTapeActive;  /**< General flag if the tape was active during the creation of the external function. */

      /**
       * @brief The data object that is stored on the tape.
       */
      ExternalFunctionData<CoDiType> *data;

      /**
       * @brief Initializes the structure also determines if the tape is currently recording. The recording state
       * may not be changed by the user until the external function is finished.
       *
       * @param[in] passiveExtFunc Parameter if the passive evaluation mode is used.
       */
      ExternalFunctionHelper(bool passiveExtFunc = false) :
        outputValues(),
        storeInputPrimals(true),
        storeOutputPrimals(true),
        isPassiveExtFunc(passiveExtFunc),
        isTapeActive(CoDiType::getGlobalTape().isActive()),
        data(nullptr) {
        data = new ExternalFunctionData<CoDiType>();
      }

      /**
       * @brief If the tape is not recording then the data structure is deleted here.
       */
      ~ExternalFunctionHelper() {
        if(!isTapeActive) {
          delete data;
        }
      }

      /**
       * @brief Disables the storing of the primal values for the reverse function.
       *
       * The pointer to x will be a null pointer.
       */
      void disableInputPrimalStore() {
        storeInputPrimals = false;
      }

      /**
       * @brief Disables the storing of the primal values for the reverse function.
       *
       * The pointer to y will be a null pointer.
       */
      void disableOutputPrimalStore() {
        storeOutputPrimals = false;
      }

      /**
       * @brief Add a CoDiPack type to the input set of the function.
       *
       * @param[in] input  The input value which is added to the input set.
       */
      void addInput(const CoDiType& input) {
        if(isTapeActive) {
          data->inputIndices.push_back(input.getGradientData());
        }

        // ignore the setting at this place and the active check
        // We might need the values for the evaluation.
        if (!isPassiveExtFunc || storeInputPrimals){
          data->inputValues.push_back(input.getValue());
        }
      }

    private:

      /**
       * @brief Helper function for the correct adding of an output value.
       *
       * The value is registered on the tape and if necessary the old primal values are stored.
       *
       * @param[in] output  The output value which is added to the output set.
       */
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

      /**
       * @brief Adds a CoDiPack type to the output set of the function.
       *
       * For a primal function call (operation mode one) the pointer to the
       * value is stored. The pointer is used in callPrimalFunc to modify the value.
       *
       * For a passive function call (operation mode two) the value is directly modified.
       *
       * @param[in,out] output  The output value which is added to the output set.
       */
      void addOutput(CoDiType& output) {
        if(isTapeActive) {
          if(isPassiveExtFunc) {
            addOutputToData(output);
          } else {
            outputValues.push_back(&output);
          }
        }
      }

      /**
       * @brief Adds data which can be used in the primal and reverse function calls.
       *
       * The data is added to a DataStore object, which is given to the primal and reverse function calls.
       *
       * @param[in] data  The data for the primal and reverse evaluation.
       *
       * @tparam Data  The type of the data for the primal and reverse evaluation.
       */
      template<typename Data>
      void addUserData(const Data& data) {
        this->data->userData.addData(data);
      }

      /**
       * @brief Helper function to get direct access to the data store.
       *
       * The direct access can be used for the more advanced storing functions of the data store.
       *
       * @return The store object which is given to the primal and reverse evaluation.
       */
      DataStore& getDataStore() {
        return this->data->userData;
      }

      /**
       * @brief Function caller for the second operation mode of the external function helper.
       *
       * All input values need to be added before this function is called. The tape will be set to passive if necessary
       * and reactivated afterwards.
       *
       * After the function is called all output values need to be added.
       *
       * @param[in]     func  The function that defines the code region which should not be recorded and is handled in
       *                      a special way during the reverse evaluation.
       * @param[in,out] args  The arguments for the function.
       *
       * @tparam FuncObj  Needs to be a function object.
       * @tparam    Args  The argument types for the function object.
       */
      template<typename FuncObj, typename ... Args>
      void callPassiveFunc(FuncObj& func, Args&& ... args) {

        if(isTapeActive) {
          CoDiType::getGlobalTape().setPassive();
        }

        func(std::forward<Args>(args)...);

        if(isTapeActive) {
          CoDiType::getGlobalTape().setActive();
        }
      }

      /**
       * @brief Function caller for the first operation mode of the external function helper.
       *
       * All input and output values as well as the user data need to be added before this function is called. All
       * primal values are extracted and given to the primal implementation.
       *
       * @param[in] func  The implementation for the primal evaluation.
       */
      void callPrimalFunc(PrimalFunc func) {
        if (!isPassiveExtFunc){
          Real* y = new Real[outputValues.size()];

          func(data->inputValues.data(), data->inputValues.size(), y, outputValues.size(), &data->userData);

          // ok now set the primal values on the output values and add them to the data for the reverse evaluation
          for(size_t i = 0; i < outputValues.size(); ++i) {
            outputValues[i]->setValue(y[i]);

            addOutputToData(*outputValues[i]);
          }

          delete [] y;
        } else {
          std::cerr << "callPrimalFunc() not available if external function helper is initialized with passive function mode enabled. Use callPassiveFunc() instead." << std::endl;
          exit(-1);
        }
      }

      /**
       * @brief This function needs to be called as the last function. It will finally add the external function to the
       * tape such that the specialized reverse implementation is called during the reverse interpretation.
       *
       * @param[in] func  The logic for the reverse implementation.
       */
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
