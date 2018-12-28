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

#include "../adjointInterface.hpp"
#include "../adjointInterfaceImpl.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Data and functions for external functions.
   *
   * The structures stores all elements for the function evaluation.
   *
   * The function itself and a handle to the user data for the
   * function is stored. The user data changes its ownership when it is provided to
   * the external function, therefore a method for deleting the data is stored.
   *
   * The data will not be deleted in the destructor as the structure is considered
   * a POD type.
   */
  struct ExternalFunction {
    /**
     * @brief Definition for the user function.
     *
     * The first parameter gives a pointer to the tape that calls the external function.
     * It can either be cast to the correct type or to the ReverseTapeInterface.
     *
     * The second parameter gives a pointer to the data from the user.
     *
     * The third paramter provides a universal adapter for the adjoint and/or primal values.
     * It can be used instead of the tape to get the adjoint values in an abstract fashion.
     * If a custiom adjoint vector is used in the tape evaluation, then this
     * interface needs to be used. The value needs to be cast to the AdjointInterface class
     * that uses the correct floating point replacement from the tape.
     */
    typedef void (*CallFunction)(void* tape, void* data, void* ra);
    /**
     * @brief Definition for the delete function of the user data.
     *
     * The first parameter gives a pointer to the tape that calls the external function.
     * It can either be cast to the correct type or to the ReverseTapeInterface.
     *
     * The second parameter gives a pointer to the data from the user.
     */
    typedef void (*DeleteFunction)(void* tape, void* data);

  private:
    /** @brief The reverse function given by the user. */
    CallFunction funcReverse;
    /** @brief The forward function given by the user. */
    CallFunction funcForward;
    /** @brief The primal function given by the user. */
    CallFunction funcPrimal;
    /** @brief The delete function for the user data. */
    DeleteFunction deleteCheckpoint;

    /** @brief The data for the function. */
    void* data;

  public:
    /**
     * @brief Needed to construct arrays.
     */
    ExternalFunction() {}
    /**
     * @brief Create the structure with all data.
     *
     * @param[in]      funcReverse  The user function which is called with the data for a reverse evaluation.
     * @param[in]      funcForward  The user function which is called with the data for a forward evaluation.
     * @param[in]       funcPrimal  The user function which is called with the data for a primal evaluation.
     * @param[in,out]         data  The data for the user function.
     * @param[in] deleteCheckpoint  The function which deletes the user data. The function handle can be NULL.
     */
    ExternalFunction(CallFunction funcReverse, CallFunction funcForward, CallFunction funcPrimal, void* data, DeleteFunction deleteCheckpoint) :
      funcReverse(funcReverse),
      funcForward(funcForward),
      funcPrimal(funcPrimal),
      deleteCheckpoint(deleteCheckpoint),
      data(data){}

    /**
     * @brief Delete the data provided by the user.
     *
     * If the function handle is NULL the data is not deleted.
     *
     * @param[in,out] tape  The tape that calls the function.
     */
    void deleteData(void* tape) {
      if (deleteCheckpoint != NULL){
        deleteCheckpoint(tape, data);
        data = NULL;
      }
    }

    /**
     * @brief Call the user function with the user data as an argument.
     *
     * @param[in,out] tape  The tape that calls the function.
     * @param[in,out]   ra  The interface to the used adjoint vector.
     */
    void evaluateReverse(void* tape, void* ra) {
      if(NULL != funcReverse) {
        funcReverse(tape, data, ra);
      }
    }

    /**
     * @brief Call the forward user function with the user data as an argument.
     *
     * @param[in,out] tape  The tape that calls the function.
     * @param[in,out]   ra  The interface to the used adjoint vector.
     */
    void evaluateForward(void* tape, void* ra) {
      if(NULL != funcForward) {
        funcForward(tape, data, ra);
      }
    }

    /**
     * @brief Call the primal user function with the user data as an argument.
     *
     * @param[in,out] tape  The tape that calls the function.
     * @param[in,out]   ra  The interface to the used adjoint vector.
     */
    void evaluatePrimal(void* tape, void* ra) {
      if(NULL != funcPrimal) {
        funcPrimal(tape, data, ra);
      }
    }
  };

  /**
   * @brief Data and functions for strongly typed external functions.
   *
   * The structures stores all elements for the function evaluation.
   *
   * The function itself and a pointer to the user data for the
   * function is stored. The user data changes its ownership when it is provided to
   * the external function, therefore a method for deleting the data is stored.
   *
   * The data will not be deleted in the destructor as the structure is considered
   * a POD type.
   *
   * @tparam Tape The type of the tape that calls the external function.
   * @tparam Data The type of the data used in the external function.
   */
  template<typename Tape, typename Data>
  class ExternalFunctionDataHelper {
  public:
    /**
     * @brief Definition for the user function.
     *
     * The first parameter gives a pointer to the tape that calls the external function.
     *
     * The second parameter gives a pointer to the data from the user.
     *
     * The third paramter provides a universal adapter for the adjoint and/or primal values.
     * It can be used instead of the tape to get the adjoint values in an abstract fashion.
     * If a custiom adjoint vector is used in the tape evaluation, then this
     * interface needs to be used.
     */
    typedef void (*CallFunction)(Tape*, Data*, AdjointInterface<typename Tape::Real, typename Tape::Index>*);
    /**
     * @brief Definition for the delete function of the user data.
     *
     * The first parameter gives a pointer to the tape that calls the external function.
     *
     * The second parameter gives a pointer to the data from the user.
     */
    typedef void (*DeleteFunction)(Tape*, Data*);

  private:
    /** @brief The reverse function given by the user. */
    CallFunction funcReverse;
    /** @brief The forward function given by the user. */
    CallFunction funcForward;
    /** @brief The primal function given by the user. */
    CallFunction funcPrimal;
    /** @brief The delete function for the user data. */
    DeleteFunction deleteData;

    /** @brief The data for the function. */
    Data* data;

  private:
    /**
     * @brief Create the structure with all data.
     *
     * @param[in]      funcReverse  The user function which is called with the data for a reverse evaluation.
     * @param[in]      funcForward  The user function which is called with the data for a forward evaluation.
     * @param[in]       funcPrimal  The user function which is called with the data for a primal evaluation.
     * @param[in,out]         data  The data for the user function.
     * @param[in] deleteCheckpoint  The function which deletes the user data. The function handle can be NULL.
     */
    ExternalFunctionDataHelper(CallFunction funcReverse, CallFunction funcForward, CallFunction funcPrimal, Data* data, DeleteFunction deleteData) :
      funcReverse(funcReverse),
      funcForward(funcForward),
      funcPrimal(funcPrimal),
      deleteData(deleteData),
      data(data){}

    /**
     * @brief Helper function which is called from the #ExternalFunction structure.
     *
     * The method calls the user function with the data stored in the data handle.
     *
     * @param[in,out] tape  The tape which calls the function
     * @param[in,out] data  The void handle from the #ExternalFunction is a handle to a #ExternalFunctionDataHelper object.
     */
    static void callFunctionReverse(void* tape, void* data, void* ra) {
      ExternalFunctionDataHelper<Tape, Data>* castData = cast(data);
      castData->funcReverse((Tape*)tape, castData->data, (AdjointInterface<typename Tape::Real, typename Tape::Index>*)ra);
    }

    /**
     * @brief Helper function which is called from the #ExternalFunction structure.
     *
     * The method calls the user function with the data stored in the data handle.
     *
     * @param[in,out] tape  The tape which calls the function
     * @param[in,out] data  The void handle from the #ExternalFunction is a handle to a #ExternalFunctionDataHelper object.
     */
    static void callFunctionForward(void* tape, void* data, void* ra) {
      ExternalFunctionDataHelper<Tape, Data>* castData = cast(data);
      castData->funcForward((Tape*)tape, castData->data, (AdjointInterface<typename Tape::Real, typename Tape::Index>*)ra);
    }

    /**
     * @brief Helper function which is called from the #ExternalFunction structure.
     *
     * The method calls the user function with the data stored in the data handle.
     *
     * @param[in,out] tape  The tape which calls the function
     * @param[in,out] data  The void handle from the #ExternalFunction is a handle to a #ExternalFunctionDataHelper object.
     */
    static void callFunctionPrimal(void* tape, void* data, void* ra) {
      ExternalFunctionDataHelper<Tape, Data>* castData = cast(data);
      castData->funcPrimal((Tape*)tape, castData->data, (AdjointInterface<typename Tape::Real, typename Tape::Index>*)ra);
    }

    /**
     * @brief Helper function which is called from the #ExternalFunction structure.
     *
     * The method deletes the user data for the user function and deletes the  the user function with the data stored in the data handle.
     *
     * @param[in,out] tape  The tape which calls the function
     * @param[in,out] data  The void handle from the #ExternalFunction is a handle to a #ExternalFunctionDataHelper object.
     */
    static void deleteFunction(void* tape, void* data) {
      ExternalFunctionDataHelper<Tape, Data>* castData = cast(data);
      castData->deleteData((Tape*)tape,castData->data);

      // delete self
      delete castData;
    }

    /**
     * @brief Cast the data handle to a ExternalFunctionDataHelper pointer.
     * @param[in,out] data  The void handle to the data.
     * @return The converted void pointer.
     */
    static ExternalFunctionDataHelper<Tape, Data>* cast(void* data) {
      return (ExternalFunctionDataHelper<Tape, Data>*)data;
    }

  public:

    /**
     * @brief Create an ExternalFunction object with strong typed data.
     *
     * @param[in] funcReverse  The reverse user function which is called with the data.
     * @param[in] funcForward  The reverse user function which is called with the data.
     * @param[in]  funcPrimal  The reverse user function which is called with the data.
     * @param[in,out]    data  The data for the user function.
     * @param[in]  deleteData  The function which deletes the user data. The function handle can be NULL.
     *
     * @return The ExternalFunction object with strong typed data.
     */
    static CODI_INLINE ExternalFunction createHandle(CallFunction funcReverse, CallFunction funcForward, CallFunction funcPrimal, Data* data, DeleteFunction deleteData) {
      ExternalFunctionDataHelper<Tape, Data>* functionHelper = new ExternalFunctionDataHelper<Tape, Data>(funcReverse, funcForward, funcPrimal, data, deleteData);
      return ExternalFunction(ExternalFunctionDataHelper<Tape, Data>::callFunctionReverse, ExternalFunctionDataHelper<Tape, Data>::callFunctionForward, ExternalFunctionDataHelper<Tape, Data>::callFunctionPrimal, functionHelper, ExternalFunctionDataHelper<Tape, Data>::deleteFunction);
    }
  };
}
