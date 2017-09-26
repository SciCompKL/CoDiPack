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
     */
    typedef void (*CallFunction)(void* tape, void* data);
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
    /** @brief The function given by the user. */
    CallFunction func;
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
     * @param[in]             func  The user function which is called by the data.
     * @param[in,out]         data  The data for the user function.
     * @param[in] deleteCheckpoint  The function which deletes the user data. The function handle can be NULL.
     */
    ExternalFunction(CallFunction func, void* data, DeleteFunction deleteCheckpoint) :
      func(func),
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
     */
    void evaluate(void* tape) {
      if(NULL != func) {
        func(tape, data);
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
     */
    typedef void (*CallFunction)(Tape*, Data*);
    /**
     * @brief Definition for the delete function of the user data.
     *
     * The first parameter gives a pointer to the tape that calls the external function.
     *
     * The second parameter gives a pointer to the data from the user.
     */
    typedef void (*DeleteFunction)(Tape*, Data*);

  private:
    /** @brief The function given by the user. */
    CallFunction func;
    /** @brief The delete function for the user data. */
    DeleteFunction deleteData;

    /** @brief The data for the function. */
    Data* data;

  private:
    /**
     * @brief Create the structure with all data.
     *
     * @param[in]             func  The user function which is called by the data.
     * @param[in,out]         data  The data for the user function.
     * @param[in] deleteCheckpoint  The function which deletes the user data. The function handle can be NULL.
     */
    ExternalFunctionDataHelper(CallFunction func, Data* data, DeleteFunction deleteData) :
      func(func),
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
    static void callFunction(void* tape, void* data) {
      ExternalFunctionDataHelper<Tape, Data>* castData = cast(data);
      castData->func((Tape*)tape, castData->data);
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
     * @param[in]       func  The user function which is called by the data.
     * @param[in,out]   data  The data for the user function.
     * @param[in] deleteData  The function which deletes the user data. The function handle can be NULL.
     *
     * @return The ExternalFunction object with strong typed data.
     */
    static CODI_INLINE ExternalFunction createHandle(CallFunction func, Data* data, DeleteFunction deleteData) {
      ExternalFunctionDataHelper<Tape, Data>* functionHelper = new ExternalFunctionDataHelper<Tape, Data>(func, data, deleteData);
      return ExternalFunction(ExternalFunctionDataHelper<Tape, Data>::callFunction, functionHelper, ExternalFunctionDataHelper<Tape, Data>::deleteFunction);
    }
  };
}
