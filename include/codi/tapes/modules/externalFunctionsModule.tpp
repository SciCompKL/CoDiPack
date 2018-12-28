/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2018 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

/*
 * In order to include this file the user has to define the preprocessor macro CHILD_VECTOR_TYPE, CHILD_VECTOR_NAME,
 * VECTOR_TYPE and TAPE_NAME.
 *
 * CHILD_VECTOR_TYPE defines the type of the nested vector for the data vector.
 * CHILD_VECTOR_NAME defines the member name of the nested vector instantiation.
 * VECTOR_TYPE defines the type of the data vector.
 *
 * All these macros are undefined at the end of the file.
 *
 * TAPE_NAME defines the type name of the tape and is not undefined at the end of the file.
 *
 * The module defines the structures extFuncVector.
 * The module defines the types ExtFuncChildVector, ExtFuncChildPosition, ExtFuncVector, ExtFuncChunk,
 * ExtFuncPosition.
 *
 * It defines the methods setExternalFunctionChunkSize, pushExternalFunctionHandle, pushExternalFunction,
 * printExtFuncStatistics from the TapeInterface and ReverseTapeInterface.
 *
 * It defines the methods getExtFuncPosition, getExtFuncZeroPosition, resetExtFunc, evaluateExtFunc, evaluateExtFuncForward as interface functions for the
 * including class.
 */

#ifndef TAPE_NAME
  #error Please define the name of the tape.
#endif

#ifndef CHILD_VECTOR_TYPE
  #error Please define the type of the child vector
#endif
#ifndef CHILD_VECTOR_NAME
  #error Please define the name of the child vector
#endif
#ifndef VECTOR_TYPE
  #error Please define the name of the chunk vector type.
#endif

  private:

  // ----------------------------------------------------------------------
  // All definitions of the module
  // ----------------------------------------------------------------------

    /** @brief The child vector for the external function data vector. */
    typedef CHILD_VECTOR_TYPE ExtFuncChildVector;

    /** @brief The position type of the external function child vector */
    typedef typename ExtFuncChildVector::Position ExtFuncChildPosition;

    /** @brief The vector for the external function data. */
    typedef VECTOR_TYPE ExtFuncVector;
    /** @brief The data for the external functions. */
    typedef typename ExtFuncVector::ChunkType ExtFuncChunk;

    /** @brief The position type of the external function module. */
    typedef typename ExtFuncVector::Position ExtFuncPosition;

    /** @brief The data for the external functions. */
    ExtFuncVector extFuncVector;

  // ----------------------------------------------------------------------
  // Private function of the module
  // ----------------------------------------------------------------------

    /**
     * @brief Private common method to add to the external function stack.
     *
     * @param[in] function The external function structure to push.
     */
    void pushExternalFunctionHandle(const ExternalFunction& function){
      extFuncVector.reserveItems(1);
      extFuncVector.setDataAndMove(function, CHILD_VECTOR_NAME.getPosition());
    }

    /**
     * @brief Function object for the deletion of external functions.
     */
    struct ExtFuncDeleter {
      TAPE_NAME& tape;

      /**
       * @brief Create the function object.
       *
       * @param[in,out]     tape  The reference to the actual tape.
       */
      ExtFuncDeleter(TAPE_NAME& tape) :
        tape(tape){}

      /**
       * @brief The operator deletes the external function object.
       *
       * @param[in]     extFunc  The external function object.
       * @param[in] endInnerPos  The position were the external function object was stored.
       */
      void operator () (ExternalFunction* extFunc, const ExtFuncChildPosition* endInnerPos) {
        CODI_UNUSED(endInnerPos);

        /* we just need to call the delete function */
        extFunc->deleteData(&tape);
      }
    };

  private:

  // ----------------------------------------------------------------------
  // Private function for the communication with the including class
  // ----------------------------------------------------------------------

    /**
     * @brief Get the position from the external function vector.
     *
     * The position contains also the positions information for all nested vectors.
     *
     * @return The position for the external function vector.
     */
    ExtFuncPosition getExtFuncPosition() const {
      return extFuncVector.getPosition();
    }

    /**
     * @brief Get the zero position from the external function vector.
     *
     * The position contains also the positions information for all nested vectors.
     *
     * @return The zero position for the external function vector.
     */
    ExtFuncPosition getExtFuncZeroPosition() const {
      return extFuncVector.getZeroPosition();
    }

    /**
     * @brief Reset the external function module to the position.
     *
     * The reset will also reset the vector and therefore all nested vectors.
     *
     * @param[in] pos  The position to which the tape is reset.
     */
    void resetExtFunc(const ExtFuncPosition& pos) {
      ExtFuncDeleter deleter(*this);

      extFuncVector.forEachReverse(getExtFuncPosition(), pos, deleter);

      // reset will be done iteratively through the vectors
      extFuncVector.reset(pos);
    }

    /**
     * @brief Evaluate a part of the external function vector.
     *
     * It has to hold start <= end.
     *
     * It calls the primal evaluation method for the statement vector.
     *
     * @param[in]       start  The starting point for the external function vector.
     * @param[in]         end  The ending point for the external function vector.
     * @param[in]        func  The function that is evaluated before and after each external function call.
     * @param[in,out]     obj  The object on which the function is evaluated.
     * @param[in,out]    args  The arguments for the evaluation.
     *
     * @tparam     Args  The types of the other arguments.
     * @tparam Function  A function object which is called with the nested start and end positions.
     * @tparam      Obj  The object on which the function is called.
     */
    template<typename Function, typename Obj, typename ... Args>
    void evaluateExtFuncPrimal(const ExtFuncPosition& start, const ExtFuncPosition &end, const Function& func,
                               Obj& obj, AdjointInterface<Real, Index>* adjointInterface, Args&&... args){
      ExtFuncChildPosition curInnerPos = start.inner;
      auto evalFunc = [&] (ExternalFunction* extFunc, const ExtFuncChildPosition* endInnerPos,
                           AdjointInterface<Real, GradientData>* adjointInterface, Args&&... args) {

        (obj.*func)(curInnerPos, *endInnerPos, std::forward<Args>(args)...);

        extFunc->evaluatePrimal(&obj, adjointInterface);

        curInnerPos = *endInnerPos;

      };
      extFuncVector.forEachForward(start, end, evalFunc, adjointInterface, std::forward<Args>(args)...);

      // Iterate over the reminder also covers the case if there have been no external functions.
      (obj.*func)(curInnerPos, end.inner, std::forward<Args>(args)...);
    }

    /**
     * @brief Evaluate a part of the external function vector.
     *
     * It has to hold start >= end.
     *
     * It calls the evaluation method for the statement vector.
     *
     * @param[in]                start  The starting point for the external function vector.
     * @param[in]                  end  The ending point for the external function vector.
     * @param[in]                 func  The function that is evaluated before and after each external function call.
     * @param[in,out]              obj  The object on which the function is evaluated.
     * @param[in,out] adjointInterface  Interface for accessing the adjoint vector for this evaluation.
     * @param[in,out]             args  The arguments for the evaluation.
     *
     * @tparam     Args  The types of the other arguments.
     * @tparam Function  A function object which is called with the nested start and end positions.
     * @tparam      Obj  The object on which the function is called.
     */
    template<typename Function, typename Obj, typename ... Args>
    CODI_INLINE void evaluateExtFunc(const ExtFuncPosition& start, const ExtFuncPosition &end, const Function& func,
                                     Obj& obj, AdjointInterface<Real, Index>* adjointInterface, Args&&... args){


      // Iterate over the reminder also covers the case if there have been no external functions.
      ExtFuncChildPosition curInnerPos = start.inner;
      auto evalFunc = [&] (ExternalFunction* extFunc, const ExtFuncChildPosition* endInnerPos,
                           AdjointInterface<Real, GradientData>* adjointInterface, Args&&... args) {

         (obj.*func)(curInnerPos, *endInnerPos, std::forward<Args>(args)...);

         extFunc->evaluateReverse(&obj, adjointInterface);

         curInnerPos = *endInnerPos;

       };
       extFuncVector.forEachReverse(start, end, evalFunc, adjointInterface, std::forward<Args>(args)...);

       // Iterate over the reminder also covers the case if there have been no external functions.
       (obj.*func)(curInnerPos, end.inner, std::forward<Args>(args)...);
    }

    /**
     * @brief Evaluate a part of the external function vector.
     *
     * It has to hold start <= end.
     *
     * It calls the forward evaluation method for the statement vector.
     *
     * @param[in]                start  The starting point for the external function vector.
     * @param[in]                  end  The ending point for the external function vector.
     * @param[in]                 func  The function that is evaluated before and after each external function call.
     * @param[in,out]              obj  The object on which the function is evaluated.
     * @param[in,out] adjointInterface  Interface for accessing the adjoint vector for this evaluation.
     * @param[in,out]             args  The arguments for the evaluation.
     *
     * @tparam     Args  The types of the other arguments.
     * @tparam Function  A function object which is called with the nested start and end positions.
     * @tparam      Obj  The object on which the function is called.
     */
    template<typename Function, typename Obj, typename ... Args>
    void evaluateExtFuncForward(const ExtFuncPosition& start, const ExtFuncPosition &end, const Function& func,
                                Obj& obj, AdjointInterface<Real, Index>* adjointInterface, Args&&... args){
      ExtFuncChildPosition curInnerPos = start.inner;
      auto evalFunc = [&] (ExternalFunction* extFunc, const ExtFuncChildPosition* endInnerPos,
                           AdjointInterface<Real, GradientData>* adjointInterface, Args&&... args) {

        (obj.*func)(curInnerPos, *endInnerPos, std::forward<Args>(args)...);

        extFunc->evaluateForward(&obj, adjointInterface);

        curInnerPos = *endInnerPos;

      };
      extFuncVector.forEachForward(start, end, evalFunc, adjointInterface, std::forward<Args>(args)...);

      // Iterate over the reminder also covers the case if there have been no external functions.
      (obj.*func)(curInnerPos, end.inner, std::forward<Args>(args)...);
    }

  public:

  // ----------------------------------------------------------------------
  // Public function from the TapeInterface and ReverseTapeInterface
  // ----------------------------------------------------------------------

    /**
     * @brief Set the size of the external function data chunks.
     *
     * @param[in] extChunkSize The new size for the external function data chunks.
     */
    void setExternalFunctionChunkSize(const size_t& extChunkSize) {
      extFuncVector.setChunkSize(extChunkSize);
    }

    /**
     * @brief Add an external function with a void handle as user data.
     *
     * The data handle provided to the tape is considered in possession of the tape. The tape will now be responsible to
     * free the handle. For this it will use the delete function provided by the user.
     *
     * @param[in]         extFunc  The external function which is called by the tape.
     * @param[in,out]        data  The data for the external function. The tape takes ownership over the data.
     * @param[in]         delData  The delete function for the data.
     * @param[in]  extFuncForward  The external function which is called by the tape for a forward evaluation [Default=null].
     * @param[in]   extFuncPrimal  The external function which is called by the tape for a primal evaluation [Default=null].
     */
    void pushExternalFunctionHandle(ExternalFunction::CallFunction extFunc,
                                    void* data,
                                    ExternalFunction::DeleteFunction delData,
                                    ExternalFunction::CallFunction extFuncForward = nullptr,
                                    ExternalFunction::CallFunction extFuncPrimal = nullptr){
      ENABLE_CHECK (OptTapeActivity, isActive()){
        pushExternalFunctionHandle(ExternalFunction(extFunc, extFuncForward, extFuncPrimal, data, delData));
      }
    }


    /**
     * @brief Add an external function with a specific data type.
     *
     * The data pointer provided to the tape is considered in possession of the tape. The tape will now be responsible to
     * free the data. For this it will use the delete function provided by the user.
     *
     * @param[in]         extFunc  The external function which is called by the tape.
     * @param[in,out]        data  The data for the external function. The tape takes ownership over the data.
     * @param[in]         delData  The delete function for the data.
     * @param[in]  extFuncForward  The external function which is called by the tape for a forward evaluation [Default=null].
     * @param[in]   extFuncPrimal  The external function which is called by the tape for a primal evaluation [Default=null].
     */
    template<typename Data>
    void pushExternalFunction(typename ExternalFunctionDataHelper<TAPE_NAME<TapeTypes>, Data>::CallFunction extFunc,
                              Data* data, typename ExternalFunctionDataHelper<TAPE_NAME<TapeTypes>, Data>::DeleteFunction delData,
                              typename ExternalFunctionDataHelper<TAPE_NAME<TapeTypes>, Data>::CallFunction extFuncForward = nullptr,
                              typename ExternalFunctionDataHelper<TAPE_NAME<TapeTypes>, Data>::CallFunction extFuncPrimal = nullptr){
      ENABLE_CHECK (OptTapeActivity, isActive()){
        pushExternalFunctionHandle(ExternalFunctionDataHelper<TAPE_NAME<TapeTypes>, Data>::createHandle(extFunc, extFuncForward, extFuncPrimal, data, delData));
      }
    }

    /**
     * @brief Adds information about the external functions.
     *
     * Adds the number of external functions.
     *
     * @param[in,out] values  The information is added to the values
     */
    void addExtFuncValues(TapeValues& values) const {
      size_t nExternalFunc = extFuncVector.getDataSize();

      values.addSection("External functions");
      values.addData("Total Number", nExternalFunc);
    }

#undef CHILD_VECTOR_TYPE
#undef CHILD_VECTOR_NAME
#undef VECTOR_TYPE
