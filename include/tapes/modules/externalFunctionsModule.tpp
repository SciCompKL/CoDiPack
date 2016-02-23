/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
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
 * It defines the methods getExtFuncPosition, resetExtFunc, evaluateExtFunc as interface functions for the
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
  // All definitons of the module
  // ----------------------------------------------------------------------

    typedef CHILD_VECTOR_TYPE ExtFuncChildVector;
    typedef typename ExtFuncChildVector::Position ExtFuncChildPosition;

    /** @brief The vector for the external function data. */
    typedef VECTOR_TYPE ExtFuncVector;
    /** @brief The data for the external functions. */
    typedef typename ExtFuncVector::ChunkType ExtFuncChunk;

    typedef typename ExtFuncVector::Position ExtFuncPosition;

    /** @brief The data for the external functions. */
    ExtFuncVector extFuncVector;

  // ----------------------------------------------------------------------
  // Private function of the module
  // ----------------------------------------------------------------------

    /*
     * Function object for the evaluation of the external functions.
     *
     * It stores the last position for the statement vector. With this
     * position it evaluates the statement vector to the position
     * where the external function was added and then calls the
     * external function.
     */
    struct ExtFuncEvaluator {
      ExtFuncChildPosition curInnerPos;

      TAPE_NAME& tape;

      ExtFuncEvaluator(ExtFuncChildPosition curInnerPos, TAPE_NAME& tape) :
        curInnerPos(curInnerPos),
        tape(tape){}

      void operator () (ExternalFunction* extFunc, ExtFuncChildPosition* endInnerPos) {
        // always evaluate the stack to the point of the external function
        tape.evalExtFuncCallback(curInnerPos, *endInnerPos);

        extFunc->evaluate();

        curInnerPos = *endInnerPos;
      }
    };

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
     * @brief Delete the data of the external function.
     * @param extFunction The external function in the vector.
     */
    static void popExternalFunction(ExternalFunction* extFunc, ExtFuncChildPosition* endInnerPos) {
      /* we just need to call the delete function */
      extFunc->deleteData();
    }

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
     * @brief Reset the external function module to the position.
     *
     * The reset will also reset the vector and therefore all nested vectors.
     *
     * @param[in] pos  The position to which the tape is reset.
     */
    void resetExtFunc(const ExtFuncPosition& pos) {
      ExternalFunction* extFunc;
      ExtFuncChildPosition* endInnerPos;
      extFuncVector.forEach(getExtFuncPosition(), pos, popExternalFunction, extFunc, endInnerPos);

      // reset will be done iteratively through the vectors
      extFuncVector.reset(pos);
    }

    /**
     * @brief Evaluate a part of the external function vector.
     *
     * It has to hold start >= end.
     *
     * It calls the evaluation method for the statement vector.
     *
     * @param[in]       start The starting point for the external function vector.
     * @param[in]         end The ending point for the external function vector.
     */
    void evaluateExtFunc(const ExtFuncPosition& start, const ExtFuncPosition &end){
      ExternalFunction* extFunc;
      ExtFuncChildPosition* endInnerPos;
      ExtFuncEvaluator evaluator(start.inner, *this);

      extFuncVector.forEach(start, end, evaluator, extFunc, endInnerPos);

      // Iterate over the reminder also covers the case if there have been no external functions.
      evalExtFuncCallback(evaluator.curInnerPos, end.inner);
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
     * @param[in] extFunc  The external function which is called by the tape.
     * @param[inout] data  The data for the external function. The tape takes ownership over the data.
     * @param[in] delData  The delete function for the data.
     */
    void pushExternalFunctionHandle(ExternalFunction::CallFunction extFunc, void* data, ExternalFunction::DeleteFunction delData){
      ENABLE_CHECK (OptTapeActivity, isActive()){
        pushExternalFunctionHandle(ExternalFunction(extFunc, data, delData));
      }
    }


    /**
     * @brief Add an external function with a specific data type.
     *
     * The data pointer provided to the tape is considered in possession of the tape. The tape will now be responsible to
     * free the data. For this it will use the delete function provided by the user.
     *
     * @param[in] extFunc  The external function which is called by the tape.
     * @param[inout] data  The data for the external function. The tape takes ownership over the data.
     * @param[in] delData  The delete function for the data.
     */
    template<typename Data>
    void pushExternalFunction(typename ExternalFunctionDataHelper<Data>::CallFunction extFunc, Data* data, typename ExternalFunctionDataHelper<Data>::DeleteFunction delData){
      ENABLE_CHECK (OptTapeActivity, isActive()){
        pushExternalFunctionHandle(ExternalFunctionDataHelper<Data>::createHandle(extFunc, data, delData));
      }
    }

    /**
     * @brief Prints statistics about the tape on the screen
     *
     * Prints information such as stored statements/adjoints and memory usage on screen.
     */
    void printExtFuncStatistics() const {
      size_t nExternalFunc = (extFuncVector.getNumChunks()-1)*extFuncVector.getChunkSize()
          +extFuncVector.getChunkUsedData(extFuncVector.getNumChunks()-1);


      std::cout << "-------------------------------------" << std::endl
                << "External functions  "                  << std::endl
                << "-------------------------------------" << std::endl
                << "  Total Number:     " << std::setw(10) << nExternalFunc << std::endl;

    }

#undef CHILD_VECTOR_TYPE
#undef CHILD_VECTOR_NAME
#undef VECTOR_TYPE
