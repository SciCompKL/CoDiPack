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

/*
 * In order to include this file the user has to define the preprocessor macro POSITION_TYPE, INDEX_HANDLER_NAME
 * RESET_FUNCTION_NAME, EVALUATE_FUNCTION_NAME and TAPE_NAME.
 *
 * POSITION_TYPE defines the type of the position structure that is used in the tape.
 * INDEX_HANDLER_NAME defines the name of the index handler.
 * RESET_FUNCTION_NAME is the name of the reset function that is implemented in the including class.
 * EVALUATE_FUNCTION_NAME is the name of the tape evaluation function that is implemented in the including class.
 *
 * All these macros are undefined at the end of the file.
 *
 * TAPE_NAME defines the type name of the tape and is not undefined at the end of the file.
 *
 * The module defines the structures adjoints, adjointSize and active that have to initialized
 * in the including class.
 * The module defines the types Position.
 *
 * It defines the methods initGradientData, destroyGradientData, setGradient, getGradient, gradient, clearAdjoints,
 * reset(Pos), reset(), evaluate(), evaluate(Pos, Pos), setActive, setPassive, isActive, printTapeBaseStatistics
 * from the TapeInterface and ReverseTapeInterface.
 *
 * It defines the methods resizeAdjoints, cleanTapeBase, swapTapeBaseModule as interface functions for the
 * including class.
 */

#ifndef TAPE_NAME
  #error Please define the name of the tape.
#endif

#ifndef POSITION_TYPE
#error Please define the position type.
#endif

#ifndef INDEX_HANDLER_NAME
#error Please define the name of the index handler.
#endif

#ifndef RESET_FUNCTION_NAME
#error Please define the name of the reset function for the tape.
#endif

#ifndef EVALUATE_FUNCTION_NAME
#error Please define the name of the evaluation function for the tape.
#endif

  // ----------------------------------------------------------------------
  // All definitons of the module
  // ----------------------------------------------------------------------

  public:

    /** @brief The position for all the different data vectors. */
    typedef POSITION_TYPE Position;

  private:

    /**
     * @brief The adjoint vector.
     *
     * The size of the adjoint vector is set according to the requested positions.
     * But the positions should not be greater than the current expression counter.
     */
    GradientValue* adjoints;

    /** @brief The current size of the adjoint vector. */
    IndexType adjointsSize;

    /**
     * @brief Determines if statements are recorded or ignored.
     */
    bool active;

  private:

  // ----------------------------------------------------------------------
  // Private function for the communication with the including class
  // ----------------------------------------------------------------------

    /**
     * @brief Helper function: Sets the adjoint vector to a new size.
     *
     * @param[in] size The new size for the adjoint vector.
     */
    void resizeAdjoints(const IndexType& size) {
      IndexType oldSize = adjointsSize;
      adjointsSize = size;

      adjoints = (GradientValue*)realloc(adjoints, sizeof(GradientValue) * (size_t)adjointsSize);

      if(NULL == adjoints) {
        throw std::bad_alloc();
      }

      for(IndexType i = oldSize; i < adjointsSize; ++i) {
        adjoints[i] = GradientValue();
      }
    }

    /**
     * @brief Helper function: Deletes all arrays
     */
    void cleanTapeBase() {
      if(NULL != adjoints) {
        free(adjoints);
        adjoints = NULL;
        adjointsSize = 0;
      }
    }

    /**
     * @brief Swap the data of the tape base module with the data of the other tape base module.
     *
     * @param[in] other  The object with the other tape base module.
     */
    void swapTapeBaseModule(TAPE_NAME<TapeTypes>& other) {
      std::swap(adjoints, other.adjoints);
      std::swap(adjointsSize, other.adjointsSize);
      std::swap(active, other.active);

      // the index handler is not swaped because it is either swaped in the recursive call to of the data vectors
      // or it is handled by the including class
    }

  public:

  // ----------------------------------------------------------------------
  // Public function from the TapeInterface and ReverseTapeInterface
  // ----------------------------------------------------------------------


    /**
     * @brief Set the index to zero.
     * @param[in] value Not used in this implementation.
     * @param[out] index The index of the active type.
     */
    CODI_INLINE void initGradientData(Real& value, IndexType& index) {
      CODI_UNUSED(value);
      index = IndexType();
    }

    /**
     * @brief The free method is called with the index on the index handler.
     * @param[in] value Not used in this implementation.
     * @param[in] index The index of the active type.
     */
    CODI_INLINE void destroyGradientData(Real& value, IndexType& index) {
      CODI_UNUSED(value);

      INDEX_HANDLER_NAME.freeIndex(index);
    }

    /**
     * @brief No check is performed because the gradient values do not exist.
     *
     * @param[in] gradientData  No used in this implementation.
     * @return always true
     */
    CODI_INLINE bool isGradientTotalZero(const GradientData& gradientData) {
      CODI_UNUSED(gradientData);

      return true;
    }


    /**
     * @brief Set the gradient value of the corresponding index.
     *
     * If the index 0 is the inactive indicator and is ignored.
     *
     * @param[in]    index  The index of the active type.
     * @param[in] gradient  The new value for the gradient.
     */
    void setGradient(IndexType& index, const GradientValue& gradient) {
      if(0 != index) {
        this->gradient(index) = gradient;
      }
    }

    /**
     * @brief Check whether the gradient data is zero.
     *
     * @param[in] index The index of the active type.
     * @return False if the gradient data is zero, otherwise returns true.
     */
    bool isActive(const IndexType& index) const{
      return (index != 0);
    }

    /**
     * @brief Get the gradient value of the corresponding index.
     *
     * @param[in] index The index of the active type.
     * @return The gradient value corresponding to the given index.
     */
    CODI_INLINE GradientValue getGradient(const IndexType& index) const {
      if(0 == index || adjointsSize <= index) {
        return GradientValue();
      } else {
        return adjoints[index];
      }
    }

    /**
     * @brief Get a reference to the gradient value of the corresponding index.
     *
     * An index of 0 will raise an codiAssert exception.
     *
     * @param[in] index The index of the active type.
     * @return The reference to the gradient data.
     */
    CODI_INLINE GradientValue& gradient(IndexType& index) {
      codiAssert(0 != index);
      codiAssert(index <= INDEX_HANDLER_NAME.getMaximumGlobalIndex());

      //TODO: Add error when index is bigger than expression count
      if(adjointsSize <= index) {
        resizeAdjoints(INDEX_HANDLER_NAME.getMaximumGlobalIndex() + 1);
      }

      return adjoints[index];
    }

    /**
     * @brief Get a constant reference to the gradient value of the corresponding index.
     *
     * @param[in] index The index of the active type.
     * @return The constant reference to the gradient data.
     */
    CODI_INLINE const GradientValue& gradient(const IndexType& index) const {
      if(adjointsSize <= index) {
        return adjoints[0];
      } else {
        return adjoints[index];
      }
    }

    /**
     * @brief Sets all adjoint/gradients to zero.
     */
    CODI_INLINE void clearAdjoints(){
      if(NULL != adjoints) {
        for(IndexType i = 0; i < adjointsSize; ++i) {
          adjoints[i] = GradientValue();
        }
      }
    }

    /**
     * @brief Reset the tape to the given position.
     *
     * @param[in] pos Reset the state of the tape to the given position.
     */
    CODI_INLINE void reset(const Position& pos) {
      clearAdjoints(getPosition(), pos);

      // reset will be done iteratively through the vectors
      RESET_FUNCTION_NAME(pos);
    }


    /**
     * @brief Reset the tape to its initial state.
     */
    CODI_INLINE void reset() {
      clearAdjoints();

      // reset will be done iteratively through the vectors
      RESET_FUNCTION_NAME(getZeroPosition());
    }

    /**
     * @brief Perform the adjoint evaluation from start to end.
     *
     * It has to hold start >= end.
     *
     * @param[in] start  The starting position for the adjoint evaluation.
     * @param[in]   end  The ending position for the adjoint evaluation.
     */
    void evaluate(const Position& start, const Position& end) {
      if(adjointsSize <= INDEX_HANDLER_NAME.getMaximumGlobalIndex()) {
        resizeAdjoints(INDEX_HANDLER_NAME.getMaximumGlobalIndex() + 1);
      }

      EVALUATE_FUNCTION_NAME(start, end);
    }

    /**
     * @brief Perform the adjoint evaluation from the current position to the initial position.
     */
    void evaluate() {
      evaluate(getPosition(), getZeroPosition());
    }

    /**
     * @brief Start recording.
     */
    CODI_INLINE void setActive(){
      active = true;
    }

    /**
     * @brief Stop recording.
     */
    CODI_INLINE void setPassive(){
      active = false;
    }

    /**
     * @brief Check if the tape is active.
     * @return true if the tape is active.
     */
    CODI_INLINE bool isActive() const {
      return active;
    }

    /**
     * @brief Prints statistics about the adjoint vector.
     *
     * Displays the number of adjoints and the allocated memory. Also
     * displays the information about the index handler.
     *
     * @param[in,out]   out  The information is written to the stream.
     * @param[in]     hLine  The horizontal line that seperates the sections of the output.
     *
     * @tparam Stream The type of the stream.
     */
    template<typename Stream>
    void printTapeBaseStatistics(Stream& out, const std::string hLine) const {

      size_t nAdjoints      = INDEX_HANDLER_NAME.getMaximumGlobalIndex() + 1;
      double memoryAdjoints = (double)nAdjoints * (double)sizeof(GradientValue) * BYTE_TO_MB;

      out << hLine
          << "Adjoint vector\n"
          << hLine
          << "  Number of Adjoints: " << std::setw(10) << nAdjoints << "\n"
          << "  Memory allocated:   " << std::setiosflags(std::ios::fixed)
                                      << std::setprecision(2)
                                      << std::setw(10)
                                      << memoryAdjoints << " MB" << "\n";
      INDEX_HANDLER_NAME.printStatistics(out, hLine);
    }

    /**
     * @brief Return the size of the adjoint vector.
     *
     * @return The size of the adjoint vector.
     */
    size_t getAdjointSize() {
      return INDEX_HANDLER_NAME.getMaximumGlobalIndex();
    }

#undef POSITION_TYPE
#undef INDEX_HANDLER_NAME
#undef RESET_FUNCTION_NAME
#undef EVALUATE_FUNCTION_NAME
