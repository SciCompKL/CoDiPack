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

#ifndef POSITION_TYPE
#error Please define the position type.
#endif

#ifndef INDEX_HANDLER_TYPE
#error Please define the type for the index handler.
#endif

#ifndef RESET_FUNCTION_NAME
#error Please define the name of the reset function for the tape.
#endif

#ifndef EVALUATE_FUNCTION_NAME
#error Please define the name of the evaluation function for the tape.
#endif

  public:

    /** @brief The position for all the different data vectors. */
    typedef POSITION_TYPE Position;

  private:

    INDEX_HANDLER_TYPE indexHandler;

    /**
     * @brief The adjoint vector.
     *
     * The size of the adjoint vector is set according to the requested positions.
     * But the positions should not be greater than the current expression counter.
     */
    Real* adjoints;
    /* @brief The current size of the adjoint vector. */
    IndexType adjointsSize;

    /**
     * @brief Determines if statements are recorded or ignored.
     */
    bool active;

    /**
     * @brief Helper function: Sets the adjoint vector to a new size.
     *
     * @param[in] size The new size for the adjoint vector.
     */
    void resizeAdjoints(const IndexType& size) {
      IndexType oldSize = adjointsSize;
      adjointsSize = size;

      adjoints = (Real*)realloc(adjoints, sizeof(Real) * (size_t)adjointsSize);

      for(IndexType i = oldSize; i < adjointsSize; ++i) {
        adjoints[i] = Real();
      }
    }

  public:

    /**
     * @brief Set the index to zero.
     * @param[in] value Not used in this implementation.
     * @param[out] index The index of the active type.
     */
    inline void initGradientData(Real& value, IndexType& index) {
      CODI_UNUSED(value);
      index = IndexType();
    }

    /**
     * @brief Does nothing.
     * @param[in] value Not used in this implementation.
     * @param[in] index Not used in this implementation.
     */
    inline void destroyGradientData(Real& value, IndexType& index) {
      CODI_UNUSED(value);

      indexHandler.freeIndex(index);
    }


    /**
     * @brief Set the gradient value of the corresponding index.
     *
     * If the index 0 is the inactive indicator and is ignored.
     *
     * @param[in]    index  The index of the active type.
     * @param[in] gradient  The new value for the gradient.
     */
    void setGradient(IndexType& index, const Real& gradient) {
      if(0 != index) {
        this->gradient(index) = gradient;
      }
    }

    /**
     * @brief Get the gradient value of the corresponding index.
     *
     * @param[in] index The index of the active type.
     * @return The gradient value corresponding to the given index.
     */
    inline Real getGradient(const IndexType& index) const {
      if(adjointsSize <= index) {
        return Real();
      } else {
        return adjoints[index];
      }
    }

    /**
     * @brief Get a reference to the gradient value of the corresponding index.
     *
     * An index of 0 will raise an assert exception.
     *
     * @param[in] index The index of the active type.
     * @return The reference to the gradient data.
     */
    inline Real& gradient(IndexType& index) {
      assert(0 != index);
      assert(index <= indexHandler.getMaximumGlobalIndex());

      //TODO: Add error when index is bigger than expression count
      if(adjointsSize <= index) {
        resizeAdjoints(indexHandler.getMaximumGlobalIndex() + 1);
      }

      return adjoints[index];
    }

    /**
     * @brief Get the current position of the tape.
     *
     * The position can be used to reset the tape to that position or to
     * evaluate only parts of the tape.
     * @return The current position of the tape.
     */
    inline Position getPosition() const {
      return getExtFuncPosition();
    }

    /**
     * @brief Sets all adjoint/gradients to zero.
     */
    inline void clearAdjoints(){
      for(IndexType i = 0; i <= indexHandler.getMaximumGlobalIndex(); ++i) {
        adjoints[i] = 0.0;
      }
    }

    /**
     * @brief Reset the tape to the given position.
     *
     * @param[in] pos Reset the state of the tape to the given position.
     */
    inline void reset(const Position& pos) {
      clearAdjoints(getPosition(), pos);

      // reset will be done iteratively through the vectors
      RESET_FUNCTION_NAME(pos);
    }


    /**
     * @brief Reset the tape to its initial state.
     */
    inline void reset() {
      reset(Position());
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
      if(adjointsSize <= indexHandler.getMaximumGlobalIndex()) {
        resizeAdjoints(indexHandler.getMaximumGlobalIndex() + 1);
      }

      EVALUATE_FUNCTION_NAME(start, end);
    }

    /**
     * @brief Perform the adjoint evaluation from the current position to the initial position.
     */
    void evaluate() {
      evaluate(getPosition(), Position());
    }

    /**
     * @brief Start recording.
     */
    inline void setActive(){
      active = true;
    }

    /**
     * @brief Stop recording.
     */
    inline void setPassive(){
      active = false;
    }

    /**
     * @brief Check if the tape is active.
     * @return true if the tape is active.
     */
    inline bool isActive(){
      return active;
    }

    void printTapeBaseStatistics(){

      size_t nAdjoints      = indexHandler.getMaximumGlobalIndex() + 1;
      double memoryAdjoints = (double)nAdjoints * (double)sizeof(Real) * BYTE_TO_MB;

      std::cout << std::endl
                << "-------------------------------------" << std::endl
                << "Adjoint vector"                        << std::endl
                << "-------------------------------------" << std::endl
                << "  Number of Adjoints: " << std::setw(10) << nAdjoints << std::endl
                << "  Memory allocated:   " << std::setiosflags(std::ios::fixed)
                                            << std::setprecision(2)
                                            << std::setw(10)
                                            << memoryAdjoints << " MB" << std::endl;
     indexHandler.printStatistics();
    }

#undef POSITION_TYPE
#undef INDEX_HANDLER_TYPE
#undef RESET_FUNCTION_NAME
#undef EVALUATE_FUNCTION_NAME
