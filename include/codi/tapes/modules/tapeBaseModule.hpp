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

#include "../chunk.hpp"
#include "../reverseTapeInterface.hpp"
#include "../../configure.h"
#include "../../tapeTypes.hpp"
#include "../../tools/tapeValues.hpp"
#include "../../typeFunctions.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * The module defines the structures adjoints, adjointSize and active that have to initialized
   * in the including class.
   * The module defines the types Position.
   *
   * It defines the methods initGradientData, destroyGradientData, getPosition, getZeroPosition, setGradient, getGradient, gradient, clearAdjoints,
   * reset(Pos), reset(), evaluate(), evaluate(Pos, Pos), evaluateForward(), evaluateForward(Pos, Pos), setActive,
   * setPassive, isActive, print Statistics from the TapeInterface and ReverseTapeInterface.
   *
   * It defines the methods resizeAdjoints, resizeAdjointsToIndexSize, cleanTapeBase, swapTapeBaseModule as interface functions for the
   * including class.
   *
   * @tparam    TapeTypes  All the types for the tape. Including the calculation type and the vector types.
   * @tparam         Tape  The full tape implementation
   */
  template<typename TapeTypes, typename Tape>
  struct TapeBaseModule : public virtual ReverseTapeInterface<typename TapeTypes::Real, typename TapeTypes::Index, typename TapeTypes::GradientValue, Tape, typename TapeTypes::Position > {

    private:
  // ----------------------------------------------------------------------
  // All definitons of the module
  // ----------------------------------------------------------------------

      CODI_INLINE_REVERSE_TAPE_TYPES(TapeTypes::BaseTypes)

      /** @brief The global position for the tape */
      typedef typename TapeTypes::Position Position;

      /** @brief Forward the GradientData from the tape types. */
      typedef typename TapeTypes::GradientData GradientData;

      /**
       * @brief Cast this class to the full.
       *
       * The full type is able to access all functions from the tape and other modules.
       *
       * @return  The full tape implemenation.
       */
      Tape& cast() {
        return *static_cast<Tape*>(this);
      }

      /**
       * @brief Cast this class to the full.
       *
       * The full type is able to access all functions from the tape and other modules.
       *
       * @return  The full tape implemenation.
       */
      const Tape& cast() const {
        return *static_cast<const Tape*>(this);
      }

    public:

      /** @brief True if the index handler of the tape provied linear indices. */
      static const bool LinearIndexHandler = TapeTypes::IndexHandler::IsLinear;

    protected:

      /**
       * @brief The adjoint vector.
       *
       * The size of the adjoint vector is set according to the requested positions.
       * But the positions should not be greater than the current expression counter.
       */
      GradientValue* adjoints;

      /** @brief The current size of the adjoint vector. */
      Index adjointsSize;

      /**
       * @brief Determines if statements are recorded or ignored.
       */
      bool active;

      /**
       * @brief Default constructor
       */
      TapeBaseModule() :
        adjoints(NULL),
        adjointsSize(0),
        active(false)
      {}

      ~TapeBaseModule() {
        cleanTapeBase();
      }

    protected:

      /**
       * @brief Initialize the TapeBaseModule.
       *
       * Called after all members of the tape have been initialized.
       */
      void initTapeBaseModule() {
        // Nothing to do
      }

    protected:

    // ----------------------------------------------------------------------
    // Protected function for the communication with the including class
    // ----------------------------------------------------------------------

      /**
      * @brief Adds information about adjoint vector.
      *
      * Adds the number of adjoint vector entries and the size of the adjoint vector.
      *
      * @param[in,out] values  The information is added to the values
      */
      void addTapeBaseValues(TapeValues& values) const {

        size_t nAdjoints      = cast().indexHandler.getMaximumGlobalIndex() + 1;
        double memoryAdjoints = (double)nAdjoints * (double)sizeof(GradientValue) * BYTE_TO_MB;

        values.addSection("Adjoint vector");
        values.addData("Number of adjoints", nAdjoints);
        values.addData("Memory allocated", memoryAdjoints, true, true);

        cast().indexHandler.addValues(values);
      }

      /**
       * @brief Helper function: Sets the adjoint vector to a new size.
       *
       * @param[in] size The new size for the adjoint vector.
       */
      void resizeAdjoints(const Index& size) {
        Index oldSize = adjointsSize;
        adjointsSize = size;

        for(Index i = adjointsSize; i < oldSize; ++i) {
          adjoints[i].~GradientValue();
        }

        adjoints = (GradientValue*)realloc((void*)adjoints, sizeof(GradientValue) * (size_t)adjointsSize);

        if(NULL == adjoints) {
          throw std::bad_alloc();
        }

        for(Index i = oldSize; i < adjointsSize; ++i) {
          new (adjoints + i) GradientValue();
        }
      }

      /**
       * @brief Resize the adjoint vector such that it fits the number of indices.
       */
      void resizeAdjointsToIndexSize() {
        if(adjointsSize <= cast().indexHandler.getMaximumGlobalIndex()) {
          resizeAdjoints(cast().indexHandler.getMaximumGlobalIndex() + 1);
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
      void swapTapeBaseModule(Tape& other) {
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
      CODI_INLINE void initGradientData(Real& value, Index& index) {
        CODI_UNUSED(value);
        index = Index();
      }

      /**
       * @brief The free method is called with the index on the index handler.
       * @param[in] value Not used in this implementation.
       * @param[in] index The index of the active type.
       */
      CODI_INLINE void destroyGradientData(Real& value, Index& index) {
        CODI_UNUSED(value);

        cast().indexHandler.freeIndex(index);
      }

      /**
       * @brief Get the current position of the tape.
       *
       * The position can be used to reset the tape to that position or to
       * evaluate only parts of the tape.
       * @return The current position of the tape.
       */
      CODI_INLINE Position getPosition() const {
        return cast().getRootVector().getPosition();
      }

      /**
       * @brief Get the initial position of the tape.
       *
       * The position can be used to reset the tape to that position or to
       * evaluate only parts of the tape.
       * @return The initial position of the tape.
       */
      CODI_INLINE Position getZeroPosition() const {
        return cast().getRootVector().getZeroPosition();
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
      void setGradient(Index& index, const GradientValue& gradient) {
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
      bool isActive(const Index& index) const{
        return (index != 0);
      }

      /**
       * @brief Get the gradient value of the corresponding index.
       *
       * @param[in] index The index of the active type.
       * @return The gradient value corresponding to the given index.
       */
      CODI_INLINE GradientValue getGradient(const Index& index) const {
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
      CODI_INLINE GradientValue& gradient(Index& index) {
        codiAssert(0 != index);
        codiAssert(index <= cast().indexHandler.getMaximumGlobalIndex());

        //TODO: Add error when index is bigger than expression count
        if(adjointsSize <= index) {
          resizeAdjoints(cast().indexHandler.getMaximumGlobalIndex() + 1);
        }

        return adjoints[index];
      }

      /**
       * @brief Get a constant reference to the gradient value of the corresponding index.
       *
       * @param[in] index The index of the active type.
       * @return The constant reference to the gradient data.
       */
      CODI_INLINE const GradientValue& gradient(const Index& index) const {
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
          for(Index i = 0; i < adjointsSize; ++i) {
            adjoints[i] = GradientValue();
          }
        }
      }

      /**
       * @brief Clear the derivative information from a value.
       *
       * The value is considered afterswards as not dependent on any input variables.
       *
       * @param[in,out] value The cleared variable.
       */
      void deactivateValue(ActiveReal<Tape>& value) {
        cast().indexHandler.freeIndex(value.getGradientData());
      }

      /**
       * @brief Reset the tape to the given position.
       *
       * @param[in] pos Reset the state of the tape to the given position.
       */
      CODI_INLINE void reset(const Position& pos) {
        cast().clearAdjoints(cast().getPosition(), pos);

        // reset will be done iteratively through the vectors
        cast().resetInternal(pos);
      }


      /**
       * @brief Reset the tape to its initial state.
       */
      CODI_INLINE void reset(bool resetAdjoints = true) {
        if(resetAdjoints) {
          cast().clearAdjoints();
        }

        cast().indexHandler.reset();

        // reset will be done iteratively through the vectors
        cast().resetInternal(cast().getZeroPosition());
      }

      /**
       * @brief Perform the adjoint evaluation from start to end with a custom adjoint vector.
       *
       * It has to hold start >= end.
       *
       * @param[in]           start  The starting position for the adjoint evaluation.
       * @param[in]             end  The ending position for the adjoint evaluation.
       * @param[in,out] adjointData  The vector for the adjoint evaluation. It has to have the size of getAdjointSize() + 1.
       *
       * @tparam AdjointData  The type needs to provide an add, multiply and comparison operation.
       */
      template<typename AdjointData>
      CODI_NO_INLINE void evaluate(const Position& start, const Position& end, AdjointData* adjointData) {
        cast().evaluateInternal(start, end, adjointData);
      }

      /**
       * @brief Perform the adjoint evaluation from start to end.
       *
       * It has to hold start >= end.
       *
       * @param[in] start  The starting position for the adjoint evaluation.
       * @param[in]   end  The ending position for the adjoint evaluation.
       */
      CODI_NO_INLINE void evaluate(const Position& start, const Position& end) {
        resizeAdjointsToIndexSize();

        evaluate(start, end, adjoints);
      }

      /**
       * @brief Perform the adjoint evaluation from the current position to the initial position.
       */
      void evaluate() {
        evaluate(cast().getPosition(), cast().getZeroPosition());
      }

      /**
       * @brief Perform the forward evaluation from start to end with a custom adjoint vector.
       *
       * It has to hold start <= end.
       *
       * @param[in]           start  The starting position for the forward evaluation.
       * @param[in]             end  The ending position for the forward evaluation.
       * @param[in,out] adjointData  The vector for the adjoint evaluation. It has to have the size of getAdjointSize() + 1.
       *
       * @tparam AdjointData  The type needs to provide an add, multiply and comparison operation.
       */
      template<typename AdjointData>
      CODI_NO_INLINE void evaluateForward(const Position& start, const Position& end, AdjointData* adjointData) {
        cast().evaluateForwardInternal(start, end, adjointData);
      }

      /**
       * @brief Perform the forward evaluation from start to end.
       *
       * It has to hold start <= end.
       *
       * @param[in] start  The starting position for the forward evaluation.
       * @param[in]   end  The ending position for the forward evaluation.
       */
      CODI_NO_INLINE void evaluateForward(const Position& start, const Position& end) {
        resizeAdjointsToIndexSize();

        evaluateForward(start, end, adjoints);
      }

      /**
       * @brief Perform the foward evaluation from the initial position to the current position.
       */
      void evaluateForward() {
        evaluateForward(cast().getZeroPosition(), cast().getPosition());
      }

      /**
       * @brief Perform the primal evaluation from start to end.
       *
       * It has to hold start <= end.
       *
       * @param[in] start  The starting position for the forward evaluation.
       * @param[in]   end  The ending position for the forward evaluation.
       */
      CODI_NO_INLINE void evaluatePrimal(const Position& start, const Position& end) {
        cast().evaluatePrimalInternal(start, end);
      }

      /**
       * @brief Perform the primal evaluation from the initial position to the current position.
       */
      void evaluatePrimal() {
        evaluatePrimal(cast().getZeroPosition(), cast().getPosition());
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
       * @brief The default passive index from the tape.
       *
       * @return Zero in the defined index type.
       */
      Index getPassiveIndex() const {
        return Index(0);
      }

      /**
       * @brief The default invalid index from the tape.
       *
       * @return Minus one in the defined index type.
       */
      Index getInvalidIndex() const {
        return Index(-1);
      }

      /**
       * @brief Prints statistics about the tape on the screen or into a stream
       *
       * Prints information such as stored statements/adjoints and memory usage on screen or into
       * the stream when an argument is provided.
       *
       * @param[in,out] out  The information is written to the stream.
       *
       * @tparam Stream The type of the stream.
       */
      template<typename Stream = std::ostream>
      void printStatistics(Stream& out = std::cout) const {

        TapeValues values = cast().getTapeValues();

        values.formatDefault(out);
      }

      /**
       * @brief Print the table header for the tape information to the stream.
       *
       * The data is written in a csv format with semicolons as an seperator.
       *
       * @param[in,out] out  The stream which is used for the printing of the table header.
       * @tparam Stream  The type of the stream.
       */
      template<typename Stream = std::ostream>
      void printTableHeader(Stream& out = std::cout) const {

        TapeValues values = cast().getTapeValues();

        values.formatHeader(out);
      }

      /**
       * @brief Print the table data of the current tape state to the stream.
       *
       * The data is written in a csv format with semicolons as an seperator.
       *
       * @param[in,out] out  The stream which is used for the printing of the table data.
       * @tparam Stream  The type of the stream.
       */
      template<typename Stream = std::ostream>
      void printTableRow(Stream& out = std::cout) const {

        TapeValues values = cast().getTapeValues();

        values.formatRow(out);
      }

      /**
       * @brief Return the size of the adjoint vector.
       *
       * @return The size of the adjoint vector.
       */
      size_t getAdjointSize() const {
        return cast().indexHandler.getMaximumGlobalIndex();
      }

      /**
       * @brief Clear the adjoint vector and delete it.
       */
      void deleteAdjointVector() {
        cleanTapeBase();
      }
  };
}
