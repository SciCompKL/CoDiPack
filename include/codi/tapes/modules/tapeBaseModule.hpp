﻿/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *     Max Sagebaum
 *     Tim Albring
 *     Johannes Blühdorn
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
   * The module defines the structures adjoints, adjointSize and active that have to be initialized
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
  // All definitions of the module
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
       * @return  The full tape implementation.
       */
      Tape& cast() {
        return *static_cast<Tape*>(this);
      }

      /**
       * @brief Cast this class to the full.
       *
       * The full type is able to access all functions from the tape and other modules.
       *
       * @return  The full tape implementation.
       */
      const Tape& cast() const {
        return *static_cast<const Tape*>(this);
      }

    public:

      /** @brief True if the index handler of the tape provided linear indices. */
      static const bool LinearIndexHandler = TapeTypes::IndexHandler::IsLinear;

    protected:

      /**
       * @brief Determines if statements are recorded or ignored.
       */
      bool active;

      /**
       * @brief Default constructor
       */
      TapeBaseModule() :
        active(false)
      {}

      ~TapeBaseModule() {
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

      /**
       * @brief Swap the data of the tape base module with the data of the other tape base module.
       *
       * @param[in] other  The object with the other tape base module.
       */
      void swapTapeBaseModule(Tape& other) {
        std::swap(active, other.active);
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
       * @brief Check whether the gradient data is zero.
       *
       * @param[in] index The index of the active type.
       * @return False if the gradient data is zero, otherwise returns true.
       */
      bool isActive(const Index& index) const{
        return (index != 0);
      }

      /**
       * @brief Clear the derivative information from a value.
       *
       * The value is considered afterwards as not dependent on any input variables.
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
      CODI_INLINE void reset(const Position& pos, bool resetAdjoints = true) {
        if (resetAdjoints) {
          cast().clearAdjoints(cast().getPosition(), pos);
        }

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
        cast().resizeAdjointsToIndexSize();

        auto adjoints = cast().getAdjoints();
        cast().lockForUse();
        evaluate(start, end, adjoints);
        cast().unlockAfterUse();
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
        cast().resizeAdjointsToIndexSize();

        auto adjoints = cast().getAdjoints();
        cast().lockForUse();
        evaluateForward(start, end, adjoints);
        cast().unlockAfterUse();
      }

      /**
       * @brief Perform the forward evaluation from the initial position to the current position.
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
       * The data is written in a csv format with semicolons as an separator.
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
       * The data is written in a csv format with semicolons as an separator.
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
  };
}
