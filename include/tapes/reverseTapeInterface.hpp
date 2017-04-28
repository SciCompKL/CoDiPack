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

#include "../activeReal.hpp"
#include "externalFunctions.hpp"
#include "tapeInterface.hpp"


/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {
  /**
   * @brief Interface common to all reverse tapes.
   *
   * The extended interface each reverse tape has to implement. It defines functions which are used to activate
   * or deactivate the recording. The user can add external functions and register the input variable and output
   * variables for the computation.
   *
   * @tparam               Real  Floating point type of the gradients.
   * @tparam   GradientDataType  The data the tape uses to identify each active variable
   *                               and where the tape can store information about the
   *                               gradient.
   * @tparam   GradientValueType The value type that is used for the gradient calculation.
   * @tparam TapeImplementation  The implementing tape of the interface. It is needed to define the active type
   *                               for the registration of the variables.
   * @tparam           Position  Position used by the implementing tape.
   *
   */
  template <typename Real, typename GradientDataType, typename GradientValueType, typename TapeImplementation, typename Position>
  class ReverseTapeInterface : public TapeInterface<Real, GradientDataType, GradientValueType> {
  public:

    /**
     * @brief Evaluate the tape from start to end.
     *
     * The function performs the reverse evaluation of the recorded tape from
     * the start position to the end position.
     *
     * It has to hold start >= end.
     *
     * @param[in] start The starting position for the reverse evaluation.
     * @param[in]   end The ending position for the reverse evaluation.
     */
    virtual void evaluate(const Position& start, const Position& end) = 0;

    /**
     * @brief Evaluate the tape from the current position to the beginning.
     */
    virtual void evaluate() = 0;

    /**
     * @brief Declare a variable as an input variable.
     *
     * @param[in,out] value The input variable.
     */
    virtual void registerInput(ActiveReal<TapeImplementation>& value) = 0;

    /**
     * @brief Declare a variable as an output variable.
     *
     * @param[in,out] value The output variable.
     */
    virtual void registerOutput(ActiveReal<TapeImplementation>& value) = 0;

    /**
    * @brief Set the tape to active.
    *
    *  While active each operation involving active variables is stored on the tape.
    *
    */
    virtual void setActive() = 0;


    /**
    * @brief Set the tape to passive.
    *
    *  While passive no operation involving active variables is stored on the tape.
    *
    */
    virtual void setPassive() = 0;

    using TapeInterface<Real, GradientDataType, GradientValueType>::isActive; // enable both is active methods in this interface

    /**
    * @brief Get the current status if of tape.
    *
    * @return The current state. If true the tape is active.
    */
    virtual bool isActive() const = 0;

    /**
    * @brief Clears the currently stored adjoints.
    *
    * @return Sets the currently stored adjoints to zero, thereby enabling a reevaluation of the tape.
    */
    virtual void clearAdjoints() = 0;

    /**
     * @brief Reset the tape to the given position.
     *
     * The reset will clear everything the tape has recorded after the given position.
     *
     * @param[in] pos The position to which the tape is reset.
     */
    virtual void reset(const Position& pos) = 0;

    /**
     * @brief Completely reset the tape.
     *
     * The reset will clear everything the tape has recorded.
     */
    virtual void reset() = 0;

    /**
     * @brief Get the current position of the tape.
     *
     * The position can be used to reset the tape to that position,
     * as a start point or end point for the evaluation.
     *
     * @return The current position of the tape.
     */
    virtual Position getPosition() const = 0;

    /**
     * @brief Get the initial position of the tape.
     *
     * The position can be used to reset the tape to that position,
     * as a start point or end point for the evaluation.
     *
     * @return The initial position of the tape.
     */
    virtual Position getZeroPosition() const = 0;

    /**
     * @brief Print some statistics about the currently stored information.
     *
     * @param[in,out] out  The information is written to the stream.
     *
     * @tparam Stream The type of the stream.
     */
    template<typename Stream = std::ostream>
    void printStatistics(Stream& out = std::cout) const;

    /**
     * @brief Add a external function to the tape.
     *
     * The external function is called during the reverse evaluation of the tape. It can be used to
     * give special treatment to code sections which have simpler reverse implementation than the
     * AD tool.
     *
     * @param[in]       extFunc The function which is called during the reverse evaluation of the tape.
     * @param[in]    checkpoint The data argument for the function. The tape takes procession of the data and will delete it.
     * @param[in] delCheckpoint The delete function for the data.
     */
    virtual void pushExternalFunctionHandle(ExternalFunction::CallFunction extFunc, void* checkpoint, ExternalFunction::DeleteFunction delCheckpoint) = 0;

    /**
     * @brief Add a external function to the tape.
     *
     * The external function is called during the reverse evaluation of the tape. It can be used to
     * give special treatment to code sections which have simpler reverse implementation than the
     * AD tool.
     *
     * @param[in]       extFunc The function which is called during the reverse evaluation of the tape.
     * @param[in]    checkpoint The data argument for the function. The tape takes procession of the data and will delete it.
     * @param[in] delCheckpoint The delete function for the data.
     *
     * @tparam Data The data type for the data.
     */
    template<typename Data>
    void pushExternalFunction(
        typename ExternalFunctionDataHelper<TapeImplementation, Data>::CallFunction extFunc,
        Data* checkpoint,
        typename ExternalFunctionDataHelper<TapeImplementation, Data>::DeleteFunction delCheckpoint);

  };
}
