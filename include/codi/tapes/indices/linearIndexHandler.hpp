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

#include <vector>

#include "../../configure.h"
#include "../../macros.h"
#include "../../tools/tapeValues.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Provides linear increasing indices.
   *
   * New indices are generated if required. No index is generated twice.
   *
   * The handler will be reset together with the tape.
   *
   * @tparam IndexType  The type for the handled indices.
   */
  template<typename IndexType>
  class LinearIndexHandler {
    public:

      /**
       * @brief The type definition for other tapes who want to access the type.
       */
      typedef IndexType Index;

      /**
       * @brief The required position definition for a ChunkVector sequence terminator.
       *
       * Just the integer for the current statement.
       */
      typedef Index Position;

      /**
       * @brief Indicates if the index handler privides linear increasing indices.
       *
       * true for this index manager.
       */
      static const bool IsLinear = true;

    private:

      /**
       * @brief The amount of indices that is not given to statements.
       *
       * The indices 1 to zeroState can be used by the tape.
       */
      Index zeroState;

      /**
       * @brief The current count for the indices or statements.
       */
      Index count;
    public:


      /**
       * @brief Create a handler that starts with the index one.
       */
      LinearIndexHandler(Index zeroState) :
        zeroState(zeroState),
        count(zeroState) {}

      /**
       * @brief Swap the contents of this linear index handler with the contents of the other
       *        linear index handler.
       *
       * On standard containers the default std::swap method is used.
       *
       * @param[in,out] other  The other linear index handler.
       */
      void swap(LinearIndexHandler<Index>& other) {
        std::swap(zeroState, other.zeroState);
        std::swap(count, other.count);
      }

      /**
       * @brief Free the index that is given to the method.
       *
       * Just zeros the index
       *
       * @param[in,out] index  The index that is freed. It is set to zero in the method.
       */
      CODI_INLINE void freeIndex(Index& index) const {
        index = 0;
      }

      /**
       * @brief Generate a new index.
       *
       * The indices are linear increasing.
       *
       * @return The new index that can be used.
       */
      CODI_INLINE Index createIndex() {
        return ++count;
      }

      /**
       * @brief No check needs to be performed. A new index is always generated.
       *
       * @param[out] index Will be set to a new index.
       */
      CODI_INLINE void assignIndex(Index& index) {
        index = this->createIndex();
      }

      /**
       * @brief Get the maximum global index.
       *
       * @return The maximum index that was used during the lifetime of this index handler.
       */
      CODI_INLINE Index getMaximumGlobalIndex() const {
        return count;
      }

      /**
       * @brief Get the current maximum index.
       *
       * @return The current maximum index that is in use.
       */
      CODI_INLINE Index getCurrentIndex() const {
        return count;
      }

      /**
       * @brief The needed getPosition method for a ChunkVector sequence terminator.
       *
       * The method returns the current state of the count value.
       * @return The current value of count.
       */
      CODI_INLINE Position getPosition() const {
        return count;
      }

      /**
       * @brief The needed getZeroPosition method for a ChunkVector sequence terminator.
       *
       * The method returns the zero state of terminator.
       * @return The zero state of the terminator.
       */
      CODI_INLINE Position getZeroPosition() const {
        return zeroState;
      }

      /**
       * @brief The needed reset method for a ChunkVector sequence terminator.
       *
       * The method sets count to the value from the argument.
       *
       * @param pos The new value of count.
       */
      CODI_INLINE void reset(const Position& pos) {
        codiAssert(pos >= zeroState);

        count = pos;
      }

      /**
       * @brief Reset the position to the zero state.
       */
      CODI_INLINE void reset() {
        count = zeroState;
      }

      /**
       * @brief Resets the linear index handler to the initial state.
       */
      CODI_INLINE void resetHard() {
        count = zeroState;
      }

      /**
       * @ brief There are no statistics for this handler.
       *
       * @param[in,out] values  The values where the information is added to.
       */
      void addValues(TapeValues& values) const {
        CODI_UNUSED(values);

        // Do nothing
      }

    /**
     * @brief There are no chunks, that need to be iterated.
     *
     * @param  function  The function called for each chunk.
     * @param recursive  If also the chunks of the nested vectors should be iterated.
     * @param      args  The pointers are used as the arguments for the function.
     *
     * @tparam  Args  The data types for the arguments of the function.
     */
      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachChunkForward(FunctionObject& function, bool recursive, Args &... args) {
        CODI_UNUSED(function);
        CODI_UNUSED(recursive);
        CODI_UNUSED_VAR(args...);

        // Do nothing
      }

      /**
       * @brief Reverse stack evaluation of the tape.
       *
       * The function is called on the accumulated pointers and ranges from all previous vectors.
       *
       * The function call is
       * \code{.cpp}
       * func(start, end, <other arguments>, startDataPos, endDataPos, pointerChunkItem1, pointerChunkItem2, etc.);
       * \endcode
       *
       * It has to hold start >= end.
       *
       * @param    start  The start point for the stack interpretation.
       * @param      end  The end point for the stack interpretation.
       * @param function  The function called for each valid range.
       * @param     args  Pointers and ranges from other chunks vectors and additional arguments for the
       *                  function.
       *
       * @tparam  Args  The data types for the arguments.
       */
      template<typename Function, typename ... Args>
      CODI_INLINE void evaluateReverse(const Position& start, const Position& end,const Function& function,
                                       Args&&... args) {

        function(start, end, std::forward<Args>(args)...);
      }

      /**
       * @brief Forward stack evaluation of the tape.
       *
       * The function is called on the accumulated pointers and ranges from all previous vectors.
       *
       * The function call is
       * \code{.cpp}
       * func(start, end, <other arguments>, startDataPos, endDataPos, pointerChunkItem1, pointerChunkItem2, etc.);
       * \endcode
       *
       * It has to hold start <= end.
       *
       * @param    start  The start point for the stack interpretation.
       * @param      end  The end point for the stack interpretation.
       * @param function  The function called for each valid range.
       * @param     args  Pointers and ranges from other chunks vectors and additional arguments for the
       *                  function.
       *
       * @tparam  Args  The data types for the arguments.
       */
      template<typename Function, typename ... Args>
      CODI_INLINE void evaluateForward(const Position& start, const Position& end,const Function& function,
                                       Args&&... args) {

        function(start, end, std::forward<Args>(args)...);
      }
  };
}
