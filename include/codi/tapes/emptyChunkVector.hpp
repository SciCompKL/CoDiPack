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

#include "../macros.h"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Implementation for a terminal sequence chunk vector
   *
   * This interface provides the basic implementation for a terminal point
   * in a chain of chunk vectors.
   */
  struct EmptyChunkVector {
    /**
     * @brief Position without any data.
     */
    struct Position {
        /**
         * @brief Returns always false since the position has not data.
         * @param[in] o  The other position.
         * @return Always false.
         */
        bool operator != (const Position& o) {
          CODI_UNUSED(o);
          return false;
        }

        /**
         * @brief Returns always true since the position has not data.
         * @param[in] o  The other position.
         * @return Always true.
         */
        bool operator == (const Position& o) {
          CODI_UNUSED(o);
          return true;
        }
    };

    /**
     * @brief Do nothing.
     *
     * @param[in,out] other  The other empty chunk vector.
     */
    void swap(EmptyChunkVector& other) {
      CODI_UNUSED(other);
      //Do nothing
    }

    /**
     * @brief Do nothing.
     */
    CODI_INLINE void resetHard() {
      // Do nothing
    }

    /**
     * @brief Empty position.
     * @return Empty position
     */
    CODI_INLINE Position getPosition() const {
      return Position();
    }

    /**
     * @brief Empty position.
     * @return Empty position
     */
    CODI_INLINE Position getZeroPosition() const {
      return Position();
    }

    /**
     * @brief Will do nothing.
     */
    CODI_INLINE void reset(const Position& pos) const {
      CODI_UNUSED(pos);
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
     * The position arguments are omitted from the function call.
     *
     * The function call is
     * \code{.cpp}
     * func(<other arguments>, startDataPos, endDataPos, pointerChunkItem1, pointerChunkItem2, etc.);
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

      CODI_UNUSED(start);
      CODI_UNUSED(end);

      function(std::forward<Args>(args)...);
    }

    /**
     * @brief Forward stack evaluation of the tape.
     *
     * The function is called on the accumulated pointers and ranges from all previous vectors.
     * The position arguments are omitted from the function call.
     *
     * The function call is
     * \code{.cpp}
     * func(<other arguments>, startDataPos, endDataPos, pointerChunkItem1, pointerChunkItem2, etc.);
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

      CODI_UNUSED(start);
      CODI_UNUSED(end);

      function(std::forward<Args>(args)...);
    }
  };
}
