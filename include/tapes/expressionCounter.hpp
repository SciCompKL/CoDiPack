/**
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

#pragma once

namespace codi {

  /**
   * @brief Terminator for a sequence of ChunkVectors, which counts the number of statements.
   *
   * The structure is a helper structure for the ChunkTape. It terminates the ChunkVector sequence
   * and counts the number of expressions which have been recorded on the tape.
   *
   * @tparam IndexType The type of the index for the counting.
   */
  template <typename IndexType>
  struct ExpressionCounter {

    /**
     * @brief The needed position definition for a ChunkVector sequence terminator.
     *
     * Just the integer for the current statement.
     */
    typedef IndexType Position;

    /**
     * @brief The current count for the statements.
     */
    IndexType count;

    /**
     * @brief The needed getPosition method for a ChunkVector sequence terminator.
     *
     * The method returns the current state of the count value.
     * @return The current value of count.
     */
    inline Position getPosition() {
      return count;
    }

    /**
     * @brief The needed reset method for a ChunkVector sequence terminator.
     *
     * The method sets count to the value from the argument.
     *
     * @param pos The new value of count.
     */
    inline void reset(const Position& pos) {
      count = pos;
    }
  };
}
