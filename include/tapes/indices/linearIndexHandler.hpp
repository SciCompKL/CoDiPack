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

#pragma once

#include <vector>

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
   * @tparam Index  The type for the handled indices.
   */
  template<typename Index>
  class LinearIndexHandler {
    public:

      /**
       * @brief The type definition for other tapes who want to access the type.
       */
      typedef Index IndexType;

      /**
       * @brief The required position definition for a ChunkVector sequence terminator.
       *
       * Just the integer for the current statement.
       */
      typedef Index Position;

    private:

      /**
       * @brief The current count for the indices or statements.
       */
      Index count;

    public:


      /**
       * @brief Create a handler that starts with the index one.
       */
      LinearIndexHandler() :
        count() {}

      /**
       * @brief Free the index that is given to the method.
       *
       * Just zeros the index
       *
       * @param[inout] index  The index that is freed. It is set to zero in the method.
       */
      inline void freeIndex(Index& index) const {
        index = 0;
      }

      /**
       * @brief Generate a new index.
       *
       * The indices are linear increasing.
       *
       * @return The new index that can be used.
       */
      inline Index createIndex() {
        return ++count;
      }

      /**
       * @brief No check needs to be performed. A new index is always generated.
       *
       * @param[out] index Will be set to a new index.
       */
      inline void checkIndex(Index& index) {
        index = this->createIndex();
      }

      /**
       * @brief Get the maximum global index.
       *
       * @return The maximum index that was used during the lifetime of this index handler.
       */
      inline Index getMaximumGlobalIndex() const {
        return count;
      }

      /**
       * @brief Get the current maximum index.
       *
       * @return The current maximum index that is in use.
       */
      inline Index getCurrentIndex() const {
        return count;
      }

      /**
       * @brief The needed getPosition method for a ChunkVector sequence terminator.
       *
       * The method returns the current state of the count value.
       * @return The current value of count.
       */
      inline Position getPosition() const {
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

      /**
       * @ brief There are no statistics for this handler.
       * @param[in,out] out  The information is written to the stream.
       *
       * @tparam Stream The type of the stream.
       */
      template<typename Stream>
      void printStatistics(Stream& out, const std::string hLine) const {
        CODI_UNUSED(out);
        // Do nothing
      }
  };
}
