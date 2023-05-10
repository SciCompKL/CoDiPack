/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <cstdlib>
#include <iostream>

#include "../../config.h"
#include "../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Empty Position with no nested data.
  struct EmptyPosition {
    public:

      /// Always false
      bool operator!=(EmptyPosition const& o) const {
        CODI_UNUSED(o);
        return false;
      }

      /// Always true
      bool operator==(EmptyPosition const& o) const {
        CODI_UNUSED(o);
        return true;
      }

      /// Always false
      bool operator<(EmptyPosition const& o) const {
        CODI_UNUSED(o);
        return false;
      }

      /// Always true
      bool operator<=(EmptyPosition const& o) const {
        CODI_UNUSED(o);
        return true;
      }

      /// Always false
      bool operator>(EmptyPosition const& o) const {
        CODI_UNUSED(o);
        return false;
      }

      /// Always true
      bool operator>=(EmptyPosition const& o) const {
        CODI_UNUSED(o);
        return true;
      }

      /// Stream output
      friend std::ostream& operator<<(std::ostream& stream, EmptyPosition const& pos) {
        CODI_UNUSED(pos);
        stream << "[]";
        return stream;
      }
  };

  /**
   * @brief Position with one index for e.g. array access.
   *
   * Used for data that is allocated en bloc, e.g. BlockData.
   *
   * @tparam T_NestedPosition  Position implementation
   */
  template<typename T_NestedPosition>
  struct ArrayPosition {
    public:

      using NestedPosition = CODI_DD(T_NestedPosition, EmptyPosition);  ///< See ArrayPosition

      size_t data;  ///< Array position index.

      NestedPosition inner;  ///< Position of nested data

      /// Constructor
      ArrayPosition() : data(0), inner() {}

      /// Constructor
      ArrayPosition(size_t const& data, NestedPosition const& inner) : data(data), inner(inner) {}

      /// Operator != also compares with nested data
      bool operator!=(ArrayPosition const& o) const {
        return data != o.data || this->inner != o.inner;
      }

      /// Operator == also compares with nested data
      bool operator==(ArrayPosition const& o) const {
        return data == o.data && this->inner == o.inner;
      }

      /// Operator < also compares with nested data
      bool operator<(ArrayPosition const& o) const {
        return data < o.data || (data == o.data && inner < o.inner);
      }

      /// Operator <= also compares with nested data
      bool operator<=(ArrayPosition const& o) const {
        return data < o.data || (data == o.data && inner <= o.inner);
      }

      /// Operator > also compares with nested data
      bool operator>(ArrayPosition const& o) const {
        return o < *this;
      }

      /// Operator >= also compares with nested data
      bool operator>=(ArrayPosition const& o) const {
        return o <= *this;
      }

      /// Stream output
      friend std::ostream& operator<<(std::ostream& stream, ArrayPosition const& pos) {
        stream << "[" << pos.data << ", " << pos.inner << "]";
        return stream;
      }
  };

  /**
   * @brief Position with two indices for e.g. chunked data access.
   *
   * Used for data that is allocated with multiple chunks, e.g. ChunkData.
   * `chunk` is the major index and identifies the chunk and `Base::data` is the secondary index which refers to the
   * position in the chunk.
   *
   * For `p1 < p2` it is enough that `p1.chunk < p2.chunk`, only in the equality case `Base::data` needs to be checked.
   *
   * @tparam T_NestedPosition  Position implementation
   */
  template<typename T_NestedPosition>
  struct ChunkPosition : public ArrayPosition<T_NestedPosition> {
    public:

      using NestedPosition = CODI_DD(T_NestedPosition, EmptyPosition);  ///< See ChunkPosition
      using Base = ArrayPosition<NestedPosition>;                       ///< Base abbreviation

      size_t chunk;  ///< Chunk position index

      /// Constructor
      ChunkPosition() : Base(), chunk(0) {}

      /// Constructor
      ChunkPosition(size_t const& chunk, size_t const& data, NestedPosition const& inner)
          : Base(data, inner), chunk(chunk) {}

      /// Operator != also compares with nested data
      bool operator!=(ChunkPosition const& o) const {
        return chunk != o.chunk || Base::operator!=(static_cast<Base const&>(o));
      }

      /// Operator == also compares with nested data
      bool operator==(ChunkPosition const& o) const {
        return chunk == o.chunk && Base::operator==(static_cast<Base const&>(o));
      }

      /// Operator < also compares with nested data
      bool operator<(ChunkPosition const& o) const {
        return chunk < o.chunk || (chunk == o.chunk && Base::operator<(static_cast<Base const&>(o)));
      }

      /// Operator <= also compares with nested data
      bool operator<=(ChunkPosition const& o) const {
        return chunk < o.chunk || (chunk == o.chunk && Base::operator<=(static_cast<Base const&>(o)));
      }

      /// Operator > also compares with nested data
      bool operator>(ChunkPosition const& o) const {
        return o < *this;
      }

      /// Operator >= also compares with nested data
      bool operator>=(ChunkPosition const& o) const {
        return o <= *this;
      }

      /// Stream output
      friend std::ostream& operator<<(std::ostream& stream, ChunkPosition const& pos) {
        stream << "[" << pos.chunk << ", " << pos.data << ", " << pos.inner << "]";
        return stream;
      }
  };
}
