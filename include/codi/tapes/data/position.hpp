#pragma once

#include <cstdlib>
#include <iostream>

#include "../../aux/macros.hpp"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /// Empty Position with no nested data.
  struct EmptyPosition {
    public:

      /// Always false
      bool operator != (EmptyPosition const& o) const {
        CODI_UNUSED(o);
        return false;
      }

      /// Always true
      bool operator == (EmptyPosition const& o) const {
        CODI_UNUSED(o);
        return true;
      }

      /// Always false
      bool operator < (EmptyPosition const& o) const {
        CODI_UNUSED(o);
        return false;
      }

      /// Always true
      bool operator <= (EmptyPosition const& o) const {
        CODI_UNUSED(o);
        return true;
      }

      /// Always false
      bool operator > (EmptyPosition const& o) const {
        CODI_UNUSED(o);
        return false;
      }

      /// Always true
      bool operator >= (EmptyPosition const& o) const {
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
   * Used for data that is allocated en block e.g. BlockData
   *
   * @tparam _NestedPosition  Position implementation
   */
  template<typename _NestedPosition>
  struct ArrayPosition {
    public:

      using NestedPosition = CODI_DD(_NestedPosition, EmptyPosition); ///< See ArrayPosition

      size_t data;  ///< Array position index.

      NestedPosition inner; ///< Position of nested data

      /// Constructor
      ArrayPosition() :
        data(0),
        inner() {}

      /// Constructor
      ArrayPosition(size_t const& data, NestedPosition const& inner) :
        data(data),
        inner(inner) {}

      /// Operator != also compares with nested data
      bool operator != (ArrayPosition const& o) const {
        return data != o.data || this->inner != o.inner;
      }

      /// Operator == also compares with nested data
      bool operator == (ArrayPosition const& o) const {
        return data == o.data && this->inner == o.inner;
      }

      /// Operator < also compares with nested data
      bool operator < (ArrayPosition const& o) const {
        return data < o.data || (data == o.data && inner < o.inner);
      }

      /// Operator <= also compares with nested data
      bool operator <= (ArrayPosition const& o) const {
        return data < o.data || (data == o.data && inner <= o.inner);
      }

      /// Operator > also compares with nested data
      bool operator > (ArrayPosition const& o) const {
        return o < *this;
      }

      /// Operator >= also compares with nested data
      bool operator >= (ArrayPosition const& o) const {
        return o <= *this;
      }

      /// Stream output
      friend std::ostream& operator<<(std::ostream& stream, ArrayPosition const& pos) {
        stream << "[" << pos.data << ", " << pos.inner << "]";
        return stream;
      }
  };

  /**
   * @brief Position with two nested indices for e.g. chunked data access.
   *
   * Used for data that is allocated with multiple chunk e.g. ChunkData.
   * `chunk` is the major one and provides the chunk index and `Base::data` is the secondary one which provides the
   * position in the chunk.
   *
   * For `p1 < p2` it is enough that `p1.chunk < p2.chunk`, only in the equality case `Base::data` needs to be checked.
   *
   * @tparam _NestedPosition  Position implementation
   */
  template<typename _NestedPosition>
  struct ChunkPosition : public ArrayPosition<_NestedPosition> {
    public:

      using NestedPosition = CODI_DD(_NestedPosition, EmptyPosition); ///< See ChunkPosition
      using Base = ArrayPosition<NestedPosition>; ///< Base abbreviation

      size_t chunk; ///< Chunk position index

      /// Constructor
      ChunkPosition() :
        Base(),
        chunk(0) {}

      /// Constructor
      ChunkPosition(size_t const& chunk, size_t const& data, NestedPosition const& inner) :
        Base(data, inner),
        chunk(chunk) {}

      /// Operator != also compares with nested data
      bool operator != (ChunkPosition const& o) const {
        return chunk != o.chunk || Base::operator !=(static_cast<Base const&>(o));
      }

      /// Operator == also compares with nested data
      bool operator == (ChunkPosition const& o) const {
        return chunk == o.chunk && Base::operator ==(static_cast<Base const&>(o));
      }

      /// Operator < also compares with nested data
      bool operator < (ChunkPosition const& o) const {
        return chunk < o.chunk || (chunk == o.chunk && Base::operator <(static_cast<Base const&>(o)));
      }

      /// Operator <= also compares with nested data
      bool operator <= (ChunkPosition const& o) const {
        return chunk < o.chunk || (chunk == o.chunk && Base::operator <=(static_cast<Base const&>(o)));
      }

      /// Operator > also compares with nested data
      bool operator > (ChunkPosition const& o) const {
        return o < *this;
      }

      /// Operator >= also compares with nested data
      bool operator >= (ChunkPosition const& o) const {
        return o <= *this;
      }

      /// Stream output
      friend std::ostream& operator<<(std::ostream& stream, ChunkPosition const& pos) {
        stream << "[" << pos.chunk << ", " << pos.data << ", " << pos.inner << "]";
        return stream;
      }
  };
}
