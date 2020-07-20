#pragma once

#include <cstdlib>
#include <iostream>

#include "../../aux/macros.h"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  struct EmptyPosition {
    public:

      bool operator != (EmptyPosition const& o) const {
        CODI_UNUSED(o);
        return false;
      }

      bool operator == (EmptyPosition const& o) const {
        CODI_UNUSED(o);
        return true;
      }

      bool operator < (EmptyPosition const& o) const {
        CODI_UNUSED(o);
        return false;
      }

      bool operator <= (EmptyPosition const& o) const {
        CODI_UNUSED(o);
        return true;
      }

      bool operator > (EmptyPosition const& o) const {
        CODI_UNUSED(o);
        return false;
      }

      bool operator >= (EmptyPosition const& o) const {
        CODI_UNUSED(o);
        return true;
      }

      friend std::ostream& operator<<(std::ostream& stream, EmptyPosition const& pos) {
        CODI_UNUSED(pos);
        stream << "[]";
        return stream;
      }
  };

  template<typename _NestedPosition>
  struct ArrayPosition {
    public:

      using NestedPosition = DECLARE_DEFAULT(_NestedPosition, EmptyPosition);

      size_t data;

      NestedPosition inner;

      ArrayPosition() :
        data(0),
        inner() {}

      ArrayPosition(size_t const& data, NestedPosition const& inner) :
        data(data),
        inner(inner) {}

      bool operator != (ArrayPosition const& o) const {
        return data != o.data || this->inner != o.inner;
      }

      bool operator == (ArrayPosition const& o) const {
        return data == o.data && this->inner == o.inner;
      }

      bool operator < (ArrayPosition const& o) const {
        return data < o.data || (data == o.data && inner < o.inner);
      }

      bool operator <= (ArrayPosition const& o) const {
        return data < o.data || (data == o.data && inner <= o.inner);
      }

      bool operator > (ArrayPosition const& o) const {
        return o < *this;
      }

      bool operator >= (ArrayPosition const& o) const {
        return o <= *this;
      }

      friend std::ostream& operator<<(std::ostream& stream, ArrayPosition const& pos) {
        stream << "[" << pos.data << ", " << pos.inner << "]";
        return stream;
      }
  };

  template<typename _NestedPosition>
  struct ChunkPosition : public ArrayPosition<_NestedPosition> {
    public:

      using NestedPosition = DECLARE_DEFAULT(_NestedPosition, EmptyPosition);
      using Base = ArrayPosition<NestedPosition>;

      size_t chunk;

      ChunkPosition() :
        Base(),
        chunk(0) {}

      ChunkPosition(size_t const& chunk, size_t const& data, NestedPosition const& inner) :
        Base(data, inner),
        chunk(chunk) {}

      bool operator != (ChunkPosition const& o) const {
        return chunk != o.chunk || Base::operator !=(static_cast<Base const&>(o));
      }

      bool operator == (ChunkPosition const& o) const {
        return chunk == o.chunk && Base::operator ==(static_cast<Base const&>(o));
      }

      bool operator < (ChunkPosition const& o) const {
        return chunk < o.chunk || (chunk == o.chunk && Base::operator <(static_cast<Base const&>(o)));
      }

      bool operator <= (ChunkPosition const& o) const {
        return chunk < o.chunk || (chunk == o.chunk && Base::operator <=(static_cast<Base const&>(o)));
      }

      bool operator > (ChunkPosition const& o) const {
        return o < *this;
      }

      bool operator >= (ChunkPosition const& o) const {
        return o <= *this;
      }

      friend std::ostream& operator<<(std::ostream& stream, ChunkPosition const& pos) {
        stream << "[" << pos.chunk << ", " << pos.data << ", " << pos.inner << "]";
        return stream;
      }
  };
}
