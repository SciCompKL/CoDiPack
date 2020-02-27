#pragma once

#include "../../aux/macros.h"
#include "../../config.h"
#include "chunk.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename ChunkData>
  struct PointerStore {
    public:

      template<typename FuncObj, typename ... Args>
      void call(FuncObj& func, Args&&... args);

      template<typename Nested, typename ... Args>
      CODI_INLINE void callNestedForward(Nested* nested, size_t& start, size_t const& end, Args&&... args);

      template<typename Nested, typename ... Args>
      CODI_INLINE void callNestedReverse(Nested* nested, size_t& start, size_t const& end, Args&&... args);

      void setPointers(const size_t& dataPos, ChunkData* chunk);
  };

  template<typename _Data1>
  struct PointerStore<Chunk1<_Data1> > {
    public:

      using Data1 = DECLARE_DEFAULT(_Data1, int);

      using Chunk = Chunk1<Data1>;

      Data1* p1;

      template<typename FuncObj, typename ... Args>
      void call(FuncObj& func, Args&&... args) {
        func(p1, std::forward<Args>(args)...);
      }

      template<typename Nested, typename ... Args>
      CODI_INLINE void callNestedForward(Nested* nested, size_t& start, size_t const& end, Args&&... args) {
        nested->evaluateForward(std::forward<Args>(args)..., start, end, p1);
      }

      template<typename Nested, typename ... Args>
      CODI_INLINE void callNestedReverse(Nested* nested, size_t& start, size_t const& end, Args&&... args) {
        nested->evaluateReverse(std::forward<Args>(args)..., start, end, p1);
      }

      void setPointers(size_t const& dataPos, Chunk* chunk) {
        chunk->dataPointer(dataPos, p1);
      }
  };

  template<typename _Data1, typename _Data2>
  struct PointerStore<Chunk2<_Data1, _Data2> > {
    public:

      using Data1 = DECLARE_DEFAULT(_Data1, int);
      using Data2 = DECLARE_DEFAULT(_Data2, int);

      using Chunk = Chunk2<Data1, Data2>;

      Data1* p1;
      Data2* p2;

      template<typename FuncObj, typename ... Args>
      void call(FuncObj& func, Args&&... args) {
        func(p1, p2, std::forward<Args>(args)...);
      }

      template<typename Nested, typename ... Args>
      CODI_INLINE void callNestedForward(Nested* nested, size_t& start, size_t const& end, Args&&... args) {
        nested->evaluateForward(std::forward<Args>(args)..., start, end, p1, p2);
      }

      template<typename Nested, typename ... Args>
      CODI_INLINE void callNestedReverse(Nested* nested, size_t& start, size_t const& end, Args&&... args) {
        nested->evaluateReverse(std::forward<Args>(args)..., start, end, p1, p2);
      }

      void setPointers(size_t const& dataPos, Chunk* chunk) {
        chunk->dataPointer(dataPos, p1, p2);
      }
  };

  template<typename _Data1, typename _Data2, typename _Data3>
  struct PointerStore<Chunk3<_Data1, _Data2, _Data3> > {
    public:

      using Data1 = DECLARE_DEFAULT(_Data1, int);
      using Data2 = DECLARE_DEFAULT(_Data2, int);
      using Data3 = DECLARE_DEFAULT(_Data3, int);

      using Chunk = Chunk3<Data1, Data2, Data3>;

      Data1* p1;
      Data2* p2;
      Data3* p3;

      template<typename FuncObj, typename ... Args>
      void call(FuncObj& func, Args&&... args) {
        func(p1, p2, p3, std::forward<Args>(args)...);
      }

      template<typename Nested, typename ... Args>
      CODI_INLINE void callNestedForward(Nested* nested, size_t& start, size_t const& end, Args&&... args) {
        nested->evaluateForward(std::forward<Args>(args)..., start, end, p1, p2, p3);
      }

      template<typename Nested, typename ... Args>
      CODI_INLINE void callNestedReverse(Nested* nested, size_t& start, size_t const& end, Args&&... args) {
        nested->evaluateReverse(std::forward<Args>(args)..., start, end, p1, p2, p3);
      }

      void setPointers(size_t const& dataPos, Chunk* chunk) {
        chunk->dataPointer(dataPos, p1, p2, p3);
      }
  };

  template<typename _Data1, typename _Data2, typename _Data3, typename _Data4>
  struct PointerStore<Chunk4<_Data1, _Data2, _Data3, _Data4> > {
    public:

      using Data1 = DECLARE_DEFAULT(_Data1, int);
      using Data2 = DECLARE_DEFAULT(_Data2, int);
      using Data3 = DECLARE_DEFAULT(_Data3, int);
      using Data4 = DECLARE_DEFAULT(_Data4, int);

      using Chunk = Chunk4<Data1, Data2, Data3, Data4>;

      Data1* p1;
      Data2* p2;
      Data3* p3;
      Data4* p4;

      template<typename FuncObj, typename ... Args>
      void call(FuncObj& func, Args&&... args) {
        func(p1, p2, p3, p4, std::forward<Args>(args)...);
      }

      template<typename Nested, typename ... Args>
      CODI_INLINE void callNestedForward(Nested* nested, size_t& start, size_t const& end, Args&&... args) {
        nested->evaluateForward(std::forward<Args>(args)..., start, end, p1, p2, p3, p4);
      }

      template<typename Nested, typename ... Args>
      CODI_INLINE void callNestedReverse(Nested* nested, size_t& start, size_t const& end, Args&&... args) {
        nested->evaluateReverse(std::forward<Args>(args)..., start, end, p1, p2, p3, p4);
      }

      void setPointers(size_t const& dataPos, Chunk* chunk) {
        chunk->dataPointer(dataPos, p1, p2, p3, p4);
      }
  };
}
