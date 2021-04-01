#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "chunk.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Inserts data pointers at the back of all arguments in the nested call hierarchy.
   *
   * Used in DataInterface implementations for the generalized call to object functions. This class is mainly used in
   * the evaluate and forEach method implementations of these classes. They have to call the evaluation functions on the
   * nested DataInterfaces but they do not know about the data layout in the chunks. This class is the bridge for this
   * call. It stores the data pointers to the chunks and inserts these data pointers in the call to the nested
   * objects.
   *
   * First a call to setPointers has to made and then either call(), callNestedForward() or callNestedReverse().
   *
   * \code{.cpp}
   *   Chunk2<double, int> data(100);
   *
   *   PointerStore<Chunk2<double, int>> ps;
   *
   *   ps.setPointers(0, data);
   *   // one of
   *   ps.call(func, user); // will call func(p1 (double*), p2 (int*), user);
   *   ps.callNestedForward(nested, 0, 10, func, user); // will call nested->evaluateForward(0, 10, func, user,
   *                                                    //                                   p1 (double*), p2 (int*));
   *   ps.callNestedReverse(nested, 10, 0, func, user); // will call nested->evaluateReverse(10, 0, func, user,
   *                                                    //                                   p1 (double*), p2 (int*));
   * \endcode
   *
   * @tparam _ChunkData  Implementation of ChunkBase.
   */
  template<typename _ChunkData>
  struct PointerStore {
    public:

      using ChunkData = CODI_DD(_ChunkData, ChunkBase); ///< See PointerStore.

      /// Calls func(pointers, args...);
      template<typename FuncObj, typename... Args>
      void call(FuncObj& func, Args&&... args);

      /// Calls nested->evaluateForward(args..., start, end, pointers);
      template<typename Nested, typename... Args>
      CODI_INLINE void callNestedForward(Nested* nested, size_t& start, size_t const& end, Args&&... args);

      /// Calls nested->evaluateReverse(args..., start, end, pointers);
      template<typename Nested, typename... Args>
      CODI_INLINE void callNestedReverse(Nested* nested, size_t& start, size_t const& end, Args&&... args);

      /// Sets the internal pointers to the data of the chunk. Afterwards on of the call functions can be called.
      void setPointers(const size_t& dataPos, ChunkData* chunk);
  };

  /**
   * @brief Pointer store for Chunk1 data.
   *
   * See PointerStore for details.
   */
  template<typename _Data1>
  struct PointerStore<Chunk1<_Data1> > {
    public:

      using Data1 = CODI_DD(_Data1, int);  ///< Data entry 1.

      using Chunk = Chunk1<Data1>;  ///< Template specialization type.

    private:

      Data1* p1;  ///< Internal pointer store.

    public:

      /// \copydoc PointerStore::call
      template<typename FuncObj, typename... Args>
      void call(FuncObj& func, Args&&... args) {
        func(p1, std::forward<Args>(args)...);
      }

      /// \copydoc PointerStore::callNestedForward
      template<typename Nested, typename... Args>
      CODI_INLINE void callNestedForward(Nested* nested, size_t& start, size_t const& end, Args&&... args) {
        nested->evaluateForward(std::forward<Args>(args)..., start, end, p1);
      }

      /// \copydoc PointerStore::callNestedReverse
      template<typename Nested, typename... Args>
      CODI_INLINE void callNestedReverse(Nested* nested, size_t& start, size_t const& end, Args&&... args) {
        nested->evaluateReverse(std::forward<Args>(args)..., start, end, p1);
      }

      /// \copydoc PointerStore::setPointers
      void setPointers(size_t const& dataPos, Chunk* chunk) {
        chunk->dataPointer(dataPos, p1);
      }
  };

  /**
   * @brief Pointer store for Chunk2 data.
   *
   * See PointerStore for details.
   */
  template<typename _Data1, typename _Data2>
  struct PointerStore<Chunk2<_Data1, _Data2> > {
    public:

      using Data1 = CODI_DD(_Data1, int);  ///< Data entry 1.
      using Data2 = CODI_DD(_Data2, int);  ///< Data entry 2.

      using Chunk = Chunk2<Data1, Data2>;  ///< Template specialization type.

    private:
      Data1* p1;  ///< Internal pointer store.
      Data2* p2;  ///< Internal pointer store.

    public:

      /// \copydoc PointerStore::call
      template<typename FuncObj, typename... Args>
      void call(FuncObj& func, Args&&... args) {
        func(p1, p2, std::forward<Args>(args)...);
      }

      /// \copydoc PointerStore::callNestedForward
      template<typename Nested, typename... Args>
      CODI_INLINE void callNestedForward(Nested* nested, size_t& start, size_t const& end, Args&&... args) {
        nested->evaluateForward(std::forward<Args>(args)..., start, end, p1, p2);
      }

      /// \copydoc PointerStore::callNestedReverse
      template<typename Nested, typename... Args>
      CODI_INLINE void callNestedReverse(Nested* nested, size_t& start, size_t const& end, Args&&... args) {
        nested->evaluateReverse(std::forward<Args>(args)..., start, end, p1, p2);
      }

      /// \copydoc PointerStore::setPointers
      void setPointers(size_t const& dataPos, Chunk* chunk) {
        chunk->dataPointer(dataPos, p1, p2);
      }
  };

  /**
   * @brief Pointer store for Chunk3 data.
   *
   * See PointerStore for details.
   */
  template<typename _Data1, typename _Data2, typename _Data3>
  struct PointerStore<Chunk3<_Data1, _Data2, _Data3> > {
    public:

      using Data1 = CODI_DD(_Data1, int);  ///< Data entry 1.
      using Data2 = CODI_DD(_Data2, int);  ///< Data entry 2.
      using Data3 = CODI_DD(_Data3, int);  ///< Data entry 3.

      using Chunk = Chunk3<Data1, Data2, Data3>;  ///< Template specialization type.

    private:
      Data1* p1;  ///< Internal pointer store.
      Data2* p2;  ///< Internal pointer store.
      Data3* p3;  ///< Internal pointer store.

    public:

      /// \copydoc PointerStore::call
      template<typename FuncObj, typename... Args>
      void call(FuncObj& func, Args&&... args) {
        func(p1, p2, p3, std::forward<Args>(args)...);
      }

      /// \copydoc PointerStore::callNestedForward
      template<typename Nested, typename... Args>
      CODI_INLINE void callNestedForward(Nested* nested, size_t& start, size_t const& end, Args&&... args) {
        nested->evaluateForward(std::forward<Args>(args)..., start, end, p1, p2, p3);
      }

      /// \copydoc PointerStore::callNestedReverse
      template<typename Nested, typename... Args>
      CODI_INLINE void callNestedReverse(Nested* nested, size_t& start, size_t const& end, Args&&... args) {
        nested->evaluateReverse(std::forward<Args>(args)..., start, end, p1, p2, p3);
      }

      /// \copydoc PointerStore::setPointers
      void setPointers(size_t const& dataPos, Chunk* chunk) {
        chunk->dataPointer(dataPos, p1, p2, p3);
      }
  };

  /**
   * @brief Pointer store for Chunk4 data.
   *
   * See PointerStore for details.
   */
  template<typename _Data1, typename _Data2, typename _Data3, typename _Data4>
  struct PointerStore<Chunk4<_Data1, _Data2, _Data3, _Data4> > {
    public:

      using Data1 = CODI_DD(_Data1, int);  ///< Data entry 1.
      using Data2 = CODI_DD(_Data2, int);  ///< Data entry 2.
      using Data3 = CODI_DD(_Data3, int);  ///< Data entry 3.
      using Data4 = CODI_DD(_Data4, int);  ///< Data entry 4.

      using Chunk = Chunk4<Data1, Data2, Data3, Data4>;  ///< Template specialization type.

    private:
      Data1* p1;  ///< Internal pointer store.
      Data2* p2;  ///< Internal pointer store.
      Data3* p3;  ///< Internal pointer store.
      Data4* p4;  ///< Internal pointer store.

    public:

      /// \copydoc PointerStore::call
      template<typename FuncObj, typename... Args>
      void call(FuncObj& func, Args&&... args) {
        func(p1, p2, p3, p4, std::forward<Args>(args)...);
      }

      /// \copydoc PointerStore::callNestedForward
      template<typename Nested, typename... Args>
      CODI_INLINE void callNestedForward(Nested* nested, size_t& start, size_t const& end, Args&&... args) {
        nested->evaluateForward(std::forward<Args>(args)..., start, end, p1, p2, p3, p4);
      }

      /// \copydoc PointerStore::callNestedReverse
      template<typename Nested, typename... Args>
      CODI_INLINE void callNestedReverse(Nested* nested, size_t& start, size_t const& end, Args&&... args) {
        nested->evaluateReverse(std::forward<Args>(args)..., start, end, p1, p2, p3, p4);
      }

      /// \copydoc PointerStore::setPointers
      void setPointers(size_t const& dataPos, Chunk* chunk) {
        chunk->dataPointer(dataPos, p1, p2, p3, p4);
      }
  };
}
