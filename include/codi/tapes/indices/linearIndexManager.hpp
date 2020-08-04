#pragma once

#include <vector>

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../data/dataInterface.hpp"
#include "indexManagerInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Index>
  struct LinearIndexManager : public IndexManagerInterface<_Index>, public DataInterface<> {
    public:

      using Index = CODI_DECLARE_DEFAULT(_Index, int);

      static bool constexpr AssignNeedsStatement = false;
      static bool constexpr IsLinear = true;

      using Position = Index;
      using NestedData = void;
      using InternalPosHandle = size_t;

    private:

      Index zeroState;
      Index count;

    public:

      // sections with respect to interfaces

      LinearIndexManager(Index zeroState) :
        zeroState(zeroState),
        count(zeroState) {}

      /*******************************************************************************
       * Section: Implementation of IndexHandlerInterface
       *
       * Description: TODO
       *
       */

      void addToTapeValues(TapeValues& values) const {
        values.addLongEntry("Max. live indices", getLargestAssignedIndex());
      }

      CODI_INLINE void freeIndex(Index& index) const {
        index = 0;
      }

      CODI_INLINE bool assignIndex(Index& index) {
        CODI_ENABLE_CHECK(Config::OverflowCheck, count > count + 1) {
          CODI_EXCEPTION("Overflow in linear index handler. Use a larger index type or an reuse index manager.");
        }
        count += 1;
        index = count;
        return true;
      }

      CODI_INLINE bool assignUnusedIndex(Index& index) {
        return assignIndex(index);
      }

      CODI_INLINE void copyIndex(Index& lhs, Index const& rhs) {
        lhs = rhs;
      }

      CODI_INLINE Index getLargestAssignedIndex() const {
        return count;
      }

      /*******************************************************************************
       * Section: Implementation of DataInterface
       *
       * Description: TODO
       *
       */

      CODI_INLINE size_t getDataSize() const {
        return 0;
      }

      CODI_INLINE Position getPosition() const {
        return count;
      }

      CODI_INLINE size_t getPushedDataCount(InternalPosHandle const& startPos) {
        return count - startPos;
      }

      CODI_INLINE Position getZeroPosition() const {
        return zeroState;
      }

      CODI_INLINE void pushData() {}

      CODI_INLINE InternalPosHandle reserveItems(size_t const& items) {
        CODI_UNUSED(items);

        return count;
      }

      void resize(size_t const& totalSize) { CODI_UNUSED(totalSize); }

      CODI_INLINE void resetTo(Position const& pos) {
        codiAssert(pos >= zeroState);

        count = pos;
      }

      CODI_INLINE void reset() {
        count = zeroState;
      }

      CODI_INLINE void resetHard() {
        count = zeroState;
      }

      void setNested(NestedData* v) { CODI_UNUSED(v); }

      void swap(LinearIndexManager<Index>& other) {
        std::swap(zeroState, other.zeroState);
        std::swap(count, other.count);
      }

      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void evaluateForward(Position const& start, Position const& end, FunctionObject function,
                                       Args&&... args) {

        function(std::forward<Args>(args)..., start, end);
      }

      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void evaluateReverse(Position const& start, Position const& end, FunctionObject function,
                                       Args&&... args) {

        function(std::forward<Args>(args)..., start, end);
      }

      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachChunk(FunctionObject& function, bool recursive, Args&&... args) {
        CODI_UNUSED(function, recursive, args...);
        // Do nothing
      }

      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachForward(Position const& start, Position const& end, FunctionObject function,
                                      Args&&... args) {
        CODI_UNUSED(start, end, function, args...);
      }

      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachReverse(Position const& start, Position const& end, FunctionObject function,
                                      Args&&... args) {
        CODI_UNUSED(start, end, function, args...);
      }
  };
}
