#pragma once

#include <vector>

#include "../../aux/macros.h"
#include "../../config.h"
#include "chunk.hpp"
#include "dataInterface.hpp"
#include "position.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  struct EmptyVector : public DataInterface<> {

      using NestedVector = void;

      using Position = EmptyPosition;

      /*******************************************************************************
       * Section: Misc functions
       *
       * Description: TODO
       *
       */

      CODI_INLINE Position extractPosition(Position const& pos) const {
        return pos;
      }

      CODI_INLINE size_t getDataSize() const {
        return 0;
      }

      CODI_INLINE Position getPosition() const {
        return Position();
      }
      CODI_INLINE Position getZeroPosition() const {
        return Position();
      }

      CODI_INLINE void pushData() {}
      CODI_INLINE void reserveItems(size_t const& items) { CODI_UNUSED(items); }
      void resize(size_t const& totalSize) { CODI_UNUSED(totalSize); }

      void reset() {}
      void resetHard() {}
      void resetTo(Position const& pos) { CODI_UNUSED(pos); }

      void setNested(NestedVector* v) { CODI_UNUSED(v); }

      void swap(DataInterface& other) { CODI_UNUSED(other); }

      /*******************************************************************************
       * Section: Iterator functions
       *
       * Description: TODO
       *
       */

      template<typename Function, typename ... Args>
      CODI_INLINE void evaluateForward(Position const& start, Position const& end,Function const& function,
                                       Args&&... args) {
        CODI_UNUSED(start, end);
        function(std::forward<Args>(args)...);
      }

      template<typename Function, typename ... Args>
      CODI_INLINE void evaluateReverse(Position const& start, Position const& end,Function const& function,
                                       Args&&... args) {
        CODI_UNUSED(start, end);
        function(std::forward<Args>(args)...);
      }

      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachChunk(FunctionObject& function, bool recursive, Args&&... args) {
        CODI_UNUSED(function, recursive, args...);
      }

      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachForward(Position const& start, Position const& end, FunctionObject& function,
                                      Args&&... args) {
        CODI_UNUSED(start, end, function, args...);
      }

      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachReverse(Position const& start, Position const& end, FunctionObject& function,
                                      Args&&... args) {
        CODI_UNUSED(start, end, function, args...);
      }
  };
}
