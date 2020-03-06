#pragma once

#include <vector>

#include "../../aux/macros.h"
#include "../../config.h"
#include "../aux/tapeValues.hpp"
#include "chunk.hpp"

/** \copydoc codi::Namespace */
namespace codi {


  template<typename _NestedVector = void, typename _InternalPosHandle = size_t>
  struct DataInterface {

      using NestedVector = DECLARE_DEFAULT(_NestedVector, DataInterface);
      using InternalPosHandle = DECLARE_DEFAULT(_InternalPosHandle, size_t);

      using Position = ANY;

      /*******************************************************************************
       * Section: Misc functions
       *
       * Description: TODO
       *
       */

      void addToTapeValues(TapeValues& values) const;

      template<typename TargetPosition>
      CODI_INLINE TargetPosition extractPosition(Position const& pos) const;

      CODI_INLINE size_t getDataSize() const;
      CODI_INLINE Position getPosition() const;
      CODI_INLINE size_t getPushedDataCount(InternalPosHandle const& startPos);
      CODI_INLINE Position getZeroPosition() const;


      template<typename ... Data>
      CODI_INLINE void pushData(Data const& ... data);
      CODI_INLINE InternalPosHandle reserveItems(size_t const& items);
      void resize(size_t const& totalSize);

      void reset();
      void resetHard();
      void resetTo(Position const& pos);

      void setNested(NestedVector* v);

      void swap(DataInterface& other);

      /*******************************************************************************
       * Section: Iterator functions
       *
       * Description: TODO
       *
       */

      template<typename Function, typename ... Args>
      CODI_INLINE void evaluateForward(Position const& start, Position const& end,Function const& function,
                                       Args&&... args);

      template<typename Function, typename ... Args>
      CODI_INLINE void evaluateReverse(Position const& start, Position const& end,Function const& function,
                                       Args&&... args);

      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachChunk(FunctionObject& function, bool recursive, Args&&... args);

      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachForward(Position const& start, Position const& end, FunctionObject& function,
                                      Args&&... args);

      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachReverse(Position const& start, Position const& end, FunctionObject& function,
                                      Args&&... args);
  };
}
