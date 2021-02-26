#pragma once

#include <vector>

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "chunk.hpp"
#include "dataInterface.hpp"
#include "position.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief No data is stored in this DataInterface implementation. It is used to terminate the recursive nature of the
   * DataInterface design.
   *
   * See DataInterface documentation for details.
   */
  struct EmptyData : public DataInterface<> {
    public:

      using NestedData = void;  ///< No nested data

      using Position = EmptyPosition;    ///< No positional data
      using InternalPosHandle = size_t;  ///< Will always be zero.

      /*******************************************************************************/
      /// @name Adding items

      /// \copydoc DataInterface::pushData
      CODI_INLINE void pushData() {}

      /// \copydoc DataInterface::reserveItems
      CODI_INLINE InternalPosHandle reserveItems(size_t const& items) {
        CODI_UNUSED(items);

        return 0;
      }

      /*******************************************************************************/
      /// @name Size management

      /// \copydoc DataInterface::resize
      void resize(size_t const& totalSize) {
        CODI_UNUSED(totalSize);
      }

      /// \copydoc DataInterface::reset
      void reset() {}

      /// \copydoc DataInterface::resetHard
      void resetHard() {}

      /// \copydoc DataInterface::resetTo
      void resetTo(Position const& pos) {
        CODI_UNUSED(pos);
      }

      /*******************************************************************************/
      /// @name Position functions

      /// \copydoc DataInterface::getDataSize
      CODI_INLINE size_t getDataSize() const {
        return 0;
      }

      /// \copydoc DataInterface::getPosition
      CODI_INLINE Position getPosition() const {
        return Position();
      }

      /// \copydoc DataInterface::getPushedDataCount
      CODI_INLINE size_t getPushedDataCount(InternalPosHandle const& startPos) {
        CODI_UNUSED(startPos);

        return 0;
      }

      /// \copydoc DataInterface::getZeroPosition
      CODI_INLINE Position getZeroPosition() const {
        return Position();
      }

      /*******************************************************************************/
      /// @name Misc functions

      /// \copydoc DataInterface::addToTapeValues
      void addToTapeValues(TapeValues& values) const {
        CODI_UNUSED(values);
      }

      /// \copydoc DataInterface::extractPosition
      CODI_INLINE Position extractPosition(Position const& pos) const {
        return pos;
      }

      /// \copydoc DataInterface::setNested <br><br>
      /// Implementation: Does not need to be called.
      void setNested(NestedData* v) {
        CODI_UNUSED(v);
      }

      /// \copydoc DataInterface::swap
      void swap(DataInterface& other) {
        CODI_UNUSED(other);
      }

      /*******************************************************************************/
      /// @name Iterator functions

      /// \copydoc DataInterface::evaluateForward <br><br>
      /// Implementation: Calls the function object with all arguments except start and end.
      template<typename FunctionObject, typename... Args>
      CODI_INLINE void evaluateForward(Position const& start, Position const& end, FunctionObject function,
                                       Args&&... args) {
        CODI_UNUSED(start, end);
        function(std::forward<Args>(args)...);
      }

      /// \copydoc DataInterface::evaluateReverse <br><br>
      /// Implementation: Calls the function object with all arguments except start and end.
      template<typename FunctionObject, typename... Args>
      CODI_INLINE void evaluateReverse(Position const& start, Position const& end, FunctionObject function,
                                       Args&&... args) {
        CODI_UNUSED(start, end);
        function(std::forward<Args>(args)...);
      }

      /// \copydoc DataInterface::forEachChunk
      template<typename FunctionObject, typename... Args>
      CODI_INLINE void forEachChunk(FunctionObject& function, bool recursive, Args&&... args) {
        CODI_UNUSED(function, recursive, args...);
      }

      /// \copydoc DataInterface::forEachForward
      template<typename FunctionObject, typename... Args>
      CODI_INLINE void forEachForward(Position const& start, Position const& end, FunctionObject function,
                                      Args&&... args) {
        CODI_UNUSED(start, end, function, args...);
      }

      /// \copydoc DataInterface::forEachReverse
      template<typename FunctionObject, typename... Args>
      CODI_INLINE void forEachReverse(Position const& start, Position const& end, FunctionObject function,
                                      Args&&... args) {
        CODI_UNUSED(start, end, function, args...);
      }
  };
}
