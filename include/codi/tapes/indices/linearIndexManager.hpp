#pragma once

#include <vector>

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../data/dataInterface.hpp"
#include "indexManagerInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Identifiers are created in a linear fashion. Each assign creates a new index which is counted up.
   *
   * A simple assign optimization is implemented here. Since each index is only bound to one primal value, the index
   * can simply be copied.
   *
   * TODO: Paper ref.
   *
   * Since this index manager is tightly coupled to the statements it also implements a simple data interface. It just
   * provides the current maximum index as positional information and adds this positional information in the evaluate
   * routines.
   *
   * @tparam _Index   Type for the identifier usually an integer type.
   */
  template<typename _Index>
  struct LinearIndexManager : public IndexManagerInterface<_Index>, public DataInterface<> {
    public:

      using Index = CODI_DECLARE_DEFAULT(_Index, int); ///< See LinearIndexManager
      using Base = IndexManagerInterface<Index>; ///< Base class abbreviation

      /*******************************************************************************/
      /// @name IndexManagerInterface: Constants
      /// @{

      static bool constexpr CopyNeedsStatement = false;  ///< Copy optimization implemented
      static bool constexpr IsLinear = true; ///< Tightly coupled to statements.

      /// @}
      /*******************************************************************************/
      /// @name DataInterface: Type declaration
      /// @{

      using Position = Index; ///< The current maximum index
      using NestedData = void; ///< Terminator
      using InternalPosHandle = size_t;  ///< The current maximum index

      /// @}

    private:

      Index zeroState;
      Index count;

    public:

      /// Constructor
      LinearIndexManager(Index zeroState) :
        zeroState(zeroState),
        count(zeroState) {}

      /*******************************************************************************/
      /// @name IndexManagerInterface: Methods
      /// @{

      /// \copydoc IndexManagerInterface::addToTapeValues <br><br>
      /// Implementation: Adds maximum live indices.
      void addToTapeValues(TapeValues& values) const {
        values.addLongEntry("Max. live indices", getLargestAssignedIndex());
      }

      /// \copydoc IndexManagerInterface::freeIndex <br><br>
      /// Implementation: Freed indices are ignored.
      CODI_INLINE void freeIndex(Index& index) const {
        index = Base::UnusedIndex;
      }

      /// \copydoc IndexManagerInterface::assignIndex
      CODI_INLINE bool assignIndex(Index& index) {
        CODI_ENABLE_CHECK(Config::OverflowCheck, count > count + 1) {
          CODI_EXCEPTION("Overflow in linear index handler. Use a larger index type or an reuse index manager.");
        }
        count += 1;
        index = count;
        return true;
      }

      /// \copydoc IndexManagerInterface::assignUnusedIndex
      CODI_INLINE bool assignUnusedIndex(Index& index) {
        return assignIndex(index);
      }

      /// \copydoc IndexManagerInterface::copyIndex
      CODI_INLINE void copyIndex(Index& lhs, Index const& rhs) {
        lhs = rhs;
      }

      /// \copydoc IndexManagerInterface::getLargestAssignedIndex
      CODI_INLINE Index getLargestAssignedIndex() const {
        return count;
      }

      /// @}
      /*******************************************************************************/
      /// @name DataInterface: Methods
      /// @{


      /// \copydoc DataInterface::extractPosition
      template<typename TargetPosition>
      CODI_INLINE TargetPosition extractPosition(Position const& pos) const {
        return pos;  // Last in line needs to be this position.
      }

      /// \copydoc DataInterface::getDataSize
      CODI_INLINE size_t getDataSize() const {
        return 0;
      }

      /// \copydoc DataInterface::getPosition
      CODI_INLINE Position getPosition() const {
        return count;
      }

      /// \copydoc DataInterface::getPushedDataCount
      CODI_INLINE size_t getPushedDataCount(InternalPosHandle const& startPos) {
        return count - startPos;
      }

      /// \copydoc DataInterface::getZeroPosition
      CODI_INLINE Position getZeroPosition() const {
        return zeroState;
      }

      /// \copydoc DataInterface::pushData
      CODI_INLINE void pushData() {}

      /// \copydoc DataInterface::reserveItems
      CODI_INLINE InternalPosHandle reserveItems(size_t const& items) {
        CODI_UNUSED(items);

        return count;
      }

      /// \copydoc DataInterface::resize
      void resize(size_t const& totalSize) { CODI_UNUSED(totalSize); }

      /// \copydoc DataInterface::resetTo
      CODI_INLINE void resetTo(Position const& pos) {
        codiAssert(pos >= zeroState);

        count = pos;
      }

      /// \copydoc DataInterface::reset
      CODI_INLINE void reset() {
        count = zeroState;
      }

      /// \copydoc DataInterface::resetHard
      CODI_INLINE void resetHard() {
        count = zeroState;
      }

      /// \copydoc DataInterface::setNested
      void setNested(NestedData* v) { CODI_UNUSED(v); }

      /// \copydoc DataInterface::swap
      void swap(LinearIndexManager<Index>& other) {
        std::swap(zeroState, other.zeroState);
        std::swap(count, other.count);
      }

      /// \copydoc DataInterface::evaluateForward
      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void evaluateForward(Position const& start, Position const& end, FunctionObject function,
                                       Args&&... args) {

        function(std::forward<Args>(args)..., start, end);
      }

      /// \copydoc DataInterface::evaluateReverse
      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void evaluateReverse(Position const& start, Position const& end, FunctionObject function,
                                       Args&&... args) {

        function(std::forward<Args>(args)..., start, end);
      }

      /// \copydoc DataInterface::forEachChunk
      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachChunk(FunctionObject& function, bool recursive, Args&&... args) {
        CODI_UNUSED(function, recursive, args...);
        // Do nothing
      }

      /// \copydoc DataInterface::forEachForward
      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachForward(Position const& start, Position const& end, FunctionObject function,
                                      Args&&... args) {
        CODI_UNUSED(start, end, function, args...);
      }

      /// \copydoc DataInterface::forEachReverse
      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachReverse(Position const& start, Position const& end, FunctionObject function,
                                      Args&&... args) {
        CODI_UNUSED(start, end, function, args...);
      }
  };
}
