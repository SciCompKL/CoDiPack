/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 *  - SciComp, TU Kaiserslautern:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <vector>

#include "../../misc/eventSystem.hpp"
#include "../../misc/macros.hpp"
#include "../../config.h"
#include "../data/dataInterface.hpp"
#include "indexManagerInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Identifiers are created in a linear fashion. Each assign creates a new index which is counted up.
   *
   * A simple copy optimization is implemented here. Since each index is only bound to one primal value, the index
   * can simply be copied.
   *
   * Mathematical and implementational details are explained in \ref SBG2021Index.
   *
   * Since this index manager is tightly coupled to the statements, it also implements a simple data interface. It just
   * provides the current maximum index as positional information and adds this positional information in the evaluate
   * routines.
   *
   * @tparam T_Index   Type for the identifier, usually an integer type.
   */
  template<typename T_Index>
  struct LinearIndexManager : public IndexManagerInterface<T_Index>, public DataInterface<> {
    public:

      using Index = CODI_DD(T_Index, int);        ///< See LinearIndexManager.
      using Base = IndexManagerInterface<Index>;  ///< Base class abbreviation.

      /*******************************************************************************/
      /// @name IndexManagerInterface: Constants
      /// @{

      static bool constexpr CopyNeedsStatement = false;  ///< Copy optimization is implemented.
      static bool constexpr IsLinear = true;             ///< Tightly coupled to statements.

      /// @}
      /*******************************************************************************/
      /// @name DataInterface: Type declaration
      /// @{

      using Position = Index;            ///< Positions coincide with indices.
      using NestedData = void;           ///< Terminator, no further nested data.
      using InternalPosHandle = size_t;  ///< Internal positions coincide with positions.

      /// @}

    private:

      Index reservedIndices;  ///< The largest index that is reserved and cannot be assigned to active AD variables.
      Index count;            ///< The current maximum index.

    public:

      /// Constructor
      LinearIndexManager(Index reservedIndices) : reservedIndices(reservedIndices), count(reservedIndices) {}

      /*******************************************************************************/
      /// @name IndexManagerInterface: Methods
      /// @{

      /// \copydoc IndexManagerInterface::addToTapeValues <br><br>
      /// Implementation: Adds maximum live indices.
      void addToTapeValues(TapeValues& values) const {
        values.addLongEntry("Max. live indices", getLargestCreatedIndex());
      }

      /// \copydoc IndexManagerInterface::freeIndex <br><br>
      /// Implementation: Freed indices are ignored.
      template<typename Tape>
      CODI_INLINE void freeIndex(Index& index) const {
        index = Base::InactiveIndex;
      }

      /// \copydoc IndexManagerInterface::assignIndex
      template<typename Tape>
      CODI_INLINE bool assignIndex(Index& index) {
        if (CODI_ENABLE_CHECK(Config::OverflowCheck, count > count + 1)) {
          CODI_EXCEPTION("Overflow in linear index handler. Use a larger index type or a reuse index manager.");
        }
        count += 1;
        EventSystem<Tape>::template notifyListeners<Event::IndexAssign>(count);
        index = count;
        return true;
      }

      /// \copydoc IndexManagerInterface::assignUnusedIndex
      template<typename Tape>
      CODI_INLINE bool assignUnusedIndex(Index& index) {
        return assignIndex<Tape>(index);
      }

      /// \copydoc IndexManagerInterface::copyIndex
      template<typename Tape>
      CODI_INLINE void copyIndex(Index& lhs, Index const& rhs) {
        lhs = rhs;
      }

      /// \copydoc IndexManagerInterface::getLargestCreatedIndex
      /// The following properties are specific to the LinearIndexManager:
      /// 1. tape resets reset the largest created index to zero,
      /// 2. the largest created index coincides with the largest assigned index.
      CODI_INLINE Index getLargestCreatedIndex() const {
        return count;
      }

      /// @}
      /*******************************************************************************/
      /// @name DataInterface: Methods
      /// @{

      /// \copydoc DataInterface::extractPosition
      template<typename TargetPosition>
      CODI_INLINE TargetPosition extractPosition(Position const& pos) const {
        return pos;  // This is a terminator, no further recursion.
      }

      /// \copydoc DataInterface::getDataSize
      CODI_INLINE size_t getDataSize() const {
        return 0;
      }

      /// \copydoc DataInterface::getPosition
      /// The current position coincides with the largest created index.
      CODI_INLINE Position getPosition() const {
        return count;
      }

      /// \copydoc DataInterface::getPushedDataCount
      CODI_INLINE size_t getPushedDataCount(InternalPosHandle const& startPos) {
        return count - startPos;
      }

      /// \copydoc DataInterface::getZeroPosition
      /// The zero position coincides with the smallest index that may be assigned to active AD variables.
      CODI_INLINE Position getZeroPosition() const {
        return reservedIndices;
      }

      /// \copydoc DataInterface::pushData
      CODI_INLINE void pushData() {}

      /// \copydoc DataInterface::reserveItems
      CODI_INLINE InternalPosHandle reserveItems(size_t const& items) {
        CODI_UNUSED(items);

        return count;
      }

      /// \copydoc DataInterface::resize
      void resize(size_t const& totalSize) {
        CODI_UNUSED(totalSize);
      }

      /// \copydoc DataInterface::resetTo
      CODI_INLINE void resetTo(Position const& pos) {
        codiAssert(pos >= reservedIndices);

        count = pos;
      }

      /// \copydoc DataInterface::reset
      CODI_INLINE void reset() {
        count = reservedIndices;
      }

      /// \copydoc DataInterface::resetHard
      CODI_INLINE void resetHard() {
        count = reservedIndices;
      }

      /// \copydoc DataInterface::setNested
      void setNested(NestedData* v) {
        CODI_UNUSED(v);
      }

      /// \copydoc DataInterface::swap
      void swap(LinearIndexManager<Index>& other) {
        std::swap(reservedIndices, other.reservedIndices);
        std::swap(count, other.count);
      }

      /// \copydoc DataInterface::evaluateForward
      template<typename FunctionObject, typename... Args>
      CODI_INLINE void evaluateForward(Position const& start, Position const& end, FunctionObject function,
                                       Args&&... args) {
        function(std::forward<Args>(args)..., start, end);
      }

      /// \copydoc DataInterface::evaluateReverse
      template<typename FunctionObject, typename... Args>
      CODI_INLINE void evaluateReverse(Position const& start, Position const& end, FunctionObject function,
                                       Args&&... args) {
        function(std::forward<Args>(args)..., start, end);
      }

      /// \copydoc DataInterface::forEachChunk
      template<typename FunctionObject, typename... Args>
      CODI_INLINE void forEachChunk(FunctionObject& function, bool recursive, Args&&... args) {
        CODI_UNUSED(function, recursive, args...);
        // Do nothing.
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
