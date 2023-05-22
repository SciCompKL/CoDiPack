/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <vector>

#include "../../config.h"
#include "../../misc/macros.hpp"
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

      using NestedData = void;  ///< No nested data.

      using Position = EmptyPosition;    ///< No positional data.
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

      /// \copydoc codi::DataInterface::resetTo
      void resetTo(Position const& pos) {
        CODI_UNUSED(pos);
      }

      /// \copydoc DataInterface::erase
      void erase(Position const& start, Position const& end, bool recursive = true) {
        CODI_UNUSED(start, end, recursive);
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

      /// \copybrief DataInterface::getDataPointers
      CODI_INLINE void getDataPointers(InternalPosHandle const& startPos) {
        CODI_UNUSED(startPos);
      }

      /*******************************************************************************/
      /// @name Misc functions

      /// \copydoc DataInterface::addToTapeValues
      void addToTapeValues(TapeValues& values) const {
        CODI_UNUSED(values);
      }

      /// \copydoc DataInterface::extractPosition
      template<typename = void>
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
