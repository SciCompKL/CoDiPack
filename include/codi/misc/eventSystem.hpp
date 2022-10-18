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

#include <list>
#include <map>

#include "../tapes/interfaces/fullTapeInterface.hpp"
#include "../config.h"
#include "macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  namespace Events {
    enum class Direction {
      Forward,
      Reverse
    };

    enum class Endpoint {
      Begin,
      End
    };
  }

  template<typename T_Tape>
  struct EventSystem {
    public:
      using Tape = CODI_DD(T_Tape, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));
      using Real = typename Tape::Real;
      using Gradient = typename Tape::Gradient;
      using Index = typename Tape::Identifier;
      using Position = typename Tape::Position;

    private:
      enum class Event {
        /* AD workflow */
        TapeStartRecording,
        TapeStopRecording,
        TapeRegisterInput,
        TapeRegisterOutput,
        TapeEvaluate,
        TapeReset,
        /* preaccumulation */
        PreaccStart,
        PreaccFinish,
        PreaccAddInput,
        PreaccAddOutput,
        /* statement events */
        /* ... */
        /* index management events */
        IndexCreate,
        IndexAssign,
        IndexFree,
        IndexCopy
      };

      using Callback = void*;
      static std::map<Event, std::list<std::pair<Callback, void*>>> listeners;

      template<typename TypedCallback>
      static CODI_INLINE void internalRegisterListener(Event event, TypedCallback callback, void* customData) {
        listeners[event].push_back(std::make_pair((void*)callback, customData));
      }

      template<typename TypedCallback, typename... Args>
      static CODI_INLINE void internalNotifyListeners(bool const& enabled, Event event, Args&&... args) {
        if (enabled) {
          for (auto const& listener : listeners[event]) {
            ((TypedCallback) listener.first)(std::forward<Args>(args)..., listener.second);
          }
        }
      }

    public:

      /*******************************************************************************/
      /// @name AD workflow
      /// @{

      static CODI_INLINE void registerTapeStartRecordingListener(void (*callback)(Tape&, void*), void* customData = nullptr) {
        internalRegisterListener(Event::TapeStartRecording, callback, customData);
      }

      static CODI_INLINE void notifyTapeStartRecordingListeners(Tape& tape) {
        internalNotifyListeners<void (*)(Tape&, void*)>(Config::ADWorkflowEvents, Event::TapeStartRecording, tape);
      }

      static CODI_INLINE void registerTapeStopRecordingListener(void (*callback)(Tape&, void*), void* customData = nullptr) {
        internalRegisterListener(Event::TapeStopRecording, callback, customData);
      }

      static CODI_INLINE void notifyTapeStopRecordingListeners(Tape& tape) {
        internalNotifyListeners<void (*)(Tape&, void*)>(Config::ADWorkflowEvents, Event::TapeStopRecording, tape);
      }

      template<typename Lhs>
      static CODI_INLINE void registerTapeRegisterInputListener(void (*callback)(Tape&, Lhs&, void*), void* customData = nullptr) {
        internalRegisterListener(Event::TapeRegisterInput, callback, customData);
      }

      template<typename Lhs>
      static CODI_INLINE void notifyTapeRegisterInputListeners(Tape& tape, Lhs& value) {
        internalNotifyListeners<void (*)(Tape&, Lhs&, void*)>(Config::ADWorkflowEvents, Event::TapeRegisterInput, tape, value);
      }

      template<typename Lhs>
      static CODI_INLINE void registerTapeRegisterOutputListener(void (*callback)(Tape&, Lhs&, void*), void* customData = nullptr) {
        internalRegisterListener(Event::TapeRegisterOutput, callback, customData);
      }

      template<typename Lhs>
      static CODI_INLINE void notifyTapeRegisterOutputListeners(Tape& tape, Lhs& value) {
        internalNotifyListeners<void (*)(Tape&, Lhs&, void*)>(Config::ADWorkflowEvents, Event::TapeRegisterOutput, tape, value);
      }

      template<typename Adjoint>
      static CODI_INLINE void registerTapeEvaluateListener(void (*callback)(Tape&, Position const&, Position const&, Adjoint*, Events::Direction, Events::Endpoint, void*), void* customData = nullptr) {
        internalRegisterListener(Event::TapeEvaluate, callback, customData);
      }

      template<typename Adjoint>
      static CODI_INLINE void notifyTapeEvaluateListeners(Tape& tape, Position const& start, Position const& end, Adjoint* adjoint, Events::Direction direction, Events::Endpoint endpoint) {
        internalNotifyListeners<void (*)(Tape&, Position const&, Position const&, Adjoint*, Events::Direction, Events::Endpoint, void*)>(Config::ADWorkflowEvents, Event::TapeEvaluate, tape, start, end, adjoint, direction, endpoint);
      }

      static CODI_INLINE void registerTapeResetListener(void (*callback)(Tape&, Position const&, bool, void*), void* customData = nullptr) {
        internalRegisterListener(Event::TapeReset, callback, customData);
      }

      static CODI_INLINE void notifyTapeResetListeners(Tape& tape, Position const& position, bool clearAdjoints) {
        internalNotifyListeners<void (*)(Tape&, Position const&, bool, void*)>(Config::ADWorkflowEvents, Event::TapeReset, tape, position, clearAdjoints);
      }

      /// @}
      /*******************************************************************************/
      /// @name Preaccumulation
      /// @{

      static CODI_INLINE void registerPreaccStartListener(void (*callback)(Tape&, void*), void* customData = nullptr) {
        internalRegisterListener(Event::PreaccStart, callback, customData);
      }

      static CODI_INLINE void notifyPreaccStartListeners(Tape& tape) {
        internalNotifyListeners<void (*)(Tape&, void*)>(Config::PreaccEvents, Event::PreaccStart, tape);
      }

      static CODI_INLINE void registerPreaccFinishListener(void (*callback)(Tape&, void*), void* customData = nullptr) {
        internalRegisterListener(Event::PreaccFinish, callback, customData);
      }

      static CODI_INLINE void notifyPreaccFinishListeners(Tape& tape) {
        internalNotifyListeners<void (*)(Tape&, void*)>(Config::PreaccEvents, Event::PreaccFinish, tape);
      }

      template<typename Lhs>
      static CODI_INLINE void registerPreaccAddInputListener(void (*callback)(Tape&, Lhs&, void*), void* customData = nullptr) {
        internalRegisterListener(Event::PreaccAddInput, callback, customData);
      }

      template<typename Lhs>
      static CODI_INLINE void notifyPreaccAddInputListeners(Tape& tape, Lhs& value) {
        internalNotifyListeners<void (*)(Tape&, Lhs&, void*)>(Config::PreaccEvents, Event::PreaccAddInput, tape, value);
      }

      template<typename Lhs>
      static CODI_INLINE void registerPreaccAddOutputListener(void (*callback)(Tape&, Lhs&, void*), void* customData = nullptr) {
        internalRegisterListener(Event::PreaccAddOutput, callback, customData);
      }

      template<typename Lhs>
      static CODI_INLINE void notifyPreaccAddOutputListeners(Tape& tape, Lhs& value) {
        internalNotifyListeners<void (*)(Tape&, Lhs&, void*)>(Config::PreaccEvents, Event::PreaccAddOutput, tape, value);
      }

      /// @}
      /*******************************************************************************/
      /// @name Index handling
      /// @{

      static CODI_INLINE void registerIndexAssignListener(void (*callback)(Index const&, void*), void* customData = nullptr) {
        internalRegisterListener(Event::IndexAssign, callback, customData);
      }

      static CODI_INLINE void notifyIndexAssignListeners(Index const& index) {
        internalNotifyListeners<void (*)(Index const&, void*)>(Config::IndexEvents, Event::IndexAssign, index);
      }

      static CODI_INLINE void registerIndexFreeListener(void (*callback)(Index const&, void*), void* customData = nullptr) {
        internalRegisterListener(Event::IndexFree, callback, customData);
      }

      static CODI_INLINE void notifyIndexFreeListeners(Index const& index) {
        internalNotifyListeners<void (*)(Index const&, void*)>(Config::IndexEvents, Event::IndexFree, index);
      }

      static CODI_INLINE void registerIndexCopyListener(void (*callback)(Index const&, void*), void* customData = nullptr) {
        internalRegisterListener(Event::IndexCopy, callback, customData);
      }

      static CODI_INLINE void notifyIndexCopyListeners(Index const& index) {
        internalNotifyListeners<void (*)(Index const&, void*)>(Config::IndexEvents, Event::IndexCopy, index);
      }

      /// @}
  };

  template<typename Tape>
  std::map<typename EventSystem<Tape>::Event, std::list<std::pair<typename EventSystem<Tape>::Callback, void*>>> EventSystem<Tape>::listeners;
}
