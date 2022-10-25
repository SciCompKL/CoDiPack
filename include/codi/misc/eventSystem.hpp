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

    enum class Statement {
      Expression,
      Copy,
      Passive
    };
    enum class Reset {
      Full,
      Hard,
      To
    };
  }

  template<typename T_Tape>
  struct EventSystem {
    public:
      using Tape = CODI_DD(T_Tape, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));
      using Real = typename Tape::Real;
      using Gradient = typename Tape::Gradient;
      using Identifier = typename Tape::Identifier;
      using Index = typename Tape::Identifier;
      using Position = typename Tape::Position;
      using VectorAccess = VectorAccessInterface<Real, Identifier>;

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
        StatementPrimal,
        StatementStoreOnTape,
        StatementEvaluate,
        /* index management events */
        IndexCreate,
        IndexAssign,
        IndexFree,
        IndexCopy
      };

      using Callback = void*;
      using EventListenerMap = std::map<Event, std::list<std::pair<Callback, void*>>>;
      static EventListenerMap listeners;

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

      static CODI_INLINE void registerTapeRegisterInputListener(void (*callback)(Tape&, Real&, Identifier&, void*), void* customData = nullptr) {
        internalRegisterListener(Event::TapeRegisterInput, callback, customData);
      }

      static CODI_INLINE void notifyTapeRegisterInputListeners(Tape& tape, Real& value, Identifier& identifier) {
        internalNotifyListeners<void (*)(Tape&, Real&, Identifier&, void*)>(Config::ADWorkflowEvents, Event::TapeRegisterInput, tape, value, identifier);
      }

      static CODI_INLINE void registerTapeRegisterOutputListener(void (*callback)(Tape&, Real&, Identifier&, void*), void* customData = nullptr) {
        internalRegisterListener(Event::TapeRegisterOutput, callback, customData);
      }

      static CODI_INLINE void notifyTapeRegisterOutputListeners(Tape& tape, Real& value, Identifier& identifier) {
        internalNotifyListeners<void (*)(Tape&, Real&, Identifier&, void*)>(Config::ADWorkflowEvents, Event::TapeRegisterOutput, tape, value, identifier);
      }

      static CODI_INLINE void registerTapeEvaluateListener(void (*callback)(Tape&, Position const&, Position const&, VectorAccess*, Events::Direction, Events::Endpoint, void*), void* customData = nullptr) {
        internalRegisterListener(Event::TapeEvaluate, callback, customData);
      }

      static CODI_INLINE void notifyTapeEvaluateListeners(Tape& tape, Position const& start, Position const& end, VectorAccess* adjoint, Events::Direction direction, Events::Endpoint endpoint) {
        internalNotifyListeners<void (*)(Tape&, Position const&, Position const&, VectorAccess*, Events::Direction, Events::Endpoint, void*)>(Config::ADWorkflowEvents, Event::TapeEvaluate, tape, start, end, adjoint, direction, endpoint);
      }

      static CODI_INLINE void registerTapeResetListener(void (*callback)(Tape&, Position const&, Events::Reset, bool, void*), void* customData = nullptr) {
        internalRegisterListener(Event::TapeReset, callback, customData);
      }

      static CODI_INLINE void notifyTapeResetListeners(Tape& tape, Position const& position, Events::Reset kind, bool clearAdjoints) {
        internalNotifyListeners<void (*)(Tape&, Position const&, Events::Reset, bool, void*)>(Config::ADWorkflowEvents, Event::TapeReset, tape, position, kind, clearAdjoints);
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

      static CODI_INLINE void registerPreaccAddInputListener(void (*callback)(Tape&, Real const&, Identifier const&, void*), void* customData = nullptr) {
        internalRegisterListener(Event::PreaccAddInput, callback, customData);
      }

      static CODI_INLINE void notifyPreaccAddInputListeners(Tape& tape, Real const& value, Identifier const& identifier) {
        internalNotifyListeners<void (*)(Tape&, Real const&, Identifier const&, void*)>(Config::PreaccEvents, Event::PreaccAddInput, tape, value, identifier);
      }

      static CODI_INLINE void registerPreaccAddOutputListener(void (*callback)(Tape&, Real&, Identifier&, void*), void* customData = nullptr) {
        internalRegisterListener(Event::PreaccAddOutput, callback, customData);
      }

      static CODI_INLINE void notifyPreaccAddOutputListeners(Tape& tape, Real& value, Identifier& identifier) {
        internalNotifyListeners<void (*)(Tape&, Real&, Identifier&, void*)>(Config::PreaccEvents, Event::PreaccAddOutput, tape, value, identifier);
      }

      /// @}
      /*******************************************************************************/
      /// @name Statements
      /// @{

      static CODI_INLINE void registerStatementPrimalListener(void (*callback)(Tape&, Real const&, Identifier const&, Real const&, Events::Statement, void*), void* customData = nullptr) {
        internalRegisterListener(Event::StatementPrimal, callback, customData);
      }

      static CODI_INLINE void notifyStatementPrimalListeners(Tape& tape, Real const& lhsValue, Identifier const& lhsIdentifier, Real const& rhsValue, Events::Statement statement) {
        internalNotifyListeners<void (*)(Tape&, Real const&, Identifier const&, Real const&, Events::Statement, void*)>(Config::StatementEvents, Event::StatementPrimal, tape, lhsValue, lhsIdentifier, rhsValue, statement);
      }

      static CODI_INLINE void registerStatementStoreOnTapeListener(void (*callback)(Tape&, Identifier const&, Real const&, size_t, Identifier const*, Real const*, void*), void* customData = nullptr) {
        internalRegisterListener(Event::StatementStoreOnTape, callback, customData);
      }

      static CODI_INLINE void notifyStatementStoreOnTapeListeners(Tape& tape, Identifier const& lhsIdentifier, Real const& rhsValue, size_t numActiveVariables, Identifier const* rhsIdentifiers, Real const* jacobians) {
        internalNotifyListeners<void (*)(Tape&, Identifier const&, Real const&, size_t, Identifier const*, Real const*, void*)>(Config::StatementEvents, Event::StatementStoreOnTape, tape, lhsIdentifier, rhsValue, numActiveVariables, rhsIdentifiers, jacobians);
      }

      template<typename Adjoint>
      static CODI_INLINE void registerStatementEvaluateListener(void (*callback)(Tape&, Adjoint const*, Identifier const&, size_t, Identifier const*, Real const*, Events::Direction, void*), void* customData = nullptr) {
        internalRegisterListener<void (*)(Tape&, Adjoint const*, Identifier const&, size_t, Identifier const*, Real const*, Events::Direction, void*)>(Event::StatementEvaluate, callback, customData);
      }

      template<typename Adjoint>
      static CODI_INLINE void notifyStatementEvaluateListeners(Tape& tape, Adjoint const* adjoints, Identifier const& lhsIdentifier, size_t numRhsVariables, Identifier const* rhsIdentifiers, Real const* jacobians, Events::Direction direction) {
        internalNotifyListeners<void (*)(Tape&, Adjoint const*, Identifier const&, size_t, Identifier const*, Real const*, Events::Direction, void*)>(Config::StatementEvents, Event::StatementEvaluate, tape, adjoints, lhsIdentifier, numRhsVariables, rhsIdentifiers, jacobians, direction);
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
  typename EventSystem<Tape>::EventListenerMap EventSystem<Tape>::listeners;
}
