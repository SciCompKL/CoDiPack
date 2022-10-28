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

  namespace EventHints {
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
  struct EventSystemBase {
    public:
      using Tape = CODI_DD(T_Tape, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));
      using Real = typename Tape::Real;
      using Identifier = typename Tape::Identifier;

    protected:
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
      static CODI_INLINE void internalRegisterListener(bool const& enabled, Event event, TypedCallback callback, void* customData) {
        if (enabled) {
          listeners[event].push_back(std::make_pair((void*)callback, customData));
        }
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
      /// @name Statements
      /// @{

      static CODI_INLINE void registerStatementPrimalListener(void (*callback)(Tape&, Real const&, Identifier const&, Real const&, EventHints::Statement, void*), void* customData = nullptr) {
        internalRegisterListener(Config::StatementEvents, Event::StatementPrimal, callback, customData);
      }

      static CODI_INLINE void notifyStatementPrimalListeners(Tape& tape, Real const& lhsValue, Identifier const& lhsIdentifier, Real const& rhsValue, EventHints::Statement statement) {
        internalNotifyListeners<void (*)(Tape&, Real const&, Identifier const&, Real const&, EventHints::Statement, void*)>(Config::StatementEvents, Event::StatementPrimal, tape, lhsValue, lhsIdentifier, rhsValue, statement);
      }

      /// @}
  };

  template<typename Tape>
  typename EventSystemBase<Tape>::EventListenerMap EventSystemBase<Tape>::listeners;

  template<typename T_Tape>
  struct EventSystem : public EventSystemBase<T_Tape> {
    public:
      using Tape = CODI_DD(T_Tape, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));
      using Real = typename Tape::Real;
      using Gradient = typename Tape::Gradient;
      using Identifier = typename Tape::Identifier;
      using Index = typename Tape::Identifier;
      using Position = typename Tape::Position;
      using VectorAccess = VectorAccessInterface<Real, Identifier>;

      using Base = EventSystemBase<Tape>;
      using Event = typename Base::Event;

      /*******************************************************************************/
      /// @name AD workflow
      /// @{

      static CODI_INLINE void registerTapeStartRecordingListener(void (*callback)(Tape&, void*), void* customData = nullptr) {
        Base::template internalRegisterListener(Config::ADWorkflowEvents, Event::TapeStartRecording, callback, customData);
      }

      static CODI_INLINE void notifyTapeStartRecordingListeners(Tape& tape) {
        Base::template internalNotifyListeners<void (*)(Tape&, void*)>(Config::ADWorkflowEvents, Event::TapeStartRecording, tape);
      }

      static CODI_INLINE void registerTapeStopRecordingListener(void (*callback)(Tape&, void*), void* customData = nullptr) {
        Base::template internalRegisterListener(Config::ADWorkflowEvents, Event::TapeStopRecording, callback, customData);
      }

      static CODI_INLINE void notifyTapeStopRecordingListeners(Tape& tape) {
        Base::template internalNotifyListeners<void (*)(Tape&, void*)>(Config::ADWorkflowEvents, Event::TapeStopRecording, tape);
      }

      static CODI_INLINE void registerTapeRegisterInputListener(void (*callback)(Tape&, Real&, Identifier&, void*), void* customData = nullptr) {
        Base::template internalRegisterListener(Config::ADWorkflowEvents, Event::TapeRegisterInput, callback, customData);
      }

      static CODI_INLINE void notifyTapeRegisterInputListeners(Tape& tape, Real& value, Identifier& identifier) {
        Base::template internalNotifyListeners<void (*)(Tape&, Real&, Identifier&, void*)>(Config::ADWorkflowEvents, Event::TapeRegisterInput, tape, value, identifier);
      }

      static CODI_INLINE void registerTapeRegisterOutputListener(void (*callback)(Tape&, Real&, Identifier&, void*), void* customData = nullptr) {
        Base::template internalRegisterListener(Config::ADWorkflowEvents, Event::TapeRegisterOutput, callback, customData);
      }

      static CODI_INLINE void notifyTapeRegisterOutputListeners(Tape& tape, Real& value, Identifier& identifier) {
        Base::template internalNotifyListeners<void (*)(Tape&, Real&, Identifier&, void*)>(Config::ADWorkflowEvents, Event::TapeRegisterOutput, tape, value, identifier);
      }

      static CODI_INLINE void registerTapeEvaluateListener(void (*callback)(Tape&, Position const&, Position const&, VectorAccess*, EventHints::Direction, EventHints::Endpoint, void*), void* customData = nullptr) {
        Base::template internalRegisterListener(Config::ADWorkflowEvents, Event::TapeEvaluate, callback, customData);
      }

      static CODI_INLINE void notifyTapeEvaluateListeners(Tape& tape, Position const& start, Position const& end, VectorAccess* adjoint, EventHints::Direction direction, EventHints::Endpoint endpoint) {
        Base::template internalNotifyListeners<void (*)(Tape&, Position const&, Position const&, VectorAccess*, EventHints::Direction, EventHints::Endpoint, void*)>(Config::ADWorkflowEvents, Event::TapeEvaluate, tape, start, end, adjoint, direction, endpoint);
      }

      static CODI_INLINE void registerTapeResetListener(void (*callback)(Tape&, Position const&, EventHints::Reset, bool, void*), void* customData = nullptr) {
        Base::template internalRegisterListener(Config::ADWorkflowEvents, Event::TapeReset, callback, customData);
      }

      static CODI_INLINE void notifyTapeResetListeners(Tape& tape, Position const& position, EventHints::Reset kind, bool clearAdjoints) {
        Base::template internalNotifyListeners<void (*)(Tape&, Position const&, EventHints::Reset, bool, void*)>(Config::ADWorkflowEvents, Event::TapeReset, tape, position, kind, clearAdjoints);
      }

      /// @}
      /*******************************************************************************/
      /// @name Preaccumulation
      /// @{

      static CODI_INLINE void registerPreaccStartListener(void (*callback)(Tape&, void*), void* customData = nullptr) {
        Base::template internalRegisterListener(Config::PreaccEvents, Event::PreaccStart, callback, customData);
      }

      static CODI_INLINE void notifyPreaccStartListeners(Tape& tape) {
        Base::template internalNotifyListeners<void (*)(Tape&, void*)>(Config::PreaccEvents, Event::PreaccStart, tape);
      }

      static CODI_INLINE void registerPreaccFinishListener(void (*callback)(Tape&, void*), void* customData = nullptr) {
        Base::template internalRegisterListener(Config::PreaccEvents, Event::PreaccFinish, callback, customData);
      }

      static CODI_INLINE void notifyPreaccFinishListeners(Tape& tape) {
        Base::template internalNotifyListeners<void (*)(Tape&, void*)>(Config::PreaccEvents, Event::PreaccFinish, tape);
      }

      static CODI_INLINE void registerPreaccAddInputListener(void (*callback)(Tape&, Real const&, Identifier const&, void*), void* customData = nullptr) {
        Base::template internalRegisterListener(Config::PreaccEvents, Event::PreaccAddInput, callback, customData);
      }

      static CODI_INLINE void notifyPreaccAddInputListeners(Tape& tape, Real const& value, Identifier const& identifier) {
        Base::template internalNotifyListeners<void (*)(Tape&, Real const&, Identifier const&, void*)>(Config::PreaccEvents, Event::PreaccAddInput, tape, value, identifier);
      }

      static CODI_INLINE void registerPreaccAddOutputListener(void (*callback)(Tape&, Real&, Identifier&, void*), void* customData = nullptr) {
        Base::template internalRegisterListener(Config::PreaccEvents, Event::PreaccAddOutput, callback, customData);
      }

      static CODI_INLINE void notifyPreaccAddOutputListeners(Tape& tape, Real& value, Identifier& identifier) {
        Base::template internalNotifyListeners<void (*)(Tape&, Real&, Identifier&, void*)>(Config::PreaccEvents, Event::PreaccAddOutput, tape, value, identifier);
      }

      /// @}
      /*******************************************************************************/
      /// @name Statements
      /// @{

      static CODI_INLINE void registerStatementStoreOnTapeListener(void (*callback)(Tape&, Identifier const&, Real const&, size_t, Identifier const*, Real const*, void*), void* customData = nullptr) {
        Base::template internalRegisterListener(Config::PreaccEvents, Event::StatementStoreOnTape, callback, customData);
      }

      static CODI_INLINE void notifyStatementStoreOnTapeListeners(Tape& tape, Identifier const& lhsIdentifier, Real const& rhsValue, size_t numActiveVariables, Identifier const* rhsIdentifiers, Real const* jacobians) {
        Base::template internalNotifyListeners<void (*)(Tape&, Identifier const&, Real const&, size_t, Identifier const*, Real const*, void*)>(Config::StatementEvents, Event::StatementStoreOnTape, tape, lhsIdentifier, rhsValue, numActiveVariables, rhsIdentifiers, jacobians);
      }

      static CODI_INLINE void registerStatementEvaluateListener(void (*callback)(Tape&, Identifier const&, size_t, Real const*, void*), void* customData = nullptr) {
        Base::template internalRegisterListener<void (*)(Tape&, Identifier const&, size_t, Real const*, void*)>(Config::PreaccEvents, Event::StatementEvaluate, callback, customData);
      }

      static CODI_INLINE void notifyStatementEvaluateListeners(Tape& tape, Identifier const& lhsIdentifier, size_t sizeLhsAdjoint, Real const* lhsAdjoint) {
        Base::template internalNotifyListeners<void (*)(Tape&, Identifier const&, size_t, Real const*, void*)>(Config::StatementEvents, Event::StatementEvaluate, tape, lhsIdentifier, sizeLhsAdjoint, lhsAdjoint);
      }

      /// @}
      /*******************************************************************************/
      /// @name Index handling
      /// @{

      static CODI_INLINE void registerIndexAssignListener(void (*callback)(Index const&, void*), void* customData = nullptr) {
        Base::template internalRegisterListener(Config::IndexEvents, Event::IndexAssign, callback, customData);
      }

      static CODI_INLINE void notifyIndexAssignListeners(Index const& index) {
        Base::template internalNotifyListeners<void (*)(Index const&, void*)>(Config::IndexEvents, Event::IndexAssign, index);
      }

      static CODI_INLINE void registerIndexFreeListener(void (*callback)(Index const&, void*), void* customData = nullptr) {
        Base::template internalRegisterListener(Config::IndexEvents, Event::IndexFree, callback, customData);
      }

      static CODI_INLINE void notifyIndexFreeListeners(Index const& index) {
        Base::template internalNotifyListeners<void (*)(Index const&, void*)>(Config::IndexEvents, Event::IndexFree, index);
      }

      static CODI_INLINE void registerIndexCopyListener(void (*callback)(Index const&, void*), void* customData = nullptr) {
        Base::template internalRegisterListener(Config::IndexEvents, Event::IndexCopy, callback, customData);
      }

      static CODI_INLINE void notifyIndexCopyListeners(Index const& index) {
        Base::template internalNotifyListeners<void (*)(Index const&, void*)>(Config::IndexEvents, Event::IndexCopy, index);
      }

      /// @}

  };

  /* specialization for the forward evaluation "tape" is identical to EventSystemBase */

  template<typename Real, typename Gradient>
  struct ForwardEvaluation;

  template<typename Real, typename Gradient>
  struct EventSystem<ForwardEvaluation<Real, Gradient>> : public EventSystemBase<ForwardEvaluation<Real, Gradient>> {

  };

}
