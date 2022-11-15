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

#include "../config.h"
#include "../tapes/interfaces/fullTapeInterface.hpp"
#include "macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * Namespace for
   */
  namespace EventHints {
    /// Classify tape evaluations.
    enum class EvaluationKind {
      Primal,
      Forward,
      Reverse
    };

    /// Distinguish between beginning and end of tape evaluations.
    enum class Endpoint {
      Begin,
      End
    };

    /// Classify statements.
    enum class Statement {
      Expression,
      Copy,
      Passive
    };

    /// Characterize the tape reset.
    enum class Reset {
      Full,
      Hard,
      To
    };
  }

  /**
   * @brief Basis class for the CoDiPack event system.
   *
   * CoDiPack provides an event system that can be used to, e.g.,
   * - gain insight into CoDiPack's internal workflow,
   * - debug AD, from the coarse AD workflow to the individual statements,
   * - insert custom functionality into CoDiPack's workflow,
   * - monitor the performance of individual AD constructs.
   *
   * For this, a set of events is defined and for each event, custom callbacks can be registered. As the event occurs,
   * CoDiPack invokes the custom callbacks and passes them details about the event itself and related AD data.
   *
   * A callback is registered by a register*Listener call and is subsequently invoked by CoDiPack by the corresponding
   * notify*Listeners call. Please refer to the individual register*Listener functions for the callback signatures.
   * Each callback can be associated with custom data that is provided to each callback invocation. This can be used,
   * e.g., to register the same callback function multiple times with different data. The functions are written such
   * that the required callback signatures should be printed by code completion tools and IDE tooltips.
   *
   * The event system is a tape-specific, global entitiy that is shared by all tapes of the same type. Different tape
   * types use different event system, e.g., second order types have different event systems for outer and inner tapes.
   *
   * This base class defines methods for the StatementPrimal event that is common to forward and reverse tapes.
   *
   * @tparam T_Tape Tape type associated with the event system.
   */
  template<typename T_Tape>
  struct EventSystemBase {
    public:
      /// See EventSystemBase.
      using Tape = CODI_DD(T_Tape, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));
      using Real = typename Tape::Real;              ///< Floating point type the tape is based on.
      using Identifier = typename Tape::Identifier;  ///< Identifier type used by the tape.

    protected:
      /// Full set of events.
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
        StatementEvaluatePrimal,
        /* index management events */
        IndexCreate,
        IndexAssign,
        IndexFree,
        IndexCopy
      };

      using Callback = void*;  ///< Internal, typeless callback storage.
      /// Map that links events and registered callbacks.
      using EventListenerMap = std::map<Event, std::list<std::pair<Callback, void*>>>;

      /**
       * @brief Access the static EventListenerMap.
       *
       * Both tapes and event systems are static entities in CoDiPack, but the tape depends on the event system. We
       * ensure with an initialize-on-first-use pattern that the event system is available when needed.
       */
      static EventListenerMap& getListeners() {
        static EventListenerMap* const listeners = new EventListenerMap;
        return *listeners;
      }

      /**
       * @brief Internal method for callback registration.
       *
       * Stores the callback together with customData in the event entry of the static EventListenerMap.
       *
       * @param enabled     Whether or not the event is active, obtained from Config.
       * @param event       The event for which we register a callback.
       * @param callback    The callback to register.
       * @param customData  A pointer to custom data that is provided whenever the callback is invoked.
       * @tparam TypedCallback  Type of the callback to register.
       */
      template<typename TypedCallback>
      static CODI_INLINE void internalRegisterListener(bool const& enabled, Event event, TypedCallback callback,
                                                       void* customData) {
        if (enabled) {
          getListeners()[event].push_back(std::make_pair((void*)callback, customData));
        }
      }

      /**
       * @brief Internal method for callback invocation.
       *
       * Invokes all callbacks stored for the event entry in the static EventListenerMap.
       * Passes custom data to the callback.
       *
       * @param enabled  Whether or not the event is active, obtained from Config.
       * @param event    The event for which we register a callback.
       * @param args     Arguments for the callback.
       * @tparam TypedCallback  Type of the callback to invoke.
       * @tparam Args           Types of the callback arguments.
       */
      template<typename TypedCallback, typename... Args>
      static CODI_INLINE void internalNotifyListeners(bool const& enabled, Event event, Args&&... args) {
        if (enabled) {
          for (auto const& listener : getListeners()[event]) {
            ((TypedCallback)listener.first)(std::forward<Args>(args)..., listener.second);
          }
        }
      }

    public:

      /*******************************************************************************/
      /// @name Statements
      /// @{

      /**
       * @brief Register callbacks for StatementPrimal events.
       *
       * See notifyStatementPrimalListeners for the callback parameters and event order.
       *
       * @param callback    Callback to be invoked.
       * @param customData  Optional. Custom data that should be linked with the callback, otherwise nullptr.
       */
      static CODI_INLINE void registerStatementPrimalListener(void (*callback)(Tape&, Real const&, Identifier const&,
                                                                               Real const&, EventHints::Statement,
                                                                               void*),
                                                              void* customData = nullptr) {
        internalRegisterListener(Config::StatementEvents, Event::StatementPrimal, callback, customData);
      }

      /**
       * @brief Invoke callbacks for StatementPrimal events.
       *
       * A StatementPrimal event is triggered whenever a code statement x = f(y) with some ActiveType x is executed.
       * The callback is invoked before the assigment is executed.
       *
       * @param tape           Reference to the tape.
       * @param lhsValue       Value of the left hand side before the assignment.
       * @param lhsIdentifier  Identifier (or gradient in forward mode) of the left hand side before the assignment.
       * @param newValue       Value of the right hand side that is assigned.
       * @param statement      Type of statement.
       */
      static CODI_INLINE void notifyStatementPrimalListeners(Tape& tape, Real const& lhsValue,
                                                             Identifier const& lhsIdentifier, Real const& newValue,
                                                             EventHints::Statement statement) {
        internalNotifyListeners<void (*)(Tape&, Real const&, Identifier const&, Real const&, EventHints::Statement,
                                         void*)>(Config::StatementEvents, Event::StatementPrimal, tape, lhsValue,
                                                 lhsIdentifier, newValue, statement);
      }

      /// @}
  };

  /**
   * @brief Full EventSystem implementation for reverse tapes.
   *
   * See EventSystemBase for a general description of the event system.
   *
   * This class implements methods for all of the remaining events.
   *
   * @tparam T_Tape  Tape type associated with the event system.
   */
  template<typename T_Tape>
  struct EventSystem : public EventSystemBase<T_Tape> {
    public:
      /// See EventSystem.
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

      /**
       * @brief Register callbacks for TapeStartRecording events.
       *
       * See notifyTapeStartRecordingListeners for the callback parameters and event order.
       *
       * @param callback    Callback to be invoked.
       * @param customData  Optional. Custom data that should be linked with the callback, otherwise nullptr.
       */
      static CODI_INLINE void registerTapeStartRecordingListener(void (*callback)(Tape&, void*),
                                                                 void* customData = nullptr) {
        Base::template internalRegisterListener(Config::ADWorkflowEvents, Event::TapeStartRecording, callback,
                                                customData);
      }

      /**
       * @brief Invoke callbacks for TapeStartRecording events.
       *
       * A TapeStartRecording event is triggered whenever a tape of the associated type is about to be set active,
       * that is, after the setActive call but before the internal status change.
       *
       * @param tape  Reference to the tape.
       */
      static CODI_INLINE void notifyTapeStartRecordingListeners(Tape& tape) {
        Base::template internalNotifyListeners<void (*)(Tape&, void*)>(Config::ADWorkflowEvents,
                                                                       Event::TapeStartRecording, tape);
      }

      /**
       * @brief Register callbacks for TapeStopRecording events.
       *
       * See notifyTapeStopRecordingListeners for the callback parameters and event order.
       *
       * @param callback    Callback to be invoked.
       * @param customData  Optional. Custom data that should be linked with the callback, otherwise nullptr.
       */
      static CODI_INLINE void registerTapeStopRecordingListener(void (*callback)(Tape&, void*),
                                                                void* customData = nullptr) {
        Base::template internalRegisterListener(Config::ADWorkflowEvents, Event::TapeStopRecording, callback,
                                                customData);
      }

      /**
       * @brief Invoke callbacks for TapeStopRecording events.
       *
       * A TapeStopRecording event is triggered whenever a tape of the associated type is about to be set passive,
       * that is, after the setPassive call but before the internal status change.
       *
       * @param tape  Reference to the tape.
       */
      static CODI_INLINE void notifyTapeStopRecordingListeners(Tape& tape) {
        Base::template internalNotifyListeners<void (*)(Tape&, void*)>(Config::ADWorkflowEvents,
                                                                       Event::TapeStopRecording, tape);
      }

      /**
       * @brief Register callbacks for TapeRegisterInput events.
       *
       * See notifyTapeRegisterInputListeners for the callback parameters and event order.
       *
       * @param callback    Callback to be invoked.
       * @param customData  Optional. Custom data that should be linked with the callback, otherwise nullptr.
       */
      static CODI_INLINE void registerTapeRegisterInputListener(void (*callback)(Tape&, Real&, Identifier&, void*),
                                                                void* customData = nullptr) {
        Base::template internalRegisterListener(Config::ADWorkflowEvents, Event::TapeRegisterInput, callback,
                                                customData);
      }

      /**
       * @brief Invoke callbacks for TapeRegisterInput events.
       *
       * A registerInput event occurs whenever registerInput calls are made to the tape, after the internal input
       * registration is completed.
       *
       * @param tape        Reference to the tape.
       * @param value       The value of the input that has been registered.
       * @param identifier  The new identifier of the input that has been registered.
       */
      static CODI_INLINE void notifyTapeRegisterInputListeners(Tape& tape, Real& value, Identifier& identifier) {
        Base::template internalNotifyListeners<void (*)(Tape&, Real&, Identifier&, void*)>(
            Config::ADWorkflowEvents, Event::TapeRegisterInput, tape, value, identifier);
      }

      /**
       * @brief Register callbacks for TapeRegisterOutput events.
       *
       * See notifyTapeRegisterOutputListeners for the callback parameters and event order.
       *
       * @param callback    Callback to be invoked.
       * @param customData  Optional. Custom data that should be linked with the callback, otherwise nullptr.
       */
      static CODI_INLINE void registerTapeRegisterOutputListener(void (*callback)(Tape&, Real&, Identifier&, void*),
                                                                 void* customData = nullptr) {
        Base::template internalRegisterListener(Config::ADWorkflowEvents, Event::TapeRegisterOutput, callback,
                                                customData);
      }

      /**
       * @brief Invoke callbacks for TapeRegisterOutput events.
       *
       * A registerInput event occurs whenever registerOutput calls are made to the tape, after the internal output
       * registration is completed.
       *
       * @param tape        Reference to the tape.
       * @param value       The value of the output that has been registered.
       * @param identifier  The new identifier of the output that has been registered.
       */
      static CODI_INLINE void notifyTapeRegisterOutputListeners(Tape& tape, Real& value, Identifier& identifier) {
        Base::template internalNotifyListeners<void (*)(Tape&, Real&, Identifier&, void*)>(
            Config::ADWorkflowEvents, Event::TapeRegisterOutput, tape, value, identifier);
      }

      /**
       * @brief Register callbacks for TapeEvaluate events.
       *
       * See notifyTapeEvaluateListeners for the callback parameters and event order.
       *
       * @param callback    Callback to be invoked.
       * @param customData  Optional. Custom data that should be linked with the callback, otherwise nullptr.
       */
      static CODI_INLINE void registerTapeEvaluateListener(void (*callback)(Tape&, Position const&, Position const&,
                                                                            VectorAccess*, EventHints::EvaluationKind,
                                                                            EventHints::Endpoint, void*),
                                                           void* customData = nullptr) {
        Base::template internalRegisterListener(Config::ADWorkflowEvents, Event::TapeEvaluate, callback, customData);
      }

      /**
       * @brief Invoke callbacks for TapeRegisterOutput events.
       *
       * Since a tape evaluation is not a singular event but a process, TapeEvaluate events occur both prior to and
       * after it, with correspondingly different endpoint parameter. The type of tape evaluation is indicated by the
       * evalKind parameter.
       *
       * @param tape       Reference to the tape.
       * @param start      Starting position of the evaluation.
       * @param end        End position of the evaluation.
       * @param adjoint    Vector access interface that provides access to adjoint variables.
       * @param evalKind   Indicates whether the evaluation is primal, forward, or reverse.
       * @param endpoint   Indicates whether the event occurs before or after the tape evaluation.
       */
      static CODI_INLINE void notifyTapeEvaluateListeners(Tape& tape, Position const& start, Position const& end,
                                                          VectorAccess* adjoint, EventHints::EvaluationKind evalKind,
                                                          EventHints::Endpoint endpoint) {
        Base::template internalNotifyListeners<void (*)(Tape&, Position const&, Position const&, VectorAccess*,
                                                        EventHints::EvaluationKind, EventHints::Endpoint, void*)>(
            Config::ADWorkflowEvents, Event::TapeEvaluate, tape, start, end, adjoint, evalKind, endpoint);
      }

      /**
       * @brief Register callbacks for TapeReset events.
       *
       * See notifyTapeResetListeners for the callback parameters and event order.
       *
       * @param callback    Callback to be invoked.
       * @param customData  Optional. Custom data that should be linked with the callback, otherwise nullptr.
       */
      static CODI_INLINE void registerTapeResetListener(void (*callback)(Tape&, Position const&, EventHints::Reset,
                                                                         bool, void*),
                                                        void* customData = nullptr) {
        Base::template internalRegisterListener(Config::ADWorkflowEvents, Event::TapeReset, callback, customData);
      }

      /**
       * @brief Invoke callbacks for TapeReset events.
       *
       * A TapeReset event occurs in the course of reset, resetTo, and resetHard calls made to a tape, the respective
       * origin indicated by kind. The event occurs before any internal reset takes place.
       *
       * @param tape           Reference to the tape.
       * @param position       Position to which we reset, zero position for full and hard resets.
       * @param kind           Reset kind, see codi::EventHints::Reset.
       * @param clearAdjoints  Indicates whether adjoints are zeroed.
       */
      static CODI_INLINE void notifyTapeResetListeners(Tape& tape, Position const& position, EventHints::Reset kind,
                                                       bool clearAdjoints) {
        Base::template internalNotifyListeners<void (*)(Tape&, Position const&, EventHints::Reset, bool, void*)>(
            Config::ADWorkflowEvents, Event::TapeReset, tape, position, kind, clearAdjoints);
      }

      /// @}
      /*******************************************************************************/
      /// @name Preaccumulation
      /// @{

      /**
       * @brief Register callbacks for PreaccStart events.
       *
       * See notifyPreaccFinishListeners for the callback parameters and event order.
       *
       * @param callback    Callback to be invoked.
       * @param customData  Optional. Custom data that should be linked with the callback, otherwise nullptr.
       */
      static CODI_INLINE void registerPreaccStartListener(void (*callback)(Tape&, void*), void* customData = nullptr) {
        Base::template internalRegisterListener(Config::PreaccEvents, Event::PreaccStart, callback, customData);
      }

      /**
       * @brief Invoke callbacks for PreaccStart events.
       *
       * A PreaccStart event occurs whenever a preaccumulation is started via a PreaccumulationHelper, prior to internal
       * preparations for the preaccumulation.
       *
       * @param tape  Reference to the tape.
       */
      static CODI_INLINE void notifyPreaccStartListeners(Tape& tape) {
        Base::template internalNotifyListeners<void (*)(Tape&, void*)>(Config::PreaccEvents, Event::PreaccStart, tape);
      }

      /**
       * @brief Register callbacks for PreaccFinish events.
       *
       * See notifyPreaccFinishListeners for the callback parameters and event order.
       *
       * @param callback    Callback to be invoked.
       * @param customData  Optional. Custom data that should be linked with the callback, otherwise nullptr.
       */
      static CODI_INLINE void registerPreaccFinishListener(void (*callback)(Tape&, void*), void* customData = nullptr) {
        Base::template internalRegisterListener(Config::PreaccEvents, Event::PreaccFinish, callback, customData);
      }

      /**
       * @brief Invoke callbacks for PreaccFinish events.
       *
       * A PreaccFinish event occurs whenever a preaccumulation is finished via a PreaccumulationHelper, after
       * everything related to the preaccumulation is completed. Note that there may be TapeEvaluate, TapeReset, and
       * StatementStore events due to preaccumulation after the call to finish but before this event.
       *
       * @param tape  Reference to the tape.
       */
      static CODI_INLINE void notifyPreaccFinishListeners(Tape& tape) {
        Base::template internalNotifyListeners<void (*)(Tape&, void*)>(Config::PreaccEvents, Event::PreaccFinish, tape);
      }

      /**
       * @brief Register callbacks for PreaccAddInput events.
       *
       * See notifyPreaccAddInputListeners for the callback parameters and event order.
       *
       * @param callback    Callback to be invoked.
       * @param customData  Optional. Custom data that should be linked with the callback, otherwise nullptr.
       */
      static CODI_INLINE void registerPreaccAddInputListener(void (*callback)(Tape&, Real const&, Identifier const&,
                                                                              void*),
                                                             void* customData = nullptr) {
        Base::template internalRegisterListener(Config::PreaccEvents, Event::PreaccAddInput, callback, customData);
      }

      /**
       * @brief Invoke callbacks for PreaccAddInput events.
       *
       * A PreaccAddInput event occurs when an input is added via a PreaccumulationHelper, either by passing them to
       * start or by calling addInput. The event occurs before the internal registration.
       *
       * @param tape        Reference to the tape.
       * @param value       The value of the input that is registered.
       * @param identifier  The identifier of the input that is registered.
       */
      static CODI_INLINE void notifyPreaccAddInputListeners(Tape& tape, Real const& value,
                                                            Identifier const& identifier) {
        Base::template internalNotifyListeners<void (*)(Tape&, Real const&, Identifier const&, void*)>(
            Config::PreaccEvents, Event::PreaccAddInput, tape, value, identifier);
      }

      /**
       * @brief Register callbacks for PreaccAddOutput events.
       *
       * See notifyPreaccAddOutputListeners for the callback parameters and event order.
       *
       * @param callback    Callback to be invoked.
       * @param customData  Optional. Custom data that should be linked with the callback, otherwise nullptr.
       */
      static CODI_INLINE void registerPreaccAddOutputListener(void (*callback)(Tape&, Real&, Identifier&, void*),
                                                              void* customData = nullptr) {
        Base::template internalRegisterListener(Config::PreaccEvents, Event::PreaccAddOutput, callback, customData);
      }

      /**
       * @brief Invoke callbacks for PreaccAddOutput events.
       *
       * A PreaccAddOutput event occurs when an output is added via a PreaccumulationHelper, either by passing them to
       * finish or by calling addOutput. The event occurs before the internal registration.
       *
       * @param tape        Reference to the tape.
       * @param value       The value of the output that is registered.
       * @param identifier  The identifier of the output that is registered.
       */
      static CODI_INLINE void notifyPreaccAddOutputListeners(Tape& tape, Real& value, Identifier& identifier) {
        Base::template internalNotifyListeners<void (*)(Tape&, Real&, Identifier&, void*)>(
            Config::PreaccEvents, Event::PreaccAddOutput, tape, value, identifier);
      }

      /// @}
      /*******************************************************************************/
      /// @name Statements
      /// @{

      /**
       * @brief Register callbacks for StatementStoreOnTape events.
       *
       * See notifyStatementStoreOnTapeListeners for the callback parameters and event order.
       *
       * @param callback    Callback to be invoked.
       * @param customData  Optional. Custom data that should be linked with the callback, otherwise nullptr.
       */
      static CODI_INLINE void registerStatementStoreOnTapeListener(
          void (*callback)(Tape&, Identifier const&, Real const&, size_t, Identifier const*, Real const*, void*),
          void* customData = nullptr) {
        Base::template internalRegisterListener(Config::StatementEvents, Event::StatementStoreOnTape, callback,
                                                customData);
      }

      /**
       * @brief Invoke callbacks for StatementStoreOnTape events.
       *
       * A StatementStoreOnTape event occurs whenever a statement is stored on the tape. Note that statements might not
       * be stored on the tape, for example because nothing is active, but still emit StatementPrimal events. The event
       * is emitted after storing all data on the tape, after the left hand side has been assigned an identifier but
       * before the left hand side updates its primal value.
       *
       * @param tape                Reference to the tape.
       * @param lhsIdentifier       Identifier assigned to the left hand side in the course of the statement.
       * @param newValue            Value that is assigned to the left hand side.
       * @param numActiveVariables  Number of active variables on the right hand side.
       * @param rhsIdentifiers      Pointer to numActiveVariables identifiers.
       * @param jacobians           Pointer to numActiveVariables Jacobian values.
       */
      static CODI_INLINE void notifyStatementStoreOnTapeListeners(Tape& tape, Identifier const& lhsIdentifier,
                                                                  Real const& newValue, size_t numActiveVariables,
                                                                  Identifier const* rhsIdentifiers,
                                                                  Real const* jacobians) {
        Base::template internalNotifyListeners<void (*)(Tape&, Identifier const&, Real const&, size_t,
                                                        Identifier const*, Real const*, void*)>(
            Config::StatementEvents, Event::StatementStoreOnTape, tape, lhsIdentifier, newValue, numActiveVariables,
            rhsIdentifiers, jacobians);
      }

      /**
       * @brief Register callbacks for StatementEvaluate events.
       *
       * See notifyStatementEvaluateListeners for the callback parameters and event order.
       *
       * @param callback    Callback to be invoked.
       * @param customData  Optional. Custom data that should be linked with the callback, otherwise nullptr.
       */
      static CODI_INLINE void registerStatementEvaluateListener(void (*callback)(Tape&, Identifier const&, size_t,
                                                                                 Real const*, void*),
                                                                void* customData = nullptr) {
        Base::template internalRegisterListener<void (*)(Tape&, Identifier const&, size_t, Real const*, void*)>(
            Config::StatementEvents, Event::StatementEvaluate, callback, customData);
      }

      /**
       * @brief Invoke callbacks for StatementEvaluate events.
       *
       * A StatementEvaluate occurs whenever a tape evaluates a statement.
       *
       * @param tape            Reference to the tape.
       * @param lhsIdentifier   Left hand side identifier.
       * @param sizeLhsAdjoint  Number of left hand side adjoint components.
       * @param lhsAdjoint      Pointer to the left hand side adjoint components.
       */
      static CODI_INLINE void notifyStatementEvaluateListeners(Tape& tape, Identifier const& lhsIdentifier,
                                                               size_t sizeLhsAdjoint, Real const* lhsAdjoint) {
        Base::template internalNotifyListeners<void (*)(Tape&, Identifier const&, size_t, Real const*, void*)>(
            Config::StatementEvents, Event::StatementEvaluate, tape, lhsIdentifier, sizeLhsAdjoint, lhsAdjoint);
      }

      /**
       * @brief Register callbacks for TapeStopRecording events.
       *
       * See notifyTapeStopRecording for the callback parameters and event order.
       *
       * @param callback    Callback to be invoked.
       * @param customData  Optional. Custom data that should be linked with the callback, otherwise nullptr.
       */
      static CODI_INLINE void registerStatementEvaluatePrimalListener(void (*callback)(Tape&, Identifier const&,
                                                                                       Real const&, void*),
                                                                      void* customData = nullptr) {
        Base::template internalRegisterListener<void (*)(Tape&, Identifier const&, Real const&, void*)>(
            Config::StatementEvents, Event::StatementEvaluatePrimal, callback, customData);
      }

      static CODI_INLINE void notifyStatementEvaluatePrimalListeners(Tape& tape, Identifier const& lhsIdentifier,
                                                                     Real const& lhsValue) {
        Base::template internalNotifyListeners<void (*)(Tape&, Identifier const&, Real const&, void*)>(
            Config::StatementEvents, Event::StatementEvaluatePrimal, tape, lhsIdentifier, lhsValue);
      }

      /// @}
      /*******************************************************************************/
      /// @name Index handling
      /// @{

      /**
       * @brief Register callbacks for TapeStopRecording events.
       *
       * See notifyTapeStopRecording for the callback parameters and event order.
       *
       * @param callback    Callback to be invoked.
       * @param customData  Optional. Custom data that should be linked with the callback, otherwise nullptr.
       */
      static CODI_INLINE void registerIndexAssignListener(void (*callback)(Index const&, void*),
                                                          void* customData = nullptr) {
        Base::template internalRegisterListener(Config::IndexEvents, Event::IndexAssign, callback, customData);
      }

      static CODI_INLINE void notifyIndexAssignListeners(Index const& index) {
        Base::template internalNotifyListeners<void (*)(Index const&, void*)>(Config::IndexEvents, Event::IndexAssign,
                                                                              index);
      }

      /**
       * @brief Register callbacks for TapeStopRecording events.
       *
       * See notifyTapeStopRecording for the callback parameters and event order.
       *
       * @param callback    Callback to be invoked.
       * @param customData  Optional. Custom data that should be linked with the callback, otherwise nullptr.
       */
      static CODI_INLINE void registerIndexFreeListener(void (*callback)(Index const&, void*),
                                                        void* customData = nullptr) {
        Base::template internalRegisterListener(Config::IndexEvents, Event::IndexFree, callback, customData);
      }

      static CODI_INLINE void notifyIndexFreeListeners(Index const& index) {
        Base::template internalNotifyListeners<void (*)(Index const&, void*)>(Config::IndexEvents, Event::IndexFree,
                                                                              index);
      }

      /**
       * @brief Register callbacks for TapeStopRecording events.
       *
       * See notifyTapeStopRecording for the callback parameters and event order.
       *
       * @param callback    Callback to be invoked.
       * @param customData  Optional. Custom data that should be linked with the callback, otherwise nullptr.
       */
      static CODI_INLINE void registerIndexCopyListener(void (*callback)(Index const&, void*),
                                                        void* customData = nullptr) {
        Base::template internalRegisterListener(Config::IndexEvents, Event::IndexCopy, callback, customData);
      }

      static CODI_INLINE void notifyIndexCopyListeners(Index const& index) {
        Base::template internalNotifyListeners<void (*)(Index const&, void*)>(Config::IndexEvents, Event::IndexCopy,
                                                                              index);
      }

      /// @}
  };

  /* specialization for the forward evaluation "tape" is identical to EventSystemBase */

  template<typename Real, typename Gradient>
  struct ForwardEvaluation;

  template<typename Real, typename Gradient>
  struct EventSystem<ForwardEvaluation<Real, Gradient>> : public EventSystemBase<ForwardEvaluation<Real, Gradient>> {};

}
