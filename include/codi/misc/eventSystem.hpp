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

  /* might need to be moved to an extra file */
  enum class Event {
    /* high-level events */
    StartRecording,
    StopRecording,
    /* ... */
    /* low-level events */
    /*... */
    /* index management events */
    IndexAssign
    /* ... */
  };

  /* specialize if signature is different, and with respect to "enabled", depending on how it will be toggled */
  template<Event T_event, typename T_Tape>
  struct EventTraits {
    public:
      static Event constexpr event = CODI_DD(T_event, Event::StartRecording);
      using Tape = CODI_DD(T_Tape, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));
      static bool constexpr enabled = Config::HighLevelEvents;

      using Callback = void (*)(Tape&, void*);
  };

  /* no need to specialize for StopRecording right now */

  template<typename T_Tape>
  struct EventTraits<Event::IndexAssign, T_Tape> {
    public:
      static Event constexpr event = Event::IndexAssign;
      using Tape = CODI_DD(T_Tape, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));
      static bool constexpr enabled = Config::LowLevelEvents;

      using Index = typename Tape::Identifier;

      using Callback = void (*)(Index&, void*);
  };

  template<typename T_Tape>
  struct EventSystem {
    public:
      using Tape = CODI_DD(T_Tape, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));

      using Callback = void*;

      template<Event event>
      using EventTraits = EventTraits<event, Tape>;

    private:
      static std::map<Event, std::list<std::pair<Callback, void*>>> listeners;

    public:
      template<Event event>
      static CODI_INLINE void registerListener(typename EventTraits<event>::Callback callback, void* customData = nullptr) {
        listeners[event].push_back(std::make_pair((void*)callback, customData));
      }

      template<Event event, typename... Args>
      static CODI_INLINE void notifyListeners(Args&&... args) {
        if (EventTraits<event>::enabled) {
          for (auto const& listener : listeners[event]) {
            ((typename EventTraits<event>::Callback) listener.first)(std::forward<Args...>(args...), listener.second);
          }
        }
      }
  };

  template<typename Tape>
  std::map<Event, std::list<std::pair<typename EventSystem<Tape>::Callback, void*>>> EventSystem<Tape>::listeners;

  template<Event event, typename ActiveType>
  void registerListener(typename EventTraits<event, typename ActiveType::Tape>::Callback callback, void* customData = nullptr) {
    EventSystem<typename ActiveType::Tape>::template registerListener<event>(callback, customData);
  }
}
