/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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

#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "../../config.h"
#include "../../misc/byteDataView.hpp"
#include "../../misc/macros.hpp"
#include "vectorAccessInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// All possible call types for a low level function entry.
  enum class LowLevelFunctionEntryCallKind {
    Forward,
    Reverse,
    Primal,
    Delete,
    IterateInputs,
    IterateOutputs,
    MaxElement
  };

  /**
   * @brief Low level function entry on the tape. See LowLevelFunctionTapeInterface for details.
   *
   * @tparam T_Tape        The tape on which the entry is registered.
   * @tparam T_Real        The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_Identifier  The adjoint/tangent identification type of a tape, usually chosen as ActiveType::Identifier.
   */
  template<typename T_Tape, typename T_Real, typename T_Identifier>
  struct LowLevelFunctionEntry {
      using Tape = T_Tape;                            ///< See LowLevelFunctionEntry.
      using Real = CODI_DD(T_Real, double);           ///< See LowLevelFunctionEntry.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See LowLevelFunctionEntry.

      /// Call syntax for Forward, Reverse, and Primal calls.
      using FuncEval = void (*)(Tape* tape, ByteDataView& data, VectorAccessInterface<Real, Identifier>* access);

      /// Call syntax for Delete calls.
      using FuncDel = void (*)(Tape* tape, ByteDataView& data);

      /// Callback function for the identifier iteration.
      using IterCallback = void (*)(Identifier* id, void* userData);
      /// Call syntax for IterateInputs and IterateOutputs calls.
      using FuncIterate = void (*)(Tape* tape, ByteDataView& data, IterCallback func, void* userData);

    private:

      void* functions[(size_t)LowLevelFunctionEntryCallKind::MaxElement];  ///< Array for function pointers.
      using FunctionTypes =
          std::tuple<FuncEval, FuncEval, FuncEval, FuncDel, FuncIterate, FuncIterate>;  ///< Types for function entries.

    public:

      /// Constructors.
      LowLevelFunctionEntry(FuncEval reverse = nullptr, FuncEval forward = nullptr, FuncEval primal = nullptr,
                            FuncDel del = nullptr, FuncIterate iterIn = nullptr, FuncIterate iterOut = nullptr)
          : functions{(void*)forward, (void*)reverse, (void*)primal, (void*)del, (void*)iterIn, (void*)iterOut} {}

      /// Call the function corresponding to callType with the given arguments.
      template<LowLevelFunctionEntryCallKind callType, typename... Args>
      void call(Args&&... args) const {
        using FuncType = typename std::tuple_element<(size_t)callType, FunctionTypes>::type;
        ((FuncType)functions[(size_t)callType])(std::forward<Args>(args)...);
      }

      /// Check if a function is provided for the callType.
      template<LowLevelFunctionEntryCallKind callType, typename... Args>
      bool has(Args&&... args) const {
        CODI_UNUSED(args...);
        return nullptr != functions[(size_t)callType];
      }
  };

}
