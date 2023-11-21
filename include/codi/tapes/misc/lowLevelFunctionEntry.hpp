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
  enum class LowLevelFunctionEntryCallType {
    Forward,
    Reverse,
    Primal,
    Delete,
    Count,
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
      using FuncEval = void (*)(Tape* tape, ByteDataView& fixedData, ByteDataView& dynamicData,
                                VectorAccessInterface<Real, Identifier>* access);

      /// Call syntax for Delete calls.
      using FuncDel = void (*)(Tape* tape, ByteDataView& fixedData, ByteDataView& dynamicData);

      /// @brief Call syntax for Count calls.
      ///
      /// \c fixedSize returns the allocated fixed size on the tape, dynamicSize the allocated dynamic size on the tape
      /// and \c allocatedSize the memory allocated on the system (e.g. with new).
      using FuncCount = void (*)(Tape* tape, ByteDataView& fixedData, ByteDataView& dynamicData, int& fixedSize,
                                 int& dynamicSize, int& allocatedSize);

    private:

      void* functions[(size_t)LowLevelFunctionEntryCallType::MaxElement];  ///< Array for function pointers.
      using FunctionTypes =
          std::tuple<FuncEval, FuncEval, FuncEval, FuncDel, FuncCount>;  ///< Types for function entries.

    public:

      /// Constructors.
      LowLevelFunctionEntry(FuncEval reverse = nullptr, FuncEval forward = nullptr, FuncEval primal = nullptr,
                            FuncDel del = nullptr, FuncCount count = nullptr)
          : functions{(void*)forward, (void*)reverse, (void*)primal, (void*)del, (void*)count} {}

      /// Call the function corresponding to callType with the given arguments.
      template<LowLevelFunctionEntryCallType callType, typename... Args>
      void call(Args&&... args) const {
        using FuncType = typename std::tuple_element<(size_t)callType, FunctionTypes>::type;
        ((FuncType)functions[(size_t)callType])(std::forward<Args>(args)...);
      }

      /// Check if a function is provided for the callType.
      template<LowLevelFunctionEntryCallType callType, typename... Args>
      bool has(Args&&... args) const {
        return nullptr != functions[(size_t)callType];
      }
  };

}
