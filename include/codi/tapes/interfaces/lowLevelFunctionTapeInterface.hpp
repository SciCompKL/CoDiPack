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

#include "../../config.h"
#include "../../misc/byteDataView.hpp"
#include "../../misc/macros.hpp"
#include "../../misc/temporaryMemory.hpp"
#include "../misc/lowLevelFunctionEntry.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename T_Tape>
  struct ExternalFunction;

  template<typename T_Real, typename T_Gradient, typename T_Tape, typename T_Impl>
  struct LhsExpressionInterface;

  /**
   * @brief Add functions with custom derivatives to the tape. Can, e.g, be used to optimize small recurring functions like matrix matrix multiplication.
   *
   * A low level function like the matrix matrix multiplication can be added with this interface. First, the function needs
   * to be registered with #registerLowLevelFunction. This needs to be done only once, after this, the function is generally
   * available. It can then be pushed as often as required with #pushLowLevelFunction. Each push can be accompanied by different
   * data, e.g., the specific matrices used by individual matrix matrix multiplications.
   *
   * The user can write arbitrary data into the fixed and dynamic data streams. There is no requirement on the layout,
   * but the data should be readable from the left and right. Therefore, the fixed data stream is used for data that is
   * always read and dynamic data is used for data which depends on the data in the fixed data stream.
   *
   * @tparam T_Real        The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_Gradient    The gradient type of a tape, usually chosen as ActiveType::Gradient.
   * @tparam T_Identifier  The adjoint/tangent identification type of a tape, usually chosen as ActiveType::Identifier.
   */
  template<typename T_Real, typename T_Gradient, typename T_Identifier>
  struct LowLevelFunctionTapeInterface {
    public:

      using Real = CODI_DD(T_Real, double);           ///< See LowLevelFunctionTapeInterface.
      using Gradient = CODI_DD(T_Gradient, double);   ///< See LowLevelFunctionTapeInterface.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See LowLevelFunctionTapeInterface.

      /*******************************************************************************/
      /// @name Interface definition

      /// Temporary memory that can be used for dynamic data both during the evaluation and the recording.
      TemporaryMemory& getTemporaryMemory();

      /**
       *  @brief Internal function for evaluating a low level function entry.
       *
       *  The positions \c curDynamicDataPos,curFixedDataPos and \c curLLFTokenDataPos are advanced according to
       *  \c direction.
       *
       *  @param tape                    The tape that is evaluated.
       *  @param direction               The read direction for the data pointers.
       *  @param curDynamicDataPos       Position of dynamic data.
       *  @param dynamicDataPtr          Pointer for dynamic data.
       *  @param curFixedDataPos         Position of fixed data.
       *  @param fixedDataPtr            Pointer for fixed data.
       *  @param curLLFTokenDataPos      Position of token data.
       *  @param tokenPtr                Pointer for token data.
       *  @param args                    Additional arguments for the function call.
       *  @tparam callType               The function type that is called.
       */
      template<LowLevelFunctionEntryCallType callType, typename... Args>
      static void callLowLevelFunction(LowLevelFunctionTapeInterface& tape, ByteDataView::Direction direction,
                                         size_t& curDynamicDataPos, char* dynamicDataPtr,
                                         size_t& curFixedDataPos, char* fixedDataPtr,
                                         size_t& curLLFTokenDataPos, Config::LowLevelFunctionToken* const tokenPtr,
                                         Args&&... args);

      /**
       *  @brief Push a low level function to the tape.
       *
       *  Allocates memory with the requested sizes \c fixedSize and \c dynamicSize on the respective data streams.
       *  \c fixedData and \c dynamicData are initialized for accessing this allocated memory.
       *  After the call, they can be used to write data to the data streams.
       *  \c token is the token from #registerLowLevelFunction.
       *
       *  See LowLevelFunctionTapeInterface for the expected data layout.
       */
      void pushLowLevelFunction(Config::LowLevelFunctionToken token, size_t fixedSize, size_t dynamicSize,
                                ByteDataView& fixedData, ByteDataView& dynamicData);

      /// Register a low level function on the tape.
      Config::LowLevelFunctionToken registerLowLevelFunction(
          LowLevelFunctionEntry<LowLevelFunctionTapeInterface, Real, Identifier> const& entry);
  };
}
