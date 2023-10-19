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
#include "../../misc/byteDataStore.hpp"
#include "../../misc/macros.hpp"
#include "../../misc/temporaryMemoryAllocator.hpp"
#include "../misc/lowLevelFunctionEntry.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename T_Tape>
  struct ExternalFunction;

  template<typename T_Real, typename T_Gradient, typename T_Tape, typename T_Impl>
  struct LhsExpressionInterface;

  /**
   * @brief Add small recurring function to the tape evaluation.
   *
   * A low level function like the matrix matrix multiplication can be add with this interface. First the function needs
   * to be registered with #registerLowLevelFunction. This needs to be done only once. Afterwards a low level function
   * can be pushed as often as required with #pushLowLevelFunction.
   *
   * The user can write arbitrary data into the fixed and dynamic data streams. The only requirement is that at the
   * start and end of the fixed data stream the token created with #pushLowLevelFunction needs to be written.
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

      /// An allocator for dynamic data during the evaluation on the tape or during the recording.
      TemporaryMemoryAllocator& getTemporaryMemoryAllocator();

      /**
       *  @brief Internal function for evaluating a low level function entry.
       *
       *  @param tape                    The tape that is evaluated.
       *  @param direction               The read direction for the data pointers.
       *  @param curOtherDynamicDataPos  Position of dynamic data.
       *  @param otherDynamicPtr         Pointer for dynamic data.
       *  @param curOtherFixedDataPos    Position of fixed data.
       *  @param otherFixedPtr           Pointer for fixed data.
       *  @param args                    Additional arguments for the function call.
       *  @tparam callType               The function type that is called.
       */
      template<LowLevelFunctionEntryCallType callType, typename... Args>
      static void handleLowLevelFunction(LowLevelFunctionTapeInterface& tape, ByteDataStore::Direction direction,
                                         size_t& curOtherDynamicDataPos, char const* const otherDynamicPtr,
                                         size_t& curOtherFixedDataPos, char const* const otherFixedPtr, Args&&... args);

      /**
       *  Push a low level function to the tape.
       *
       *  \c fixedSize and \c dynamicSize are allocated on the vectors. \c fixedData and \c dynamicData are initialized
       *  with the allocated data. After the call, the data for the function can be written into the data stores.
       *
       *  See LowLevelFunctionTapeInterface for the expected data layout.
       */
      void pushLowLevelFunction(size_t fixedSize, size_t dynamicSize, ByteDataStore& fixedData,
                                ByteDataStore& dynamicData);

      /// Register a low level function on the tape.
      Config::LowLevelFunctionToken registerLowLevelFunction(
          LowLevelFunctionEntry<LowLevelFunctionTapeInterface, Real, Identifier> const& entry);
  };
}
