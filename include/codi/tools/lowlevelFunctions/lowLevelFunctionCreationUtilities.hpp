/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2024 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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

#include <tuple>

#include "../../config.h"
#include "../../misc/byteDataView.hpp"
#include "../../misc/macros.hpp"
#include "storeAndRestoreActions.hpp"
#include "traits/activeArgumentStoreTraits.hpp"
#include "traits/passiveArgumentStoreTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Helper structure for storing low level functions and their arguments on a tape.
   *
   * \section storing Storing on the tape
   *
   * The process of storing the data for a low level function consists of several steps:
   *  - Determine if the tape is active and if the arguments of the function are active.
   *  - Count the size of the data that needs to be stored.
   *  - Allocate the data on the tape.
   *  - Write the data to the tape.
   *  - Evaluate the function either in a passive taping context or with different types than the CoDiPack ones.
   *  - Register the active outputs of the function.
   *
   *  In the following sections, the steps are explained. For an example implementation, see
   *  #codi::ExtFunc_matrixMatrixMultiplication.
   *
   *  \subsection activity Determine activity
   *
   *  If the tape is not active, then no low level function should be created. Otherwise, each active argument can be
   *  checked with \c ActiveStoreTrait::isActive(). If all active arguments are passive, that is, all underlying
   *  CoDiPack types are passive, then the low level function should also not be created.
   *
   *  Currently the number of low-level function arguments is limited to 64. Note that the limit applies to the number
   *  of arguments passed to the low-level function, which is in general much smaller than the number of active
   *  CoDiPack values passed to the function. Matrix-matrix multiplication, for example, has two arguments (two
   *  matrices), but each matrix consists of multiple CoDiPack values.
   *
   *  \subsection count Count size
   *
   *  The required size for storing all arguments can be determined with \c ActiveStoreTrait::countSize() and
   *  \c PassiveStoreTrait::countSize().
   *
   *  The total size is then:
   *  \code
   *  size = countActivitySize() + <size from all arguments>
   *  \endcode
   *
   *  Depending on the activity of the arguments and if the primal values of the arguments are required for the
   *  derivative computation, the required size may vary.
   *
   *  \subsection allocate Allocate data on the tape
   *
   *  A call to #LowLevelFunctionTapeInterface::pushLowLevelFunction() adds the low level function to the tape and
   *  populates the #codi::ByteDataView for the fixed and dynamic data.
   *
   *  \subsection write Write data
   *
   *  Usually, the following needs to be done:
   *   - Call #setActivity for every input argument.
   *   - Write the activity of the arguments with #storeActivity.
   *   - Call \c ActiveStoreTrait::store and \c PassiveArgumentStoreTraits::store for all arguments.
   *
   *  \subsection evaluate Evaluation of the low level function
   *
   *  There are two options.
   *   - The tape can be set to passive and the arguments with the CoDiPack types can be used for the evaluation.
   *   - The store traits can be configured such that the primal values are always extracted. These are available via
   *   \c ActiveStoreTrait::ArgumentStore::value() and can be used to call a passive version of the low level function.
   *
   *  \subsection register Register output arguments
   *
   *  After the low level function is evaluated, each active output argument needs to be registered on the tape with
   *  a call to \c ActiveStoreTrait::setExternalFunctionOutput().
   *
   *  \section restore Restoring for reverse and forward evaluation
   *
   *  The restoring process needs to read the data in the same order as it was written. The following needs to be done:
   *    - Read the activity of the arguments with #restoreActivity.
   *    - Read the data for all arguments with \c ActiveStoreTrait::restore() or \c PassiveStoreTrait::restore().
   *
   *  \section deletion Data deletion
   *
   *  Currently, it is assumed that all data is stored in the binary data stream of the tape. This data stream is
   *  automatically reset during a tape reset, which releases the data. There is currently no mechanism to automatically
   *  delete heap allocated data in the store methods.
   *
   * @tparam T_ActiveArguments The number of active input arguments.
   */
  template<size_t T_ActiveArguments>
  struct LowLevelFunctionCreationUtilities {
      static size_t constexpr ActiveArguments =
          CODI_DD(T_ActiveArguments, 1);  ///< See LowLevelFunctionCreationUtilities.

      static_assert(ActiveArguments <= 64, "More than 64 active arguments are currently not supported.");

      // clang-format off
      /// Type for the activity store. Currently limited to 64 variables.
      using ActivityStoreType =
          typename std::conditional<ActiveArguments <= 8, uint8_t,
          typename std::conditional<ActiveArguments <= 16, uint16_t,
          typename std::conditional<ActiveArguments <= 32, uint32_t,
          typename std::conditional<ActiveArguments <= 64, uint64_t,
            void>::type>::type>::type>::type;
      // clang-format on

      /// Utilities for handling an active argument.
      template<typename T>
      using ActiveStoreTrait = ActiveArgumentStoreTraits<T>;

      /// Functionality for handling a passive argument.
      template<typename T, typename S = T>
      using PassiveStoreTrait = PassiveArgumentStoreTraits<T, S>;

      /// Token size for the low level function token.
      static size_t constexpr TokenSize = sizeof(codi::Config::LowLevelFunctionToken);

      /*******************************************************************************/
      /// @name Action creation
      /// @{

      /// Collect the restore actions for an argument in a \c RestoreActions instance.
      CODI_INLINE static RestoreActions createRestoreActions(bool isInput, bool isOutput, bool isInputActive,
                                                             bool primalRequired) {
        RestoreActions actions = {};

        if (isInput && primalRequired) {
          actions |= RestoreAction::PrimalRestore;
        } else if (isOutput) {
          actions |= RestoreAction::PrimalCreate;
        }

        if (isInput && isInputActive) {
          actions |= RestoreAction::InputIdentifierRestore | RestoreAction::InputGradientCreate;
        }
        if (isOutput) {
          actions |= RestoreAction::OutputIdentifierRestore | RestoreAction::OutputGradientCreate;
        }

        return actions;
      }

      /// Collect the store actions for an argument in a \c StoreActions instance.
      CODI_INLINE static StoreActions createStoreActions(bool tapeActive, bool isInput, bool isOutput,
                                                         bool isInputActive, bool primalRequired) {
        StoreActions actions = {};

        if (tapeActive && isInput && primalRequired) {
          actions |= StoreAction::PrimalCreateOnTape;
        }

        if (isInput) {
          actions |= StoreAction::PrimalExtract;
        }

        if (tapeActive) {
          if (isInput && isInputActive) {
            actions |= StoreAction::InputIdentifierCreateAndStore;
          }
          if (isOutput) {
            actions |= StoreAction::OutputIdentifierCreate;
          }
        }

        return actions;
      }

      /// @}
      /*******************************************************************************/
      /// @name Argument activity
      /// @{

      /// Return the size of the activity structure.
      CODI_INLINE static size_t countActivitySize() {
        return sizeof(ActivityStoreType);
      }

      /// Check the activity structure for activity of a specific argument.
      CODI_INLINE static bool getActivity(ActivityStoreType const& activity, size_t arg) {
        return 0 != (activity & ((ActivityStoreType)1) << arg);
      }

      /// Restore the activity structure from the data stream.
      CODI_INLINE static void restoreActivity(ByteDataView* fixedStore, ActivityStoreType& activity) {
        activity = fixedStore->read<ActivityStoreType>();
      }

      /// Store the activity of an argument in the activity structure.
      CODI_INLINE static void setActivity(ActivityStoreType& activity, size_t arg, bool active) {
        activity |= ((ActivityStoreType)active) << arg;
      }

      /// Store the activity structure in the data stream.
      CODI_INLINE static void storeActivity(ByteDataView* fixedStore, ActivityStoreType const& activity) {
        fixedStore->write(activity);
      }

      /// @}
  };
}
