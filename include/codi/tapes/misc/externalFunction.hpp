/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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

#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../misc/exceptions.hpp"
#include "../../misc/macros.hpp"
#include "../data/position.hpp"
#include "../interfaces/externalFunctionTapeInterface.hpp"
#include "../misc/vectorAccessInterface.hpp"
#include "lowLevelFunctionEntry.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Internal untyped data for an external function.
  struct ExternalFunctionInternalData {
    protected:
      using CallFunctionUntyped = void (*)(void* tape, void* data,
                                           void* adjointInterface);    ///< Call function definition.
      using DeleteFunctionUntyped = void (*)(void* tape, void* data);  ///< Delete function definition.

      using IterCallbackUntyped = void (*)(void* id, void* userData);  ///< Untyped callback function for id iteration.
      using IterateIdsFunctionUntyped = void (*)(void* tape, void* data, IterCallbackUntyped callback,
                                                 void* userData);  ///< Iterate ids function definition.

      CallFunctionUntyped funcReverse;        ///< Reverse evaluation function pointer.
      CallFunctionUntyped funcForward;        ///< Forward evaluation function pointer.
      CallFunctionUntyped funcPrimal;         ///< Primal evaluation function pointer.
      DeleteFunctionUntyped funcDelete;       ///< User data deletion function pointer.
      IterateIdsFunctionUntyped funcIterIn;   ///< Iterate over inputs.
      IterateIdsFunctionUntyped funcIterOut;  ///< Iterate over outputs.

      void* data;  ///< User data pointer.

    public:

      /// Constructor
      ExternalFunctionInternalData()
          : funcReverse(nullptr),
            funcForward(nullptr),
            funcPrimal(nullptr),
            funcDelete(nullptr),
            funcIterIn(nullptr),
            funcIterOut(nullptr),
            data(nullptr) {}

    protected:

      /// Constructor
      ExternalFunctionInternalData(CallFunctionUntyped funcReverse, CallFunctionUntyped funcForward,
                                   CallFunctionUntyped funcPrimal, DeleteFunctionUntyped funcDelete,
                                   IterateIdsFunctionUntyped funcIterIn, IterateIdsFunctionUntyped funcIterOut,
                                   void* data)
          : funcReverse(funcReverse),
            funcForward(funcForward),
            funcPrimal(funcPrimal),
            funcDelete(funcDelete),
            funcIterIn(funcIterIn),
            funcIterOut(funcIterOut),
            data(data) {}
  };

  /**
   * @brief User-defined evaluation functions for the taping process.
   *
   * See ExternalFunctionTapeInterface for details.
   *
   * The user can provide call functions for the reverse, forward and primal evaluation of a tape. These need to be of
   * the type CallFunction which has the tree arguments:
   *  - tape: The type of the tape on which this object was registered with `registerExternalFunction`
   *          (ExternalFunctionTapeInterface)
   *  - data: User-provided data, type is known by the user.
   *  - adjointInterface: VectorAccessInterface instantiated with <typename Tape::Real, typename Tape::Identifier>
   *
   * The tape pointer can be used for general access to the tape. For each access to the gradient data, the
   * adjointInterface should be used. If no custom adjoint vectors are used in the application, then the tape pointer
   * can also be used for the gradient data access.
   *
   * The delete function is called when the entry of the tape for the external function is deleted.
   *
   * @tparam T_Tape        The associated tape type.
   */
  template<typename T_Tape>
  struct ExternalFunction : public ExternalFunctionInternalData {
    public:

      using Tape = CODI_DD(T_Tape,
                           CODI_T(ExternalFunctionTapeInterface<double, double, int>));  ///< See ExternalFunction
      using Real = typename Tape::Real;                                                  ///< See FullTapeInterface.
      using Identifier = typename Tape::Identifier;                                      ///< See FullTapeInterface.

      using VectorAccess = VectorAccessInterface<Real, Identifier>;  ///< Shortcut for
                                                                     ///< VectorAccessInterface.
      using IterCallback =
          typename LowLevelFunctionEntry<Tape, Real, Identifier>::IterCallback;  ///< Callback for iterate ids function.

      using CallFunction = void (*)(Tape* tape, void* data,
                                    VectorAccess* adjointInterface);  ///< Call function definition.
      using DeleteFunction = void (*)(Tape* tape, void* data);        ///< Delete function definition.
      using IterateIdsFunction = void (*)(Tape* tape, void* data, IterCallback callback,
                                          void* userData);  ///< Iterate ids function definition.

      /// Any arguments can be nullptr if not required.
      ExternalFunction(CallFunction funcReverse, CallFunction funcForward, CallFunction funcPrimal, void* data,
                       DeleteFunction funcDelete, IterateIdsFunction funcIterIn, IterateIdsFunction funcIterOut)
          : ExternalFunctionInternalData((ExternalFunctionInternalData::CallFunctionUntyped)funcReverse,
                                         (ExternalFunctionInternalData::CallFunctionUntyped)funcForward,
                                         (ExternalFunctionInternalData::CallFunctionUntyped)funcPrimal,
                                         (ExternalFunctionInternalData::DeleteFunctionUntyped)funcDelete,
                                         (ExternalFunctionInternalData::IterateIdsFunctionUntyped)funcIterIn,
                                         (ExternalFunctionInternalData::IterateIdsFunctionUntyped)funcIterOut, data) {}

      /// Helper function for the creation of an ExternalFunction object.
      static ExternalFunction create(CallFunction funcReverse, void* data, DeleteFunction funcDelete,
                                     CallFunction funcForward = nullptr, CallFunction funcPrimal = nullptr,
                                     IterateIdsFunction funcIterIn = nullptr,
                                     IterateIdsFunction funcIterOut = nullptr) {
        return ExternalFunction(funcReverse, funcForward, funcPrimal, data, funcDelete, funcIterIn, funcIterOut);
      }

      /// Calls the delete function if not nullptr.
      void deleteData(Tape* tape) {
        if (funcDelete != nullptr) {
          funcDelete(tape, data);
          data = nullptr;
        }
      }

      /// Calls the reverse function if not nullptr, otherwise throws a CODI_EXCEPTION.
      void evaluateReverse(Tape* tape, VectorAccess* adjointInterface) const {
        if (nullptr != funcReverse) {
          funcReverse(tape, data, adjointInterface);
        } else {
          CODI_EXCEPTION(
              "Calling an external function in reverse mode without providing a reverse evaluation function.");
        }
      }

      /// Calls the forward function if not nullptr, otherwise throws a CODI_EXCEPTION.
      void evaluateForward(Tape* tape, VectorAccess* adjointInterface) const {
        if (nullptr != funcForward) {
          funcForward(tape, data, adjointInterface);
        } else {
          CODI_EXCEPTION(
              "Calling an external function in forward mode without providing a forward evaluation function.");
        }
      }

      /// Calls the primal function if not nullptr, otherwise throws a CODI_EXCEPTION.
      void evaluatePrimal(Tape* tape, VectorAccess* adjointInterface) const {
        if (nullptr != funcPrimal) {
          funcPrimal(tape, data, adjointInterface);
        } else {
          CODI_EXCEPTION("Calling an external function in primal mode without providing a primal evaluation function.");
        }
      }

      /// Calls the iterate inputs function if not nullptr, otherwise throws a CODI_EXCEPTION.
      void iterateInputs(Tape* tape, IterCallback func, void* userData) const {
        if (nullptr != funcIterIn) {
          funcIterIn(tape, data, (ExternalFunctionInternalData::IterCallbackUntyped)func, userData);
        } else {
          CODI_EXCEPTION(
              "Calling an external function for iteration of inputs without providing an iteration function.");
        }
      }

      /// Calls the iterate inputs function if not nullptr, otherwise throws a CODI_EXCEPTION.
      void iterateOutputs(Tape* tape, IterCallback func, void* userData) const {
        if (nullptr != funcIterOut) {
          funcIterOut(tape, data, (ExternalFunctionInternalData::IterCallbackUntyped)func, userData);
        } else {
          CODI_EXCEPTION(
              "Calling an external function for iteration of outputs without providing an iteration function.");
        }
      }
  };

  /**
   *  @brief Low level function entry implementation for external functions.
   *
   *  Stores the ExternalFunction object in the byte data stream.
   *
   * @tparam T_Tape        The associated tape type.
   * @tparam T_Gradient    The gradient type of a tape, usually chosen as ActiveType::Gradient.
   * @tparam T_Identifier  The adjoint/tangent identification of a tape, usually chosen as ActiveType::Identifier.
   *
   */
  template<typename T_Tape, typename T_Real, typename T_Identifier>
  struct ExternalFunctionLowLevelEntryMapper {
      using Tape = CODI_DD(T_Tape, CODI_DEFAULT_TAPE);  ///< See ExternalFunctionLowLevelEntryMapper.
      using Real = CODI_DD(T_Real, double);             ///< See ExternalFunctionLowLevelEntryMapper.
      using Identifier = CODI_DD(T_Identifier, int);    ///< See ExternalFunctionLowLevelEntryMapper.

      /// See LowLevelFunctionEntry.
      using IterCallback = typename LowLevelFunctionEntry<T_Tape, T_Real, T_Identifier>::IterCallback;

      using ExtFunc = ExternalFunction<Tape>;  ///< Abbreviation for ExternalFunction.
      /// Shortcut for VectorAccessInterface.
      using VectorAccess = VectorAccessInterface<Real, Identifier>;

      /// Recovers the external function data and calls evaluateForward on it.
      CODI_INLINE static void forward(Tape* tape, ByteDataView& data, VectorAccess* access) {
        ExtFunc* extFunc = data.read<ExtFunc>(1);
        extFunc->evaluateForward(tape, access);
      }

      /// Recovers the external function data and calls evaluatePrimal on it.
      CODI_INLINE static void primal(Tape* tape, ByteDataView& data, VectorAccess* access) {
        ExtFunc* extFunc = data.read<ExtFunc>(1);
        extFunc->evaluatePrimal(tape, access);
      }

      /// Recovers the external function data and calls evaluateReverse on it.
      CODI_INLINE static void reverse(Tape* tape, ByteDataView& data, VectorAccess* access) {
        ExtFunc* extFunc = data.read<ExtFunc>(1);
        extFunc->evaluateReverse(tape, access);
      }

      /// Recovers the external function data and calls deleteData on it.
      CODI_INLINE static void del(Tape* tape, ByteDataView& data) {
        ExtFunc* extFunc = data.read<ExtFunc>(1);
        extFunc->deleteData(tape);
      }

      /// Iterate over the inputs of the external function. \c func is called for each input with \c userData.
      CODI_INLINE static void iterateInputs(Tape* tape, ByteDataView& data, IterCallback func, void* userData) {
        CODI_UNUSED(tape);

        ExtFunc* extFunc = data.read<ExtFunc>(1);
        extFunc->iterateInputs(tape, func, userData);
      }

      /// Iterate over the outputs of the external function. \c func is called for each output with \c userData.
      CODI_INLINE static void iterateOutputs(Tape* tape, ByteDataView& data, IterCallback func, void* userData) {
        CODI_UNUSED(tape);

        ExtFunc* extFunc = data.read<ExtFunc>(1);
        extFunc->iterateOutputs(tape, func, userData);
      }

      /// Store an external function on the tape.
      CODI_INLINE static void store(Tape& tape, Config::LowLevelFunctionToken token, ExtFunc const& extFunc) {
        ByteDataView data = {};
        tape.pushLowLevelFunction(token, sizeof(ExtFunc), data);

        data.write(extFunc);
      }

      /// Create the function entry for the tape registration.
      CODI_INLINE static LowLevelFunctionEntry<Tape, Real, Identifier> create() {
        return LowLevelFunctionEntry<Tape, Real, Identifier>(reverse, forward, primal, del, iterateInputs,
                                                             iterateOutputs);
      }
  };
}
