#pragma once

#include "../../aux/exceptions.hpp"
#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../aux/vectorAccessInterface.hpp"
#include "../data/position.hpp"
#include "../interfaces/externalFunctionTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Internal untyped data for an external function.
  struct ExternalFunctionInternalData {
    protected:
      typedef void (*CallFunctionUntyped)(void* tape, void* data,
                                          void* adjointInterface);    ///< Call function definition.
      typedef void (*DeleteFunctionUntyped)(void* tape, void* data);  ///< Delete function definition.

      CallFunctionUntyped funcReverse;   ///< Reverse evaluation function pointer.
      CallFunctionUntyped funcForward;   ///< Forward evaluation function pointer.
      CallFunctionUntyped funcPrimal;    ///< Primal evaluation function pointer.
      DeleteFunctionUntyped funcDelete;  ///< User data deletion function pointer.

      void* data;  ///< User data pointer.

    public:

      /// Constructor
      ExternalFunctionInternalData()
          : funcReverse(NULL), funcForward(NULL), funcPrimal(NULL), funcDelete(NULL), data(NULL) {}

    protected:

      /// Constructor
      ExternalFunctionInternalData(CallFunctionUntyped funcReverse, CallFunctionUntyped funcForward,
                                   CallFunctionUntyped funcPrimal, DeleteFunctionUntyped funcDelete, void* data)
          : funcReverse(funcReverse),
            funcForward(funcForward),
            funcPrimal(funcPrimal),
            funcDelete(funcDelete),
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
   */
  template<typename _Tape>
  struct ExternalFunction : public ExternalFunctionInternalData {
    public:

      using Tape = CODI_DD(_Tape,
                           CODI_T(ExternalFunctionTapeInterface<double, double, int>));  ///< See ExternalFunction

      using VectorAccess =
          VectorAccessInterface<typename Tape::Real, typename Tape::Identifier>;  ///< Shortcut for
                                                                                  ///< VectorAccessInterface.
      typedef void (*CallFunction)(Tape* tape, void* data,
                                   VectorAccess* adjointInterface);  ///< Call function definition.
      typedef void (*DeleteFunction)(Tape* tape, void* data);        ///< Delete function definition.

      /// Any arguments can be NULL if not required.
      ExternalFunction(CallFunction funcReverse, CallFunction funcForward, CallFunction funcPrimal, void* data,
                       DeleteFunction funcDelete)
          : ExternalFunctionInternalData((ExternalFunctionInternalData::CallFunctionUntyped)funcReverse,
                                         (ExternalFunctionInternalData::CallFunctionUntyped)funcForward,
                                         (ExternalFunctionInternalData::CallFunctionUntyped)funcPrimal,
                                         (ExternalFunctionInternalData::DeleteFunctionUntyped)funcDelete, data) {}

      /// Helper function for the creation of an ExternalFunction object.
      static ExternalFunction create(CallFunction funcReverse, void* data, DeleteFunction funcDelete,
                                     CallFunction funcForward = nullptr, CallFunction funcPrimal = nullptr) {
        return ExternalFunction(funcReverse, funcForward, funcPrimal, data, funcDelete);
      }

      /// Calls the delete function if not NULL.
      void deleteData(Tape* tape) {
        if (funcDelete != NULL) {
          funcDelete(tape, data);
          data = NULL;
        }
      }

      /// Calls the reverse function if not NULL, otherwise throws a CODI_EXCEPTION.
      void evaluateReverse(Tape* tape, VectorAccess* adjointInterface) {
        if (NULL != funcReverse) {
          funcReverse(tape, data, adjointInterface);
        } else {
          CODI_EXCEPTION(
              "Calling an external function in reverse mode without providing a reverse evaluation function.");
        }
      }

      /// Calls the forward function if not NULL, otherwise throws a CODI_EXCEPTION.
      void evaluateForward(Tape* tape, VectorAccess* adjointInterface) {
        if (NULL != funcForward) {
          funcForward(tape, data, adjointInterface);
        } else {
          CODI_EXCEPTION(
              "Calling an external function in forward mode without providing a forward evaluation function.");
        }
      }

      /// Calls the primal function if not NULL, otherwise throws a CODI_EXCEPTION.
      void evaluatePrimal(Tape* tape, VectorAccess* adjointInterface) {
        if (NULL != funcPrimal) {
          funcPrimal(tape, data, adjointInterface);
        } else {
          CODI_EXCEPTION("Calling an external function in primal mode without providing a primal evaluation function.");
        }
      }
  };
}
