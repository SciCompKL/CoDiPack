#pragma once

#include "../../aux/macros.hpp"
#include "../../aux/exceptions.hpp"
#include "../../config.h"
#include "../data/position.hpp"
#include "../../expressions/lhsExpressionInterface.hpp"


/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief User defined evaluation functions for the taping process.
   *
   * See ExternalFunctionTapeInterface for details.
   *
   * The user can provide call functions for the reverse, forward and primal evaluation of a tape. These need to be of
   * the type CallFunction. The user has to cast all arguments to the correct type:
   *  - tape: The type of the tape on which this object was registered with `registerExternalFunction`
   *          (ExternalFunctionTapeInterface)
   *  - data: User provided data, type is known by the user.
   *  - adjointInterface: VectorAccessInterface with <typename Tape::Real, typename Tape::Identifier>
   *
   * The tape pointer can be used for general access to the tape. For all access to the gradient data the
   * adjointInterface should be used. If in the application no custom adjoint vectors are used, then the tape pointer
   * can also be used for the gradient data access.
   *
   * The delete function is called when the entry of the tape for the external function is deleted.
   */
  struct ExternalFunction {
    public:

      typedef void (*CallFunction)(void* tape, void* data, void* adjointInterface); ///< Call function definition.
      typedef void (*DeleteFunction)(void* tape, void* data); ///< Delete function definition.

    private:
      CallFunction funcReverse;
      CallFunction funcForward;
      CallFunction funcPrimal;
      DeleteFunction funcDelete;

      void* data;

    public:
      /// Constructor
      ExternalFunction() {}

      /// Any arguments can be NULL if not required.
      ExternalFunction(CallFunction funcReverse, CallFunction funcForward, CallFunction funcPrimal, void* data, DeleteFunction funcDelete) :
        funcReverse(funcReverse),
        funcForward(funcForward),
        funcPrimal(funcPrimal),
        funcDelete(funcDelete),
        data(data){}

      /// Helper function for the creation of an ExternalFunction object.
      static ExternalFunction create(CallFunction funcReverse,
                                     void* data,
                                     DeleteFunction funcDelete,
                                     CallFunction funcForward = nullptr,
                                     CallFunction funcPrimal = nullptr) {
        return ExternalFunction(funcReverse, funcForward, funcPrimal, data, funcDelete);
      }

      /// Calls the delete function if not NULL.
      void deleteData(void* tape) {
        if (funcDelete != NULL){
          funcDelete(tape, data);
          data = NULL;
        }
      }

      /// Calls the reverse function if not NULL, otherwise throws a CODI_EXCEPTION.
      void evaluateReverse(void* tape, void* adjointInterface) {
        if(NULL != funcReverse) {
          funcReverse(tape, data, adjointInterface);
        } else {
          CODI_EXCEPTION("Calling an external function in reverse mode without providing a reverse evaluation function.");
        }
      }

      /// Calls the forward function if not NULL, otherwise throws a CODI_EXCEPTION.
      void evaluateForward(void* tape, void* adjointInterface) {
        if(NULL != funcForward) {
          funcForward(tape, data, adjointInterface);
        } else {
          CODI_EXCEPTION("Calling an external function in forward mode without providing a forward evaluation function.");
        }
      }

      /// Calls the primal function if not NULL, otherwise throws a CODI_EXCEPTION.
      void evaluatePrimal(void* tape, void* adjointInterface) {
        if(NULL != funcPrimal) {
          funcPrimal(tape, data, adjointInterface);
        } else {
          CODI_EXCEPTION("Calling an external function in primal mode without providing a primal evaluation function.");
        }
      }
  };

  // TODO: Add typed interface
}
