#pragma once

#include "../../aux/macros.h"
#include "../../aux/exceptions.hpp"
#include "../../config.h"
#include "../data/position.hpp"
#include "../../expressions/lhsExpressionInterface.hpp"


/** \copydoc codi::Namespace */
namespace codi {

  struct ExternalFunction {

    typedef void (*CallFunction)(void* tape, void* data, void* adjointInterface);
    typedef void (*DeleteFunction)(void* tape, void* data);

  private:
    CallFunction funcReverse;
    CallFunction funcForward;
    CallFunction funcPrimal;
    DeleteFunction funcDelete;

    void* data;

  public:
    ExternalFunction() {}
    ExternalFunction(CallFunction funcReverse, CallFunction funcForward, CallFunction funcPrimal, void* data, DeleteFunction funcDelete) :
      funcReverse(funcReverse),
      funcForward(funcForward),
      funcPrimal(funcPrimal),
      funcDelete(funcDelete),
      data(data){}

    static ExternalFunction create(CallFunction funcReverse,
                                   void* data,
                                   DeleteFunction funcDelete,
                                   CallFunction funcForward = nullptr,
                                   CallFunction funcPrimal = nullptr) {
      return ExternalFunction(funcReverse, funcForward, funcPrimal, data, funcDelete);
    }

    void deleteData(void* tape) {
      if (deleteCheckpoint != NULL){
        deleteCheckpoint(tape, data);
        data = NULL;
      }
    }

    void evaluateReverse(void* tape, void* adjointInterface) {
      if(NULL != funcReverse) {
        funcReverse(tape, data, adjointInterface);
      } else {
        CODI_EXCEPTION("Calling an external function in reverse mode without providing a reverse evaluation function.");
      }
    }

    void evaluateForward(void* tape, void* adjointInterface) {
      if(NULL != funcForward) {
        funcForward(tape, data, adjointInterface);
      } else {
        CODI_EXCEPTION("Calling an external function in forward mode without providing a forward evaluation function.");
      }
    }

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
