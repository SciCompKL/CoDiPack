#pragma once

#include <medi/adjointInterface.hpp>
#include <medi/adToolImplCommon.hpp>
#include <medi/adToolInterface.h>
#include <medi/ampi/ampiMisc.h>
#include <medi/ampi/op.hpp>
#include <medi/ampi/typeDefault.hpp>
#include <medi/ampi/types/indexTypeHelper.hpp>

#include "../../aux/macros.h"
#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../tapes/aux/adjointVectorAccess.hpp"
#include "../../tapes/interfaces/fullTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Type>
  struct CoDiMeDiAdjointInterfaceWrapper : public medi::AdjointInterface {
    public:

      using Type = DECLARE_DEFAULT(_Type, TEMPLATE(LhsExpressionInterface<double, double, ANY, ANY>));

      using Real = typename Type::Real;
      using Identifier = typename Type::Identifier;


      VectorAccessInterface<Real, Identifier>* codiInterface;

      int vecSize;

      CoDiMeDiAdjointInterfaceWrapper(VectorAccessInterface<Real, Identifier>* interface) :
        codiInterface(interface),
        vecSize((int)interface->getVectorSize()) {}

      int computeElements(int elements) const {
        return elements * vecSize;
      }

      int getVectorSize() const {
        return vecSize;
      }

      inline void getAdjoints(void const* i, void* a, int elements) const {
        Real* adjoints = (Real*)a;
        Identifier* indices = (Identifier*)i;

        for(int pos = 0; pos < elements; ++pos) {
          codiInterface->getAdjointVec(indices[pos], &adjoints[pos * vecSize]);
          codiInterface->resetAdjointVec(indices[pos]);
        }
      }

      inline void updateAdjoints(void const* i, void const* a, int elements) const {
        Real* adjoints = (Real*)a;
        Identifier* indices = (Identifier*)i;

        for(int pos = 0; pos < elements; ++pos) {

          codiInterface->updateAdjointVec(indices[pos], &adjoints[pos * vecSize]);
        }
      }

      inline void getPrimals(void const* i, void const* p, int elements) const {
        Real* primals = (Real*)p;
        Identifier* indices = (Identifier*)i;

        for(int pos = 0; pos < elements; ++pos) {
          primals[pos] = codiInterface->getPrimal(indices[pos]);
        }
      }

      inline void setPrimals(void const* i, void const* p, int elements) const {
        Real* primals = (Real*)p;
        Identifier* indices = (Identifier*)i;

        for(int pos = 0; pos < elements; ++pos) {
          codiInterface->setPrimal(indices[pos], primals[pos]);
        }
      }

      inline void combineAdjoints(void* b, int const elements, int const ranks) const {
        Real* buf = (Real*)b;

        for(int curRank = 1; curRank < ranks; ++curRank) {
          for(int curPos = 0; curPos < elements; ++curPos) {
            for(int dim = 0; dim < vecSize; ++dim) {

              buf[curPos * vecSize + dim] += buf[(elements * curRank + curPos) * vecSize + dim];
            }
          }
        }
      }

      inline void createPrimalTypeBuffer(void* &buf, size_t size) const {
        buf = (void*)(new Real[size * vecSize]);
      }

      inline void deletePrimalTypeBuffer(void* &b) const {
        if(NULL != b) {
          Real* buf = (Real*)b;
          delete [] buf;
          b = NULL;
        }
      }

      inline void createAdjointTypeBuffer(void* &buf, size_t size) const {
        buf = (void*)(new Real[size * vecSize]);
      }

      inline void deleteAdjointTypeBuffer(void* &b) const {
        if(NULL != b) {
          Real* buf = (Real*)b;
          delete [] buf;
          b = NULL;
        }
      }
  };

  template<typename _Type>
  struct CoDiPackReverseTool : public medi::ADToolImplCommon<CoDiPackReverseTool<_Type>, _Type::Tape::RequiresPrimalRestore, false, _Type, typename _Type::Gradient, typename _Type::Real, typename _Type::Identifier> {
    public:
      // All type definitions for the interface
      using Type = DECLARE_DEFAULT(_Type, TEMPLATE(LhsExpressionInterface<double, double, ANY, ANY>));
      using PrimalType = typename Type::Real;
      using AdjointType = void;
      using ModifiedType = Type;
      using IndexType = typename Type::Identifier;


      // Helper definition for CoDiPack
      using Tape = DECLARE_DEFAULT(typename Type::Tape, TEMPLATE(FullTapeInterface<double, double, int, ANY>));

      using OpHelper = medi::OperatorHelper<
                medi::FunctionHelper<
                    Type, Type, typename Type::PassiveReal, typename Type::Gradient, typename Type::Identifier, CoDiPackReverseTool
                >
              >;

      using Base = medi::ADToolImplCommon<CoDiPackReverseTool, Tape::RequiresPrimalRestore, false, Type, typename Type::Gradient, PrimalType, IndexType>;
    private:
      // Private structures for the implementation

      OpHelper opHelper;

    public:
      CoDiPackReverseTool(MPI_Datatype primalMpiType, MPI_Datatype adjointMpiType) :
        Base(primalMpiType, adjointMpiType),
        opHelper()
      {
        opHelper.init();
      }

      ~CoDiPackReverseTool() {
        opHelper.finalize();
      }

      // Implementation of the interface

      inline  bool isHandleRequired() const {
        // Handle creation is based on the CoDiPack tape activity. Only if the tape is recording the adjoint communication
        // needs to be evaluated.
        return getTape().isActive();
      }

      inline void startAssembly(medi::HandleBase* h) const {
        CODI_UNUSED(h);

        // No preparation required for CoDiPack
      }

      inline void addToolAction(medi::HandleBase* h) const {
        if(NULL != h) {
          getTape().pushExternalFunction(ExternalFunction::create(callHandleReverse, h, deleteHandle, callHandleForward, callHandlePrimal));
        }
      }

      medi::AMPI_Op convertOperator(medi::AMPI_Op op) const {
        return opHelper.convertOperator(op);
      }

      inline void stopAssembly(medi::HandleBase* h) const {
        CODI_UNUSED(h);

        // No preparation required for CoDiPack
      }

      static inline IndexType getIndex(Type const& value) {
        return value.getIdentifier();
      }

      static inline void registerValue(Type& value, PrimalType& oldPrimal, IndexType& index) {

        bool wasActive = getTape().isIdentifierActive(value.getIdentifier());
        value.getIdentifier() = IndexType();

        // make the value active again if it has been active before on the other processor
        if(wasActive) {
          if(Tape::LinearIndexHandling) {
            // value has been registered in createIndices
            value.getIdentifier() = index;

            // in createIndices the primal value has been set to zero. So set now the correct value
            if(Tape::HasPrimalValues) {
              getTape().setPrimal(index, value.getValue());
            }
            if(Tape::RequiresPrimalRestore) {
              oldPrimal = PrimalType(0.0);
            }
          } else {
            PrimalType primal = getTape().registerExternalFunctionOutput(value);
            if(Tape::RequiresPrimalRestore) {
              oldPrimal = primal;
            }
            index = value.getIdentifier();
          }
        } else {

          if(Tape::RequiresPrimalRestore) {
            oldPrimal = PrimalType(0.0);
          }
          if(!Tape::LinearIndexHandling) {
            index = getTape().getPassiveIndex();
          }
        }
      }

      static inline void clearIndex(Type& value) {
        IndexType oldIndex = value.getIdentifier();
        value.~Type();
        value.getIdentifier() = oldIndex;  // restore the index here so that the other side can decide of the communication was active or not
      }

      static inline void createIndex(Type& value, IndexType& index) {
        if(Tape::LinearIndexHandling) {
          getTape().registerInput(value);
          index = value.getIdentifier();
        }
      }

      static inline PrimalType getValue(Type const& value) {
        return value.getValue();
      }

      static inline void setIntoModifyBuffer(ModifiedType& modValue, Type const& value) {
        CODI_UNUSED(modValue, value);

        // CoDiPack values are send in place. No modified buffer is crated.
      }

      static inline void getFromModifyBuffer(ModifiedType const& modValue, Type& value) {
        CODI_UNUSED(modValue, value);

        // CoDiPack values are send in place. No modified buffer is crated.
      }

      static PrimalType getPrimalFromMod(ModifiedType const& modValue) {
        return modValue.value();
      }

      static void setPrimalToMod(ModifiedType& modValue, PrimalType const& value) {
        modValue.value() = value;
      }

      static void modifyDependency(ModifiedType& inval, ModifiedType& inoutval) {

        bool active = getTape().isIdentifierActive(inoutval.getIdentifier()) || getTape().isIdentifierActive(inval.getIdentifier());
        if(active) {
          inoutval.getIdentifier() = getTape().getInvalidIndex();
        } else {
          inoutval.getIdentifier() = getTape().getPassiveIndex();
        }
      }

    private:

      static void callHandleReverse(void* tape, void* h, void* ah) {
        medi::HandleBase* handle = static_cast<medi::HandleBase*>(h);
        CoDiMeDiAdjointInterfaceWrapper<Type> ahWrapper((VectorAccessInterface<PrimalType, IndexType>*)ah);
        handle->funcReverse(handle, &ahWrapper);
      }

      static void callHandleForward(void* tape, void* h, void* ah) {
        medi::HandleBase* handle = static_cast<medi::HandleBase*>(h);
        CoDiMeDiAdjointInterfaceWrapper<Type> ahWrapper((VectorAccessInterface<PrimalType, IndexType>*)ah);
        handle->funcForward(handle, &ahWrapper);
      }

      static void callHandlePrimal(void* tape, void* h, void* ah) {
        medi::HandleBase* handle = static_cast<medi::HandleBase*>(h);
        CoDiMeDiAdjointInterfaceWrapper<Type> ahWrapper((VectorAccessInterface<PrimalType, IndexType>*)ah);
        handle->funcPrimal(handle, &ahWrapper);
      }

      static void deleteHandle(void* tape, void* h) {
        CODI_UNUSED(tape);

        medi::HandleBase* handle = static_cast<medi::HandleBase*>(h);
        delete handle;
      }

      static Tape& getTape() {
        return Type::getGlobalTape();
      }
  };
}
