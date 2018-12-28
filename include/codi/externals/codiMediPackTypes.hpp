/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#pragma once

#include <medi/adjointInterface.hpp>
#include <medi/adToolImplCommon.hpp>
#include <medi/adToolInterface.h>
#include <medi/ampi/op.hpp>
#include <medi/ampi/ampiMisc.h>
#include <medi/ampi/typeDefault.hpp>
#include <medi/ampi/types/indexTypeHelper.hpp>

#include "../adjointInterface.hpp"

template<typename CoDiType>
struct CoDiMeDiAdjointInterfaceWrapper : public medi::AdjointInterface {

    typedef typename CoDiType::Real Real;
    typedef typename CoDiType::GradientData IndexType;

    codi::AdjointInterface<Real, IndexType>* codiInterface;

    int vecSize;

    CoDiMeDiAdjointInterfaceWrapper(codi::AdjointInterface<Real, IndexType>* interface) :
      codiInterface(interface),
      vecSize((int)interface->getVectorSize()) {}

    int computeElements(int elements) const {
      return elements * vecSize;
    }

    int getVectorSize() const {
      return vecSize;
    }

    inline void getAdjoints(const void* i, void* a, int elements) const {
      Real* adjoints = (Real*)a;
      IndexType* indices = (IndexType*)i;

      for(int pos = 0; pos < elements; ++pos) {
        codiInterface->getAdjointVec(indices[pos], &adjoints[pos * vecSize]);
        codiInterface->resetAdjointVec(indices[pos]);
      }
    }

    inline void updateAdjoints(const void* i, const void* a, int elements) const {
      Real* adjoints = (Real*)a;
      IndexType* indices = (IndexType*)i;

      for(int pos = 0; pos < elements; ++pos) {

        codiInterface->updateAdjointVec(indices[pos], &adjoints[pos * vecSize]);
      }
    }

    inline void getPrimals(const void* i, const void* p, int elements) const {
      Real* primals = (Real*)p;
      IndexType* indices = (IndexType*)i;

      for(int pos = 0; pos < elements; ++pos) {
        primals[pos] = codiInterface->getPrimal(indices[pos]);
      }
    }

    inline void setPrimals(const void* i, const void* p, int elements) const {
      Real* primals = (Real*)p;
      IndexType* indices = (IndexType*)i;

      for(int pos = 0; pos < elements; ++pos) {
        codiInterface->setPrimal(indices[pos], primals[pos]);
      }
    }

    inline void combineAdjoints(void* b, const int elements, const int ranks) const {
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

template<typename CoDiType>
struct CoDiPackTool : public medi::ADToolImplCommon<CoDiPackTool<CoDiType>, CoDiType::TapeType::RequiresPrimalReset, false, CoDiType, typename CoDiType::GradientValue, typename CoDiType::Real, typename CoDiType::GradientData> {
  public:
    // All type definitions for the interface
    typedef CoDiType Type;
    typedef typename CoDiType::Real PrimalType;
    typedef void AdjointType;
    typedef CoDiType ModifiedType;
    typedef typename CoDiType::GradientData IndexType;

    // Helper definition for CoDiPack
    typedef typename CoDiType::TapeType Tape;
    typedef medi::MpiTypeDefault<CoDiPackTool<CoDiType>> MediType;

    // Static structures for the interface
    static MediType* MPI_TYPE;
    static medi::AMPI_Datatype MPI_INT_TYPE;

    static MPI_Datatype MpiType;
    static MPI_Datatype ModifiedMpiType;
    static MPI_Datatype PrimalMpiType;
    static MPI_Datatype AdjointMpiType;
  private:
    // Private structures for the implemenation
    static medi::OperatorHelper<
              medi::FunctionHelper<
                  CoDiType, CoDiType, typename CoDiType::PassiveReal, typename CoDiType::GradientData, typename CoDiType::GradientValue, CoDiPackTool<CoDiType>
              >
            > operatorHelper;

    static Tape* adjointTape;

  public:
    CoDiPackTool(MPI_Datatype primalMpiType, MPI_Datatype adjointMpiType) :
      medi::ADToolImplCommon<CoDiPackTool<CoDiType>, CoDiType::TapeType::RequiresPrimalReset, false, CoDiType, typename CoDiType::GradientValue, typename CoDiType::Real, typename CoDiType::GradientData>(primalMpiType, adjointMpiType) {}

    // Implementation of the interface

    static void init() {
      initTypes();

      MPI_TYPE = new MediType();

      operatorHelper.init(MPI_TYPE);
      MPI_INT_TYPE = operatorHelper.MPI_INT_TYPE;
    }

    static void finalize() {

      operatorHelper.finalize();

      if(nullptr != MPI_TYPE) {
        delete MPI_TYPE;
        MPI_TYPE = nullptr;
      }

      finalizeTypes();
    }

    inline  bool isHandleRequired() const {
      // Handle creation is based on the CoDiPack tape activity. Only if the tape is recording the adjoint communication
      // needs to be evaluated.
      return Type::getGlobalTape().isActive();
    }

    inline void startAssembly(medi::HandleBase* h) const {
      MEDI_UNUSED(h);

      // No preperation required for CoDiPack
    }

    inline void addToolAction(medi::HandleBase* h) const {
      if(NULL != h) {
        Type::getGlobalTape().pushExternalFunctionHandle(callHandleReverse, h, deleteHandle, callHandleForward, callHandlePrimal);
      }
    }

    medi::AMPI_Op convertOperator(medi::AMPI_Op op) const {
      return operatorHelper.convertOperator(op);
    }

    inline void stopAssembly(medi::HandleBase* h) const {
      MEDI_UNUSED(h);

      // No preperation required for CoDiPack
    }

    static inline IndexType getIndex(const Type& value) {
      return value.getGradientData();
    }

    static inline void registerValue(Type& value, PrimalType& oldPrimal, IndexType& index) {

      bool wasActive = 0 != value.getGradientData();
      value.getGradientData() = IndexType();

      // make the value active again if it has been active before on the other processor
      if(wasActive) {
        if(CoDiType::TapeType::LinearIndexHandler) {
          // value has been registered in createIndices
          value.getGradientData() = index;

          // in createIndices the primal value has been set to zero. So set now the correct value
          Type::getGlobalTape().setPrimalValue(index, value.getValue());
          if(CoDiType::TapeType::RequiresPrimalReset) {
            oldPrimal = PrimalType(0.0);
          }
        } else {
          PrimalType primal = Type::getGlobalTape().registerExtFunctionOutput(value);
          if(CoDiType::TapeType::RequiresPrimalReset) {
            oldPrimal = primal;
          }
          index = value.getGradientData();
        }
      } else {

        if(CoDiType::TapeType::RequiresPrimalReset) {
          oldPrimal = PrimalType(0.0);
        }
        if(!CoDiType::TapeType::LinearIndexHandler) {
          index = Type::getGlobalTape().getPassiveIndex();
        }
      }
    }

    static inline void clearIndex(Type& value) {
      IndexType oldIndex = value.getGradientData();
      value.~Type();
      value.getGradientData() = oldIndex;  // restore the index here so that the other side can decide of the communication was active or not
    }

    static inline void createIndex(Type& value, IndexType& index) {
      if(CoDiType::TapeType::LinearIndexHandler) {
        Type::getGlobalTape().registerInput(value);
        index = value.getGradientData();
      }
    }

    static inline PrimalType getValue(const Type& value) {
      return value.getValue();
    }

    static inline void setIntoModifyBuffer(ModifiedType& modValue, const Type& value) {
      MEDI_UNUSED(modValue);
      MEDI_UNUSED(value);

      // CoDiPack values are send in place. No modified buffer is crated.
    }

    static inline void getFromModifyBuffer(const ModifiedType& modValue, Type& value) {
      MEDI_UNUSED(modValue);
      MEDI_UNUSED(value);

      // CoDiPack values are send in place. No modified buffer is crated.
    }

    static PrimalType getPrimalFromMod(const ModifiedType& modValue) {
      return modValue.value();
    }

    static void setPrimalToMod(ModifiedType& modValue, const PrimalType& value) {
      modValue.value() = value;
    }

    static void modifyDependency(ModifiedType& inval, ModifiedType& inoutval) {

      bool active = (0 != inoutval.getGradientData()) || (0 != inval.getGradientData());
      if(active) {
        inoutval.getGradientData() = Type::getGlobalTape().getInvalidIndex();
      } else {
        inoutval.getGradientData() = Type::getGlobalTape().getPassiveIndex();
      }
    }

  private:
    // Helper functions for the implementation
    static void finalizeTypes() {
      MPI_Type_free(&MpiType);
    }

    static void initTypes() {
      // create the mpi type for CoDiPack
      // this type is used in this type and the passive formulation
      // TODO: add proper type creation
      MPI_Type_contiguous(sizeof(CoDiType), MPI_BYTE, &MpiType);
      MPI_Type_commit(&MpiType);

      ModifiedMpiType = MpiType;

      MPI_Type_contiguous(sizeof(typename CoDiType::Real), MPI_BYTE, &PrimalMpiType);
      MPI_Type_commit(&PrimalMpiType);

      // Since we use the CoDiPack adjoint interface, everything is interpreted in terms of the primal computation type
      // TODO: add proper type creation
      AdjointMpiType = PrimalMpiType;
    }

    static void callHandleReverse(void* tape, void* h, void* ah) {
      adjointTape = (Tape*)tape;
      medi::HandleBase* handle = static_cast<medi::HandleBase*>(h);
      CoDiMeDiAdjointInterfaceWrapper<CoDiType> ahWrapper((codi::AdjointInterface<typename CoDiType::Real, typename CoDiType::GradientData>*)ah);
      handle->funcReverse(handle, &ahWrapper);
    }

    static void callHandleForward(void* tape, void* h, void* ah) {
      adjointTape = (Tape*)tape;
      medi::HandleBase* handle = static_cast<medi::HandleBase*>(h);
      CoDiMeDiAdjointInterfaceWrapper<CoDiType> ahWrapper((codi::AdjointInterface<typename CoDiType::Real, typename CoDiType::GradientData>*)ah);
      handle->funcForward(handle, &ahWrapper);
    }

    static void callHandlePrimal(void* tape, void* h, void* ah) {
      adjointTape = (Tape*)tape;
      medi::HandleBase* handle = static_cast<medi::HandleBase*>(h);
      CoDiMeDiAdjointInterfaceWrapper<CoDiType> ahWrapper((codi::AdjointInterface<typename CoDiType::Real, typename CoDiType::GradientData>*)ah);
      handle->funcPrimal(handle, &ahWrapper);
    }

    static void deleteHandle(void* tape, void* h) {
      MEDI_UNUSED(tape);

      medi::HandleBase* handle = static_cast<medi::HandleBase*>(h);
      delete handle;
    }


};

template<typename CoDiType> MPI_Datatype CoDiPackTool<CoDiType>::MpiType;
template<typename CoDiType> MPI_Datatype CoDiPackTool<CoDiType>::ModifiedMpiType;
template<typename CoDiType> MPI_Datatype CoDiPackTool<CoDiType>::PrimalMpiType;
template<typename CoDiType> MPI_Datatype CoDiPackTool<CoDiType>::AdjointMpiType;
template<typename CoDiType> typename CoDiPackTool<CoDiType>::MediType* CoDiPackTool<CoDiType>::MPI_TYPE;
template<typename CoDiType> medi::AMPI_Datatype CoDiPackTool<CoDiType>::MPI_INT_TYPE;
template<typename CoDiType> medi::OperatorHelper<medi::FunctionHelper<CoDiType, CoDiType, typename CoDiType::PassiveReal, typename CoDiType::GradientData, typename CoDiType::GradientValue, CoDiPackTool<CoDiType> > > CoDiPackTool<CoDiType>::operatorHelper;
template<typename CoDiType> typename CoDiType::TapeType* CoDiPackTool<CoDiType>::adjointTape;
