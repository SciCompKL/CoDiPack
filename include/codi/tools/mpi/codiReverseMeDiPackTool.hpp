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

#include <medi/adToolInterface.h>
#include <medi/ampi/ampiMisc.h>

#include <medi/adToolImplCommon.hpp>
#include <medi/adjointInterface.hpp>
#include <medi/ampi/op.hpp>
#include <medi/ampi/typeDefault.hpp>
#include <medi/ampi/types/indexTypeHelper.hpp>

#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../misc/macros.hpp"
#include "../../tapes/interfaces/fullTapeInterface.hpp"
#include "../../tapes/misc/adjointVectorAccess.hpp"

/** \copydoc codi::Namespace */
namespace codi {

#ifndef DOXYGEN_DISABLE

  template<typename T_Type>
  struct CoDiMeDiAdjointInterfaceWrapper : public medi::AdjointInterface {
    public:

      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);

      using Real = typename Type::Real;
      using Identifier = typename Type::Identifier;

      VectorAccessInterface<Real, Identifier>* codiInterface;

      int vecSize;

      CoDiMeDiAdjointInterfaceWrapper(VectorAccessInterface<Real, Identifier>* interface)
          : codiInterface(interface), vecSize((int)interface->getVectorSize()) {}

      CODI_INLINE int computeElements(int elements) const {
        return elements * vecSize;
      }

      CODI_INLINE int getVectorSize() const {
        return vecSize;
      }

      CODI_INLINE void getAdjoints(void const* i, void* a, int elements) const {
        Real* adjoints = (Real*)a;
        Identifier* indices = (Identifier*)i;

        for (int pos = 0; pos < elements; ++pos) {
          codiInterface->getAdjointVec(indices[pos], &adjoints[pos * vecSize]);
          codiInterface->resetAdjointVec(indices[pos]);
        }
      }

      CODI_INLINE void updateAdjoints(void const* i, void const* a, int elements) const {
        Real* adjoints = (Real*)a;
        Identifier* indices = (Identifier*)i;

        for (int pos = 0; pos < elements; ++pos) {
          codiInterface->updateAdjointVec(indices[pos], &adjoints[pos * vecSize]);
        }
      }

      CODI_INLINE void getPrimals(void const* i, void const* p, int elements) const {
        Real* primals = (Real*)p;
        Identifier* indices = (Identifier*)i;

        for (int pos = 0; pos < elements; ++pos) {
          primals[pos] = codiInterface->getPrimal(indices[pos]);
        }
      }

      CODI_INLINE void setPrimals(void const* i, void const* p, int elements) const {
        Real* primals = (Real*)p;
        Identifier* indices = (Identifier*)i;

        for (int pos = 0; pos < elements; ++pos) {
          codiInterface->setPrimal(indices[pos], primals[pos]);
        }
      }

      CODI_INLINE void combineAdjoints(void* b, int const elements, int const ranks) const {
        Real* buf = (Real*)b;

        for (int curRank = 1; curRank < ranks; ++curRank) {
          for (int curPos = 0; curPos < elements; ++curPos) {
            for (int dim = 0; dim < vecSize; ++dim) {
              buf[curPos * vecSize + dim] += buf[(elements * curRank + curPos) * vecSize + dim];
            }
          }
        }
      }

      CODI_INLINE void createPrimalTypeBuffer(void*& buf, size_t size) const {
        buf = (void*)(new Real[size * vecSize]);
      }

      CODI_INLINE void deletePrimalTypeBuffer(void*& b) const {
        if (nullptr != b) {
          Real* buf = (Real*)b;
          delete[] buf;
          b = nullptr;
        }
      }

      CODI_INLINE void createAdjointTypeBuffer(void*& buf, size_t size) const {
        buf = (void*)(new Real[size * vecSize]);
      }

      CODI_INLINE void deleteAdjointTypeBuffer(void*& b) const {
        if (nullptr != b) {
          Real* buf = (Real*)b;
          delete[] buf;
          b = nullptr;
        }
      }
  };

  template<typename T_Type>
  struct CoDiPackReverseTool
      : public medi::ADToolImplCommon<CoDiPackReverseTool<T_Type>, T_Type::Tape::RequiresPrimalRestore, false, T_Type,
                                      typename T_Type::Gradient, typename T_Type::Real, typename T_Type::Identifier> {
    public:

      // All type definitions for the interface.
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);
      using PrimalType = typename Type::Real;
      using AdjointType = void;
      using ModifiedType = Type;
      using IndexType = typename Type::Identifier;

      // Helper definition for CoDiPack.
      using Tape = CODI_DD(typename Type::Tape, CODI_DEFAULT_TAPE);

      using OpHelper =
          medi::OperatorHelper<medi::FunctionHelper<Type, Type, typename Type::PassiveReal, typename Type::Gradient,
                                                    typename Type::Identifier, CoDiPackReverseTool> >;

      using Base = medi::ADToolImplCommon<CoDiPackReverseTool, Tape::RequiresPrimalRestore, false, Type,
                                          typename Type::Gradient, PrimalType, IndexType>;

    private:
      // Private structures for the implementation.

      OpHelper opHelper;

    public:
      CoDiPackReverseTool(MPI_Datatype primalMpiType, MPI_Datatype adjointMpiType)
          : Base(primalMpiType, adjointMpiType), opHelper() {
        opHelper.init();
      }

      ~CoDiPackReverseTool() {
        opHelper.finalize();
      }

      // Implementation of the interface.

      CODI_INLINE bool isHandleRequired() const {
        // Handle creation is based on the CoDiPack tape activity. Only if the tape is recording the adjoint
        // communication needs to be evaluated.
        return getTape().isActive();
      }

      CODI_INLINE void startAssembly(medi::HandleBase* h) const {
        CODI_UNUSED(h);

        // No preparation required for CoDiPack.
      }

      CODI_INLINE void addToolAction(medi::HandleBase* h) const {
        if (nullptr != h) {
          getTape().pushExternalFunction(
              ExternalFunction<Tape>::create(callHandleReverse, h, deleteHandle, callHandleForward, callHandlePrimal));
        }
      }

      medi::AMPI_Op convertOperator(medi::AMPI_Op op) const {
        return opHelper.convertOperator(op);
      }

      CODI_INLINE void stopAssembly(medi::HandleBase* h) const {
        CODI_UNUSED(h);

        // No preparation required for CoDiPack.
      }

      static CODI_INLINE IndexType getIndex(Type const& value) {
        return value.getIdentifier();
      }

      static CODI_INLINE void registerValue(Type& value, PrimalType& oldPrimal, IndexType& index) {
        bool wasActive = getTape().isIdentifierActive(value.getIdentifier());
        value.getIdentifier() = IndexType();

        // Make the value active again if it has been active before on the other processor.
        if (wasActive) {
          if (Tape::LinearIndexHandling) {
            // Value has been registered in createIndices.
            value.getIdentifier() = index;

            // In createIndices the primal value has been set to zero. So set now the correct value.
            if (Tape::HasPrimalValues) {
              getTape().setPrimal(index, value.getValue());
            }
            if (Tape::RequiresPrimalRestore) {
              oldPrimal = PrimalType(0.0);
            }
          } else {
            PrimalType primal = getTape().registerExternalFunctionOutput(value);
            if (Tape::RequiresPrimalRestore) {
              oldPrimal = primal;
            }
            index = value.getIdentifier();
          }
        } else {
          if (Tape::RequiresPrimalRestore) {
            oldPrimal = PrimalType(0.0);
          }
          if (!Tape::LinearIndexHandling) {
            index = getTape().getPassiveIndex();
          }
        }
      }

      static CODI_INLINE void clearIndex(Type& value) {
        IndexType oldIndex = value.getIdentifier();
        value.~Type();
        value.getIdentifier() = oldIndex;  // Restore the index here so that the other side can decide of the
                                           // communication was active or not.
      }

      static CODI_INLINE void createIndex(Type& value, IndexType& index) {
        if (Tape::LinearIndexHandling) {
          IndexType oldIndex = value.getIdentifier();
          getTape().registerInput(value);
          index = value.getIdentifier();
          value.getIdentifier() = oldIndex;  // Restore the index here so that the other side can decide of the
                                             // communication was active or not.
        }
      }

      static CODI_INLINE PrimalType getValue(Type const& value) {
        return value.getValue();
      }

      static CODI_INLINE void setIntoModifyBuffer(ModifiedType& modValue, Type const& value) {
        CODI_UNUSED(modValue, value);

        // CoDiPack values are send in place. No modified buffer is crated.
      }

      static CODI_INLINE void getFromModifyBuffer(ModifiedType const& modValue, Type& value) {
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
        bool active = getTape().isIdentifierActive(inoutval.getIdentifier()) ||
                      getTape().isIdentifierActive(inval.getIdentifier());
        if (active) {
          inoutval.getIdentifier() = getTape().getInvalidIndex();
        } else {
          inoutval.getIdentifier() = getTape().getPassiveIndex();
        }
      }

    private:

      static void callHandleReverse(Tape* tape, void* h, VectorAccessInterface<PrimalType, IndexType>* ah) {
        CODI_UNUSED(tape);

        medi::HandleBase* handle = static_cast<medi::HandleBase*>(h);
        CoDiMeDiAdjointInterfaceWrapper<Type> ahWrapper(ah);
        handle->funcReverse(handle, &ahWrapper);
      }

      static void callHandleForward(Tape* tape, void* h, VectorAccessInterface<PrimalType, IndexType>* ah) {
        CODI_UNUSED(tape);

        medi::HandleBase* handle = static_cast<medi::HandleBase*>(h);
        CoDiMeDiAdjointInterfaceWrapper<Type> ahWrapper(ah);
        handle->funcForward(handle, &ahWrapper);
      }

      static void callHandlePrimal(Tape* tape, void* h, VectorAccessInterface<PrimalType, IndexType>* ah) {
        CODI_UNUSED(tape);

        medi::HandleBase* handle = static_cast<medi::HandleBase*>(h);
        CoDiMeDiAdjointInterfaceWrapper<Type> ahWrapper(ah);
        handle->funcPrimal(handle, &ahWrapper);
      }

      static void deleteHandle(Tape* tape, void* h) {
        CODI_UNUSED(tape);

        medi::HandleBase* handle = static_cast<medi::HandleBase*>(h);
        delete handle;
      }

      static Tape& getTape() {
        return Type::getTape();
      }
  };
#endif
}
