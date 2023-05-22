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

#include <medi/ampi/typeDefault.hpp>
#include <medi/ampi/types/indexTypeHelper.hpp>

#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

#ifndef DOXYGEN_DISABLE

  template<typename T_Type>
  struct CoDiPackForwardTool : public medi::ADToolBase<CoDiPackForwardTool<T_Type>, typename T_Type::Gradient,
                                                       typename T_Type::PassiveReal, int> {
    public:

      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);

      using PrimalType = typename Type::Real;
      using AdjointType = void;
      using ModifiedType = Type;
      using IndexType = int;

      using Base = medi::ADToolBase<CoDiPackForwardTool, typename Type::Gradient, typename Type::PassiveReal, int>;

      using OpHelper =
          medi::OperatorHelper<medi::FunctionHelper<Type, Type, typename Type::PassiveReal, typename Type::Identifier,
                                                    typename Type::Gradient, CoDiPackForwardTool<Type> > >;

    private:
      OpHelper opHelper;

    public:

      CoDiPackForwardTool(MPI_Datatype primalMpiType, MPI_Datatype adjointMpiType)
          : Base(primalMpiType, adjointMpiType), opHelper() {
        opHelper.init();
      }

      ~CoDiPackForwardTool() {
        opHelper.finalize();
      }

      CODI_INLINE bool isActiveType() const {
        return false;
      }

      CODI_INLINE bool isHandleRequired() const {
        return false;
      }

      CODI_INLINE bool isModifiedBufferRequired() const {
        return false;
      }

      CODI_INLINE bool isOldPrimalsRequired() const {
        return false;
      }

      CODI_INLINE void startAssembly(medi::HandleBase* h) const {
        CODI_UNUSED(h);
      }

      CODI_INLINE void addToolAction(medi::HandleBase* h) const {
        CODI_UNUSED(h);
      }

      CODI_INLINE void stopAssembly(medi::HandleBase* h) const {
        CODI_UNUSED(h);
      }

      medi::AMPI_Op convertOperator(medi::AMPI_Op op) const {
        return opHelper.convertOperator(op);
      }

      CODI_INLINE void createPrimalTypeBuffer(PrimalType*& buf, size_t size) const {
        buf = new PrimalType[size];
      }

      CODI_INLINE void createIndexTypeBuffer(IndexType*& buf, size_t size) const {
        buf = new IndexType[size];
      }

      CODI_INLINE void deletePrimalTypeBuffer(PrimalType*& buf) const {
        if (nullptr != buf) {
          delete[] buf;
          buf = nullptr;
        }
      }

      CODI_INLINE void deleteIndexTypeBuffer(IndexType*& buf) const {
        if (nullptr != buf) {
          delete[] buf;
          buf = nullptr;
        }
      }

      static CODI_INLINE int getIndex(Type const& value) {
        return value.getIdentifier();
      }

      static CODI_INLINE void clearIndex(Type& value) {
        value.~Type();
        value.getIdentifier() = 0;
      }

      static CODI_INLINE void createIndex(Type& value, int& index) {
        CODI_UNUSED(value);
        index = 0;
      }

      static CODI_INLINE PrimalType getValue(Type const& value) {
        return value.getValue();
      }

      static CODI_INLINE void setIntoModifyBuffer(ModifiedType& modValue, Type const& value) {
        CODI_UNUSED(modValue, value);
      }

      static CODI_INLINE void getFromModifyBuffer(ModifiedType const& modValue, Type& value) {
        CODI_UNUSED(modValue, value);
      }

      static CODI_INLINE void registerValue(Type& value, PrimalType& oldValue, int& index) {
        CODI_UNUSED(value, oldValue, index);
      }

      static PrimalType getPrimalFromMod(ModifiedType const& modValue) {
        return modValue.value();
      }

      static void setPrimalToMod(ModifiedType& modValue, PrimalType const& value) {
        modValue.value() = value;
      }

      static void modifyDependency(ModifiedType& inval, ModifiedType& inoutval) {
        CODI_UNUSED(inval, inoutval);
      }
  };
#endif
}
