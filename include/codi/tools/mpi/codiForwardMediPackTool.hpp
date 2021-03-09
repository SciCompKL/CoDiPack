#pragma once

#include <medi/adToolInterface.h>
#include <medi/ampi/ampiMisc.h>

#include <medi/ampi/typeDefault.hpp>
#include <medi/ampi/types/indexTypeHelper.hpp>

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

#ifndef DOXYGEN_DISABLE

  template<typename _Type>
  struct CoDiPackForwardTool : public medi::ADToolBase<CoDiPackForwardTool<_Type>, typename _Type::Gradient,
                                                       typename _Type::PassiveReal, int> {
    public:

      using Type = CODI_DD(_Type, CODI_T(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

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
        if (NULL != buf) {
          delete[] buf;
          buf = NULL;
        }
      }

      CODI_INLINE void deleteIndexTypeBuffer(IndexType*& buf) const {
        if (NULL != buf) {
          delete[] buf;
          buf = NULL;
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
