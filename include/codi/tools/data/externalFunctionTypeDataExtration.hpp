#pragma once

#include <complex>
#include <vector>

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../tapes/aux/vectorAccessInterface.hpp"
#include "../../traits/computationTraits.hpp"
#include "../../traits/expressionTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Type, typename = void>
  struct VectorAccessTypeWrapper : public VectorAccessInterface<CODI_ANY, CODI_ANY> {
      static_assert(false && std::is_void<_Type>::value, "Instantiation of unspecialized VectorAccessTypeWrapper.");

      using Type = CODI_DD(_Type, CODI_ANY);

      using Real = CODI_ANY;
      using Identifier = CODI_ANY;

      /*******************************************************************************/
      /// @name Static data extraction and registration functions

      static Real getValue(Type const& v);

      static Identifier getIdentifier(Type const& v);

      static void setValue(Type& v, Real const& value);

      static Real registerExternalFunctionOutput(Type& v);
  };

  template<typename _Real, typename _Identifier, typename _InnerInterface>
  struct VectorAccessTypeWrapperBase : public VectorAccessInterface<_Real, _Identifier> {
    public:

      using Real = CODI_DD(_Real, CODI_ANY);
      using Identifier = CODI_DD(_Identifier, CODI_ANY);

      using InnerInterface = CODI_DD(_InnerInterface,
                                     CODI_T(VectorAccessInterface<double, int>));

    protected:

      InnerInterface& innerInterface;

      std::vector<Real> lhs;

    public:

      VectorAccessTypeWrapperBase(InnerInterface* innerInterface)
          : innerInterface(*innerInterface), lhs(innerInterface->getVectorSize()) {}

      /*******************************************************************************/
      /// @name Misc

      size_t getVectorSize() const {
        return innerInterface.getVectorSize();
      }

      bool isLhsZero() {
        bool isZero = true;

        for (size_t curDim = 0; isZero && curDim < lhs.size(); curDim += 1) {
          isZero &= RealTraits::isTotalZero(lhs[curDim]);
        }

        return isZero;
      }

      /*******************************************************************************/
      /// @name Indirect adjoint access

      void setLhsAdjoint(Identifier const& index) {
        getAdjointVec(index, lhs.data());
        this->resetAdjointVec(index);
      }

      void updateAdjointWithLhs(Identifier const& index, Real const& jacobi) {
        for (size_t curDim = 0; curDim < lhs.size(); curDim += 1) {
          Real update = jacobi * lhs[curDim];
          this->updateAdjoint(index, curDim, update);
        }
      }

      /*******************************************************************************/
      /// @name Indirect tangent access

      void setLhsTangent(Identifier const& index) {
        updateAdjointVec(index, lhs.data());

        for (size_t curDim = 0; curDim < lhs.size(); curDim += 1) {
          lhs[curDim] = Real();
        }
      }

      void updateTangentWithLhs(Identifier const& index, Real const& jacobi) {
        for (size_t curDim = 0; curDim < lhs.size(); curDim += 1) {
          lhs[curDim] += jacobi * this->getAdjoint(index, curDim);
        }
      }

      /*******************************************************************************/
      /// @name Direct adjoint access

      void getAdjointVec(Identifier const& index, Real* const vec) {
        for (size_t curDim = 0; curDim < lhs.size(); curDim += 1) {
          vec[curDim] = this->getAdjoint(index, curDim);
        }
      }

      void updateAdjointVec(Identifier const& index, Real const* const vec) {
        for (size_t curDim = 0; curDim < lhs.size(); curDim += 1) {
          this->updateAdjoint(index, curDim, vec[curDim]);
        }
      }

      /*******************************************************************************/
      /// @name Primal access

      bool hasPrimals() {
        return innerInterface.hasPrimals();
      }
  };

#ifndef DOXYGEN_DISABLE
  template<typename _InnerType>
  struct VectorAccessTypeWrapper<std::complex<_InnerType>>
      : public VectorAccessTypeWrapperBase<
            std::complex<typename _InnerType::Real>, std::complex<typename _InnerType::Identifier>,
            VectorAccessInterface<typename _InnerType::Real, typename _InnerType::Identifier>> {
    public:

      using InnerType = CODI_DD(
          _InnerType,
          CODI_T(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
      using Type = std::complex<InnerType>;

      using InnerInterface =
          VectorAccessInterface<typename InnerType::Real,
                                typename InnerType::Identifier>;

      using Real = std::complex<typename InnerType::Real>;
      using Identifier =
          std::complex<typename InnerType::Identifier>;

      using Base = VectorAccessTypeWrapperBase<Real, Identifier, InnerInterface>;

      VectorAccessTypeWrapper(InnerInterface* innerInterface) : Base(innerInterface) {}

      /*******************************************************************************/
      /// @name Static data extraction and registration functions

      CODI_INLINE static Real getValue(Type const& v) {
        return Real(std::real(v).getValue(), std::imag(v).getValue());
      }

      CODI_INLINE static Identifier getIdentifier(Type const& v) {
        InnerType const* vArray = reinterpret_cast<InnerType const*>(&v);
        return Real(vArray[0].getIdentifier(), vArray[1].getIdentifier());
      }

      CODI_INLINE static void setValue(Type& v, Real const& value) {
        InnerType* vArray = reinterpret_cast<InnerType*>(&v);

        vArray[0].setValue(std::real(value));
        vArray[1].setValue(std::imag(value));
      }

      CODI_INLINE static Real registerExternalFunctionOutput(Type& v) {
        InnerType* vArray = reinterpret_cast<InnerType*>(&v);

        return Real(InnerType::getGlobalTape().registerExternalFunctionOutput(vArray[0]),
                    InnerType::getGlobalTape().registerExternalFunctionOutput(vArray[1]));
      }

      /*******************************************************************************/
      /// @name Direct adjoint access

      void resetAdjoint(Identifier const& index, size_t dim) {
        Base::innerInterface.resetAdjoint(std::real(index), dim);
        Base::innerInterface.resetAdjoint(std::imag(index), dim);
      }

      void resetAdjointVec(Identifier const& index) {
        Base::innerInterface.resetAdjointVec(std::real(index));
        Base::innerInterface.resetAdjointVec(std::imag(index));
      }

      Real getAdjoint(Identifier const& index, size_t dim) {
        return Real(Base::innerInterface.getAdjoint(std::real(index), dim),
                    Base::innerInterface.getAdjoint(std::imag(index), dim));
      }

      void updateAdjoint(Identifier const& index, size_t dim, Real const& adjoint) {
        Base::innerInterface.updateAdjoint(std::real(index), dim, std::real(adjoint));
        Base::innerInterface.updateAdjoint(std::imag(index), dim, std::imag(adjoint));
      }

      /*******************************************************************************/
      /// @name Primal access

      void setPrimal(Identifier const& index, Real const& primal) {
        Base::innerInterface.setPrimal(std::real(index), std::real(primal));
        Base::innerInterface.setPrimal(std::imag(index), std::imag(primal));
      }

      Real getPrimal(Identifier const& index) {
        return Real(Base::innerInterface.getPrimal(std::real(index)), Base::innerInterface.getPrimal(std::imag(index)));
      }
  };
#endif

  template<typename _Type, typename = void>
  struct ExternalFunctionTypeDataExtraction {
      using Type = CODI_DD(_Type, CODI_ANY);

      using VectorWrapper =
          VectorAccessTypeWrapper<Type>;

      using Real = typename VectorWrapper::Real;
      using Identifier =
          typename VectorWrapper::Identifier;

      CODI_INLINE static Real getValue(Type const& v) {
        return VectorWrapper::getValue(v);
      }

      CODI_INLINE static Identifier getIdentifier(Type const& v) {
        return VectorWrapper::getIdentifier(v);
      }

      CODI_INLINE static void setValue(Type& v, Real const& value) {
        VectorWrapper::setValue(v, value);
      }

      CODI_INLINE static Real registerExternalFunctionOutput(Type& v) {
        return VectorWrapper::registerExternalFunctionOutput(v);
      }

      template<typename RealInterface, typename IdentifierInterface>
      static VectorWrapper* createVectorInterfaceWrapper(
          VectorAccessInterface<RealInterface, IdentifierInterface>* accessInterface) {
        return new VectorWrapper(accessInterface);
      }

      static void destroyVectorInterfaceWrapper(VectorWrapper* wrapper) {
        delete wrapper;
      }
  };

#ifndef DOXYGEN_DISABLE
  template<typename _ActiveType>
  struct ExternalFunctionTypeDataExtraction<_ActiveType, ExpressionTraits::EnableIfLhsExpression<_ActiveType>> {
      using Type = CODI_DD(
          _ActiveType,
          CODI_T(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

      using VectorWrapper =
          VectorAccessInterface<typename Type::Real,
                                typename Type::Identifier>;

      using Real = typename Type::Real;
      using Identifier = typename Type::Identifier;

      CODI_INLINE static Real getValue(Type const& v) {
        return v.getValue();
      }

      CODI_INLINE static Identifier getIdentifier(Type const& v) {
        return v.getIdentifier();
      }

      CODI_INLINE static void setValue(Type& v, Real const& value) {
        v.setValue(value);
      }

      CODI_INLINE static Real registerExternalFunctionOutput(Type& v) {
        return Type::getGlobalTape().registerExternalFunctionOutput(v);
      }

      template<typename Real, typename Identifier>
      static VectorWrapper* createVectorInterfaceWrapper(VectorAccessInterface<Real, Identifier>* accessInterface) {
        return accessInterface;
      }

      static void destroyVectorInterfaceWrapper(VectorWrapper* wrapper) {
        CODI_UNUSED(wrapper);

        // Did not create anything.
      }
  };
#endif
}
