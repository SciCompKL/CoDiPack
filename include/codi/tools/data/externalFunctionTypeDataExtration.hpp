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

  /**
   * @brief Interface extension of VectorAccessInterface to extract also primal values and identifiers from values, as
   * well as helper functions for external function handling.
   *
   * This interface is instantiated by ExternalFunctionTypeDataExtraction. It can be specialized for arbitrary types,
   * that consist of CoDiPack types.
   *
   * In general all implementations of the interface will forward all functions calls to the VectorAccessInterface of
   * the underlying tape. For each active type in the structure the corresponding operation should be performed.
   *
   * Implementations can use VectorAccessTypeWrapperBase which implements most of the functions from
   * VectorAccessInterface.
   *
   * @tparam _Type  An arbitrary type which consists of multiple CoDiPack types.
   */
  template<typename _Type, typename = void>
  struct VectorAccessTypeWrapper : public VectorAccessInterface<CODI_ANY, CODI_ANY> {
      static_assert(false && std::is_void<_Type>::value, "Instantiation of unspecialized VectorAccessTypeWrapper.");

      using Type = CODI_DD(_Type, CODI_ANY);  ///< See VectorAccessTypeWrapperBase.

      using Real = CODI_ANY;        ///< See ExternalFunctionTypeDataExtraction::Real.
      using Identifier = CODI_ANY;  ///< See ExternalFunctionTypeDataExtraction::Identifier.

      /*******************************************************************************/
      /// @name Static data extraction and registration functions

      /// \copydoc ExternalFunctionTypeDataExtraction::getValue()
      static Real getValue(Type const& v);

      /// \copydoc ExternalFunctionTypeDataExtraction::getIdentifier()
      static Identifier getIdentifier(Type const& v);

      /// \copydoc ExternalFunctionTypeDataExtraction::setValue()
      static void setValue(Type& v, Real const& value);

      /// \copydoc ExternalFunctionTypeDataExtraction::registerExternalFunctionOutput()
      static Real registerExternalFunctionOutput(Type& v);
  };

  /// @brief Implements all methods from VectorAccessTypeWrapper, that can be implemented as aggregates of other
  /// methods.
  ///
  /// @tparam _Real  Primal value type of the combined type.
  /// @tparam _Identifier  Identifier type type of the combined type.
  /// @tparam _InnerInterface  The VectorAccessInterface of the underlying tape.
  template<typename _Real, typename _Identifier, typename _InnerInterface>
  struct VectorAccessTypeWrapperBase : public VectorAccessInterface<_Real, _Identifier> {
    public:

      using Real = CODI_DD(_Real, CODI_ANY);              ///< See ExternalFunctionTypeDataExtraction::Real.
      using Identifier = CODI_DD(_Identifier, CODI_ANY);  ///< See ExternalFunctionTypeDataExtraction::Real.

      using InnerInterface = CODI_DD(_InnerInterface,
                                     CODI_T(VectorAccessInterface<double, int>));  ///< See VectorAccessTypeWrapperBase.

    protected:

      InnerInterface& innerInterface;  ///< Reference to the accessor of the underlying tape.

      std::vector<Real> lhs;  ///< Temporary storage for indirect adjoint or tangent updates.

    public:

      /// Constructor
      VectorAccessTypeWrapperBase(InnerInterface* innerInterface)
          : innerInterface(*innerInterface), lhs(innerInterface->getVectorSize()) {}

      /*******************************************************************************/
      /// @name Misc

      /// \copydoc VectorAccessInterface::getVectorSize()
      size_t getVectorSize() const {
        return innerInterface.getVectorSize();
      }

      /// \copydoc VectorAccessInterface::isLhsZero()
      bool isLhsZero() {
        bool isZero = true;

        for (size_t curDim = 0; isZero && curDim < lhs.size(); curDim += 1) {
          isZero &= RealTraits::isTotalZero(lhs[curDim]);
        }

        return isZero;
      }

      /*******************************************************************************/
      /// @name Indirect adjoint access

      /// \copydoc VectorAccessInterface::setLhsAdjoint()
      void setLhsAdjoint(Identifier const& index) {
        getAdjointVec(index, lhs.data());
        this->resetAdjointVec(index);
      }

      /// \copydoc VectorAccessInterface::updateAdjointWithLhs()
      void updateAdjointWithLhs(Identifier const& index, Real const& jacobi) {
        for (size_t curDim = 0; curDim < lhs.size(); curDim += 1) {
          Real update = jacobi * lhs[curDim];
          this->updateAdjoint(index, curDim, update);
        }
      }

      /*******************************************************************************/
      /// @name Indirect tangent access

      /// \copydoc VectorAccessInterface::setLhsTangent()
      void setLhsTangent(Identifier const& index) {
        updateAdjointVec(index, lhs.data());

        for (size_t curDim = 0; curDim < lhs.size(); curDim += 1) {
          lhs[curDim] = Real();
        }
      }

      /// \copydoc VectorAccessInterface::updateTangentWithLhs()
      void updateTangentWithLhs(Identifier const& index, Real const& jacobi) {
        for (size_t curDim = 0; curDim < lhs.size(); curDim += 1) {
          lhs[curDim] += jacobi * this->getAdjoint(index, curDim);
        }
      }

      /*******************************************************************************/
      /// @name Direct adjoint access

      /// \copydoc VectorAccessInterface::getAdjointVec()
      void getAdjointVec(Identifier const& index, Real* const vec) {
        for (size_t curDim = 0; curDim < lhs.size(); curDim += 1) {
          vec[curDim] = this->getAdjoint(index, curDim);
        }
      }

      /// \copydoc VectorAccessInterface::updateAdjointVec()
      void updateAdjointVec(Identifier const& index, Real const* const vec) {
        for (size_t curDim = 0; curDim < lhs.size(); curDim += 1) {
          this->updateAdjoint(index, curDim, vec[curDim]);
        }
      }

      /*******************************************************************************/
      /// @name Primal access

      /// \copydoc VectorAccessInterface::hasPrimals()
      bool hasPrimals() {
        return innerInterface.hasPrimals();
      }
  };

#ifndef DOXYGEN_DISABLE
  /// Specialization of VectorAccessTypeWrapper for std::complex.
  ///
  /// @tparam The nested type of the complex data type.
  template<typename _InnerType>
  struct VectorAccessTypeWrapper<std::complex<_InnerType>>
      : public VectorAccessTypeWrapperBase<
            std::complex<typename _InnerType::Real>, std::complex<typename _InnerType::Identifier>,
            VectorAccessInterface<typename _InnerType::Real, typename _InnerType::Identifier>> {
    public:

      using InnerType = CODI_DD(
          _InnerType,
          CODI_T(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));  ///< See VectorAccessTypeWrapper.
      using Type = std::complex<InnerType>;                                     ///< See VectorAccessTypeWrapper.

      using InnerInterface =
          VectorAccessInterface<typename InnerType::Real,
                                typename InnerType::Identifier>;  ///< See VectorAccessTypeWrapperBase::InnerInterface.

      using Real = std::complex<typename InnerType::Real>;  ///< See ExternalFunctionTypeDataExtraction::Real.
      using Identifier =
          std::complex<typename InnerType::Identifier>;  ///< See ExternalFunctionTypeDataExtraction::Identifier.

      using Base = VectorAccessTypeWrapperBase<Real, Identifier, InnerInterface>;  ///< Base class abbreviation.

      /// Constructor
      VectorAccessTypeWrapper(InnerInterface* innerInterface) : Base(innerInterface) {}

      /*******************************************************************************/
      /// @name Static data extraction and registration functions

      /// \copydoc VectorAccessTypeWrapper::getValue()
      CODI_INLINE static Real getValue(Type const& v) {
        return Real(std::real(v).getValue(), std::imag(v).getValue());
      }

      /// \copydoc VectorAccessTypeWrapper::getIdentifier()
      CODI_INLINE static Identifier getIdentifier(Type const& v) {
        InnerType const* vArray = reinterpret_cast<InnerType const*>(&v);
        return Real(vArray[0].getIdentifier(), vArray[1].getIdentifier());
      }

      /// \copydoc VectorAccessTypeWrapper::setValue()
      CODI_INLINE static void setValue(Type& v, Real const& value) {
        InnerType* vArray = reinterpret_cast<InnerType*>(&v);

        vArray[0].setValue(std::real(value));
        vArray[1].setValue(std::imag(value));
      }

      /// \copydoc VectorAccessTypeWrapper::registerExternalFunctionOutput()
      CODI_INLINE static Real registerExternalFunctionOutput(Type& v) {
        InnerType* vArray = reinterpret_cast<InnerType*>(&v);

        return Real(InnerType::getGlobalTape().registerExternalFunctionOutput(vArray[0]),
                    InnerType::getGlobalTape().registerExternalFunctionOutput(vArray[1]));
      }

      /*******************************************************************************/
      /// @name Direct adjoint access

      /// \copydoc VectorAccessInterface::resetAdjoint()
      void resetAdjoint(Identifier const& index, size_t dim) {
        Base::innerInterface.resetAdjoint(std::real(index), dim);
        Base::innerInterface.resetAdjoint(std::imag(index), dim);
      }

      /// \copydoc VectorAccessInterface::resetAdjointVec()
      void resetAdjointVec(Identifier const& index) {
        Base::innerInterface.resetAdjointVec(std::real(index));
        Base::innerInterface.resetAdjointVec(std::imag(index));
      }

      /// \copydoc VectorAccessInterface::getAdjoint()
      Real getAdjoint(Identifier const& index, size_t dim) {
        return Real(Base::innerInterface.getAdjoint(std::real(index), dim),
                    Base::innerInterface.getAdjoint(std::imag(index), dim));
      }

      /// \copydoc VectorAccessInterface::updateAdjoint()
      void updateAdjoint(Identifier const& index, size_t dim, Real const& adjoint) {
        Base::innerInterface.updateAdjoint(std::real(index), dim, std::real(adjoint));
        Base::innerInterface.updateAdjoint(std::imag(index), dim, std::imag(adjoint));
      }

      /*******************************************************************************/
      /// @name Primal access

      /// \copydoc VectorAccessInterface::setPrimal()
      void setPrimal(Identifier const& index, Real const& primal) {
        Base::innerInterface.setPrimal(std::real(index), std::real(primal));
        Base::innerInterface.setPrimal(std::imag(index), std::imag(primal));
      }

      /// \copydoc VectorAccessInterface::getPrimal()
      Real getPrimal(Identifier const& index) {
        return Real(Base::innerInterface.getPrimal(std::real(index)), Base::innerInterface.getPrimal(std::imag(index)));
      }
  };
#endif

  /**
   * @brief Data handling helper for the implementation of external functions. Generalizes the access to primal,
   * identifier and adjoint data for arbitrary types that have CoDiPack tapes as inner types.
   * (E.g. std::complex<codi::RealReverse>)
   *
   * This class helps to extract data from a data type. E.g. for std::complex<codi::RealReverse> the primal value is
   * std::complex<double> and the identifier type is std::complex<int> or int[2]. Since this can be different for other
   * types like vectors or matrices, this wrapper can be used to write generalized code, that works for arbitrary types.
   *
   * The logic for each type can be implemented by specializing either this class or VectorAccessTypeWrapper for the
   * required type. The default implementation uses VectorAccessTypeWrapper to define its logic.
   *
   * Here is an example for a generalized external function routine
   * (documentation/examples/manualStatementPushTapeInterface.cpp): \snippet
   * examples/Example_20_External_function_type_data_extraction.cpp Typed external function
   *
   * @tparam _Type  An arbitrary type which consists of multiple CoDiPack types.
   */
  template<typename _Type, typename = void>
  struct ExternalFunctionTypeDataExtraction {
      using Type = CODI_DD(_Type, CODI_ANY);  ///< See ExternalFunctionTypeDataExtraction.

      using VectorWrapper =
          VectorAccessTypeWrapper<Type>;  ///< Wrapper for the VectorAccessInterface of the underlying tape.

      using Real = typename VectorWrapper::Real;  ///< Type of primal values extracted from the type with AD values.
      using Identifier =
          typename VectorWrapper::Identifier;  ///< Type of identifier values extracted from the type with AD values.

      /// Extract the primal values from a type of combined active types.
      CODI_INLINE static Real getValue(Type const& v) {
        return VectorWrapper::getValue(v);
      }

      /// Extract the identifiers from a type of combined active types.
      CODI_INLINE static Identifier getIdentifier(Type const& v) {
        return VectorWrapper::getIdentifier(v);
      }

      /// Set the primal values of a type of combined active types.
      CODI_INLINE static void setValue(Type& v, Real const& value) {
        VectorWrapper::setValue(v, value);
      }

      /// Register all active types of a combined type as external function outputs.
      CODI_INLINE static Real registerExternalFunctionOutput(Type& v) {
        return VectorWrapper::registerExternalFunctionOutput(v);
      }

      /// Create a wrapper vector access interface from the vector access interface of the underlying tape.
      template<typename RealInterface, typename IdentifierInterface>
      static VectorWrapper* createVectorInterfaceWrapper(
          VectorAccessInterface<RealInterface, IdentifierInterface>* accessInterface) {
        return new VectorWrapper(accessInterface);
      }

      /// Destroy a created a wrapper vector access interface created by createVectorInterfaceWrapper.
      static void destroyVectorInterfaceWrapper(VectorWrapper* wrapper) {
        delete wrapper;
      }
  };

#ifndef DOXYGEN_DISABLE
  /// \copydoc ExternalFunctionTypeDataExtraction
  ///
  /// Specialization for CoDiPack active types. This implementation will not create any wrapper interfaces.
  ///
  /// @tparam _ActiveType  An CoDiPack active type.
  template<typename _ActiveType>
  struct ExternalFunctionTypeDataExtraction<_ActiveType, ExpressionTraits::EnableIfLhsExpression<_ActiveType>> {
      using Type = CODI_DD(
          _ActiveType,
          CODI_T(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));  ///< See
                                                                                ///< ExternalFunctionTypeDataExtraction.

      using VectorWrapper =
          VectorAccessInterface<typename Type::Real,
                                typename Type::Identifier>;  ///< See ExternalFunctionTypeDataExtraction.

      using Real = typename Type::Real;              ///< See ExternalFunctionTypeDataExtraction.
      using Identifier = typename Type::Identifier;  ///< See ExternalFunctionTypeDataExtraction.

      /// \copydoc ExternalFunctionTypeDataExtraction::getValue()
      CODI_INLINE static Real getValue(Type const& v) {
        return v.getValue();
      }

      /// \copydoc ExternalFunctionTypeDataExtraction::getIdentifier()
      CODI_INLINE static Identifier getIdentifier(Type const& v) {
        return v.getIdentifier();
      }

      /// \copydoc ExternalFunctionTypeDataExtraction::setValue()
      CODI_INLINE static void setValue(Type& v, Real const& value) {
        v.setValue(value);
      }

      /// \copydoc ExternalFunctionTypeDataExtraction::registerExternalFunctionOutput()
      CODI_INLINE static Real registerExternalFunctionOutput(Type& v) {
        return Type::getGlobalTape().registerExternalFunctionOutput(v);
      }

      /// \copydoc ExternalFunctionTypeDataExtraction::createVectorInterfaceWrapper()
      ///
      /// Will forward the argument as return.
      template<typename Real, typename Identifier>
      static VectorWrapper* createVectorInterfaceWrapper(VectorAccessInterface<Real, Identifier>* accessInterface) {
        return accessInterface;
      }

      /// \copydoc ExternalFunctionTypeDataExtraction::destroyVectorInterfaceWrapper();
      ///
      /// Will not delete the interface.
      static void destroyVectorInterfaceWrapper(VectorWrapper* wrapper) {
        CODI_UNUSED(wrapper);

        // Did not create anything.
      }
  };
#endif
}
