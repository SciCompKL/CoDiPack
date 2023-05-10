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

#include <complex>
#include <vector>

#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../misc/macros.hpp"
#include "../../tapes/misc/vectorAccessInterface.hpp"
#include "../../traits/computationTraits.hpp"
#include "../../traits/expressionTraits.hpp"
#include "../../traits/realTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Generalized wrapper of the VectorAccessInterface for aggregated data types,
   * e.g. std::complex<codi::RealReverse>.
   *
   * This wrapper is instantiated by AggregatedTypeVectorAccessWrapperFactory. It can be specialized for arbitrary
   * types that consist of CoDiPack types.
   *
   * This class helps to write generalized external function, that handles aggregated data types. E.g. for
   * std::complex<codi::RealReverse> the primal value as well as the adjoint value are defined as std::complex<double>
   * and the identifier type is defined as std::complex<int>. Since this can be different for other types like vectors
   * or matrices, this wrapper can be used to write generalized code, that works for arbitrary aggregated types. See
   * RealTraits::DataExtraction for a generalized access to the primal and identifier data of aggregated types.
   *
   * Here is an example for a generalized external function routine
   * (documentation/examples/Example_20_Aggregated_active_type_handling.cpp): \snippet
   * examples/Example_20_Aggregated_active_type_handling.cpp Typed external function
   *
   * In general all implementations of the wrapper will forward all functions calls to the VectorAccessInterface of
   * the underlying tape. For each active type in the structure the corresponding operation should be performed.
   *
   * Implementations can use AggregatedTypeVectorAccessWrapperBase which implements most of the functions from
   * VectorAccessInterface.
   *
   * @tparam T_Type  An arbitrary type which consists of multiple CoDiPack types.
   */
  template<typename T_Type, typename = void>
  struct AggregatedTypeVectorAccessWrapper : public VectorAccessInterface<CODI_ANY, CODI_ANY> {
      CODI_STATIC_ASSERT(false && std::is_void<T_Type>::value,
                         "Instantiation of unspecialized AggregatedTypeVectorAccessWrapper.");

      using Type = CODI_DD(T_Type, CODI_ANY);  ///< See AggregatedTypeVectorAccessWrapperBase.
  };

  /**
   * @brief Implements all methods from AggregatedTypeVectorAccessWrapper, that can be implemented with combinations of
   * other methods.
   *
   * @tparam T_Real  Primal value type of the combined type.
   * @tparam T_Identifier  Identifier type type of the combined type.
   * @tparam T_InnerInterface  The VectorAccessInterface of the underlying tape.
   */
  template<typename T_Real, typename T_Identifier, typename T_InnerInterface>
  struct AggregatedTypeVectorAccessWrapperBase : public VectorAccessInterface<T_Real, T_Identifier> {
    public:

      using Real = CODI_DD(T_Real, double);           ///< See RealTraits::DataExtraction::Real.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See RealTraits::DataExtraction::Identifier.

      using InnerInterface =
          CODI_DD(T_InnerInterface,
                  CODI_T(VectorAccessInterface<double, int>));  ///< See AggregatedTypeVectorAccessWrapperBase.

    protected:

      InnerInterface& innerInterface;  ///< Reference to the accessor of the underlying tape.

      std::vector<Real> lhs;     ///< Temporary storage for indirect adjoint or tangent updates.
      std::vector<Real> buffer;  ///< Temporary storage for getAdjointVec access.

    public:

      /// Constructor
      AggregatedTypeVectorAccessWrapperBase(InnerInterface* innerInterface)
          : innerInterface(*innerInterface),
            lhs(innerInterface->getVectorSize()),
            buffer(innerInterface->getVectorSize()) {}

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
      void updateAdjointWithLhs(Identifier const& index, Real const& jacobian) {
        for (size_t curDim = 0; curDim < lhs.size(); curDim += 1) {
          Real update = jacobian * lhs[curDim];
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
      void updateTangentWithLhs(Identifier const& index, Real const& jacobian) {
        for (size_t curDim = 0; curDim < lhs.size(); curDim += 1) {
          lhs[curDim] += jacobian * this->getAdjoint(index, curDim);
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

      /// \copydoc VectorAccessInterface::getAdjointVec()
      Real const* getAdjointVec(Identifier const& index) {
        getAdjointVec(index, buffer.data());
        return buffer.data();
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

  /**
   * @brief Factory for the creation of AggregatedTypeVectorAccessWrapper instances.
   *
   * This factory is specialized for CoDiPack types to return the provided interface, thus removing the overhead of
   * a wrapped interface.
   *
   * User can specialize this factory if the default construction of AggregatedTypeVectorAccessWrapper needs to be
   * specialized.
   *
   * @tparam T_Type See AggregatedTypeVectorAccessWrapper.
   */
  template<typename T_Type, typename = void>
  struct AggregatedTypeVectorAccessWrapperFactory {
    public:
      using Type = CODI_DD(T_Type, CODI_ANY);  ///< See AggregatedTypeVectorAccessWrapperBase.

      using RType = AggregatedTypeVectorAccessWrapper<Type>;  ///< Which instances this factory creates.

      /// Instantiate a AggregatedTypeVectorAccessWrapper class.
      ///
      /// @param access  The vector access interface from underlying tape.
      template<typename Real, typename Identifier>
      static RType* create(VectorAccessInterface<Real, Identifier>* access) {
        return new RType(access);
      }

      /// Delete the AggregatedTypeVectorAccessWrapper instance created by the crate method.
      static void destroy(RType* access) {
        delete access;
      }
  };

#ifndef DOXYGEN_DISABLE
  /// Specialization of AggregatedTypeVectorAccessWrapper for std::complex.
  ///
  /// @tparam The nested type of the complex data type.
  template<typename T_InnerType>
  struct AggregatedTypeVectorAccessWrapper<std::complex<T_InnerType>>
      : public AggregatedTypeVectorAccessWrapperBase<
            std::complex<typename T_InnerType::Real>, std::complex<typename T_InnerType::Identifier>,
            VectorAccessInterface<typename T_InnerType::Real, typename T_InnerType::Identifier>> {
    public:

      using InnerType = CODI_DD(T_InnerType, CODI_DEFAULT_LHS_EXPRESSION);  ///< See AggregatedTypeVectorAccessWrapper.
      using Type = std::complex<InnerType>;                                 ///< See AggregatedTypeVectorAccessWrapper.

      using InnerInterface = VectorAccessInterface<
          typename InnerType::Real,
          typename InnerType::Identifier>;  ///< See AggregatedTypeVectorAccessWrapperBase::InnerInterface.

      using Real = std::complex<typename InnerType::Real>;              ///< See RealTraits::DataExtraction::Real.
      using Identifier = std::complex<typename InnerType::Identifier>;  ///< See RealTraits::DataExtraction::Real.

      using Base =
          AggregatedTypeVectorAccessWrapperBase<Real, Identifier, InnerInterface>;  ///< Base class abbreviation.

      /// Constructor
      AggregatedTypeVectorAccessWrapper(InnerInterface* innerInterface) : Base(innerInterface) {}

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

      /// \copydoc VectorAccessInterface::clone()
      VectorAccessInterface<Real, Identifier>* clone() const {
        return new AggregatedTypeVectorAccessWrapper(Base::innerInterface.clone());
      }
  };

  /// Specialization of AggregatedTypeVectorAccessWrapperFactory for CoDiPack active types.
  ///
  /// @tparam T_Type  A CoDiPack active type.
  template<typename T_Type>
  struct AggregatedTypeVectorAccessWrapperFactory<T_Type, ExpressionTraits::EnableIfLhsExpression<T_Type>> {
    public:
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See AggregatedTypeVectorAccessWrapperBase.

      using RType = VectorAccessInterface<typename Type::Real, typename Type::Identifier>;

      /// \copydoc AggregatedTypeVectorAccessWrapperFactory::create()
      static RType* create(RType* access) {
        return access;
      }

      /// \copydoc AggregatedTypeVectorAccessWrapperFactory::destroy()
      static void destroy(RType* access) {
        CODI_UNUSED(access);

        // Do nothing
      }
  };
#endif
}
