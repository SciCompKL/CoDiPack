/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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
#include "../../expressions/aggregate/aggregatedActiveType.hpp"
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
   * (documentation/examples/Example_20_Aggregated_active_type_handling_in_external_functions.cpp): \snippet
   * examples/Example_20_Aggregated_active_type_handling_in_external_functions.cpp Typed external function
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

  /// Specialization of AggregatedTypeVectorAccessWrapper for aggregated types.
  ///
  /// @tparam The aggregated data type.
  template<typename T_Type>
  struct AggregatedTypeVectorAccessWrapper<T_Type, RealTraits::EnableIfAggregatedActiveType<T_Type>>
      : public VectorAccessInterface<typename RealTraits::DataExtraction<T_Type>::Real,
                                     typename RealTraits::DataExtraction<T_Type>::Identifier> {
    public:

      using Type = CODI_DD(
          T_Type,
          CODI_T(AggregatedActiveType<double, CODI_ANY, CODI_ANY>));  ///< See AggregatedTypeVectorAccessWrapper.
      using InnerType = typename Type::InnerActiveType;               ///< Active inner type of the aggregate.
      static int constexpr Elements = Type::Elements;                 ///< Number of elements in the aggregate.

      /// Inner interface expected by this wrapper.
      using InnerInterface = VectorAccessInterface<typename InnerType::Real, typename InnerType::Identifier>;

      using Real = typename RealTraits::DataExtraction<T_Type>::Real;  ///< Real value of the aggregate.
      using Identifier =
          typename RealTraits::DataExtraction<T_Type>::Identifier;  ///< Identifier value of the aggregate.

      using Traits = RealTraits::AggregatedTypeTraits<Real>;  ///< Aggregated type traits.

      InnerInterface& innerInterface;  ///< Reference to inner interface.
      int lhsOffset;  ///< Offset of indirect access if this aggregated type is part of an outer aggregated type.

      std::vector<Real> buffer;  ///< Temporary storage for getAdjointVec.

      /// Constructor
      CODI_INLINE AggregatedTypeVectorAccessWrapper(InnerInterface* innerInterface)
          : innerInterface(*innerInterface), lhsOffset(0), buffer(innerInterface->getVectorSize()) {}

      /*******************************************************************************/
      /// @name Misc

      /// \copydoc VectorAccessInterface::getVectorSize()
      CODI_INLINE size_t getVectorSize() const {
        return innerInterface.getVectorSize();
      }

      /// \copydoc VectorAccessInterface::isLhsZero()
      CODI_INLINE bool isLhsZero() const {
        bool isZero = true;

        static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
          innerInterface.setActiveVariableForIndirectAccess(lhsOffset + i.value);
          isZero &= innerInterface.isLhsZero();
        });

        return isZero;
      }

      /// \copydoc codi::VectorAccessInterface::clone()
      VectorAccessInterface<Real, Identifier>* clone() const {
        return new AggregatedTypeVectorAccessWrapper(&innerInterface);
      }

      /*******************************************************************************/
      /// @name Indirect adjoint access

      /// \copydoc codi::VectorAccessInterface::setLhsAdjoint
      CODI_INLINE void setLhsAdjoint(Identifier const& index) {
        static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
          innerInterface.setActiveVariableForIndirectAccess(lhsOffset + i.value);
          innerInterface.setLhsAdjoint(index[i.value]);
        });
      }

      /// \copydoc codi::VectorAccessInterface::updateAdjointWithLhs
      CODI_INLINE void updateAdjointWithLhs(Identifier const& index, Real const& jacobian) {
        static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
          innerInterface.setActiveVariableForIndirectAccess(lhsOffset + i.value);
          innerInterface.updateAdjointWithLhs(index[i.value], Traits::template arrayAccess<i.value>(jacobian));
        });
      }

      /*******************************************************************************/
      /// @name Indirect tangent access

      /// \copydoc codi::VectorAccessInterface::setLhsTangent
      CODI_INLINE void setLhsTangent(Identifier const& index) {
        static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
          innerInterface.setActiveVariableForIndirectAccess(lhsOffset + i.value);
          innerInterface.setLhsTangent(index[i.value]);
        });
      }

      /// \copydoc codi::VectorAccessInterface::updateTangentWithLhs
      CODI_INLINE void updateTangentWithLhs(Identifier const& index, Real const& jacobian) {
        static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
          innerInterface.setActiveVariableForIndirectAccess(lhsOffset + i.value);
          innerInterface.updateTangentWithLhs(index[i.value], Traits::template arrayAccess<i.value>(jacobian));
        });
      }

      /*******************************************************************************/
      /// @name Indirect adjoint/tangent access for functions with multiple outputs

      /// \copydoc VectorAccessInterface::setActiveVariableForIndirectAccess()
      CODI_INLINE void setActiveVariableForIndirectAccess(size_t pos) {
        lhsOffset = pos * Elements;
      }

      /*******************************************************************************/
      /// @name Direct adjoint access

      /// \copydoc VectorAccessInterface::resetAdjoint()
      CODI_INLINE void resetAdjoint(Identifier const& index, size_t dim) {
        static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
          innerInterface.resetAdjoint(index[i.value], dim);
        });
      }

      /// \copydoc VectorAccessInterface::resetAdjointVec()
      CODI_INLINE void resetAdjointVec(Identifier const& index) {
        static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
          innerInterface.resetAdjointVec(index[i.value]);
        });
      }

      /// \copydoc VectorAccessInterface::getAdjoint()
      CODI_INLINE Real getAdjoint(Identifier const& index, size_t dim) {
        Real adjoint{};

        static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
          Traits::template arrayAccess<i.value>(adjoint) = innerInterface.getAdjoint(index[i.value], dim);
        });

        return adjoint;
      }

      /// \copydoc VectorAccessInterface::getAdjointVec()
      CODI_INLINE void getAdjointVec(Identifier const& index, Real* const vec) {
        for (size_t curDim = 0; curDim < innerInterface.getVectorSize(); curDim += 1) {
          vec[curDim] = this->getAdjoint(index, curDim);
        }
      }

      /// \copydoc VectorAccessInterface::getAdjointVec()
      Real const* getAdjointVec(Identifier const& index) {
        getAdjointVec(index, buffer.data());
        return buffer.data();
      }

      /// \copydoc VectorAccessInterface::updateAdjoint()
      CODI_INLINE void updateAdjoint(Identifier const& index, size_t dim, Real const& adjoint) {
        static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
          innerInterface.updateAdjoint(index[i.value], dim, Traits::template arrayAccess<i.value>(adjoint));
        });
      }

      /// \copydoc VectorAccessInterface::updateAdjointVec()
      CODI_INLINE void updateAdjointVec(Identifier const& index, Real const* const vec) {
        for (size_t curDim = 0; curDim < innerInterface.getVectorSize(); curDim += 1) {
          this->updateAdjoint(index, curDim, vec[curDim]);
        }
      }

      /*******************************************************************************/
      /// @name Primal access

      /// \copydoc VectorAccessInterface::hasPrimals()
      CODI_INLINE bool hasPrimals() {
        return innerInterface.hasPrimals();
      }

      /// \copydoc VectorAccessInterface::setPrimal()
      CODI_INLINE void setPrimal(Identifier const& index, Real const& primal) {
        static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
          innerInterface.setPrimal(index[i.value], Traits::template arrayAccess<i.value>(primal));
        });
      }

      /// \copydoc VectorAccessInterface::getPrimal()
      CODI_INLINE Real getPrimal(Identifier const& index) {
        Real primal{};

        static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
          Traits::template arrayAccess<i.value>(primal) = innerInterface.getPrimal(index[i.value]);
        });

        return primal;
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
