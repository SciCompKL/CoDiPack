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

#include <array>
#include <cmath>
#include <complex>

#include "../config.h"
#include "../misc/compileTimeLoop.hpp"
#include "../misc/macros.hpp"
#include "computationTraits.hpp"
#include "expressionTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename T_Real, typename T_InnerActiveType, typename T_Impl, bool T_isStatic>
  struct AggregatedActiveTypeBase;

  /// Traits for values that can be used as real values, e.g. double, float, codi::RealReverse etc..
  namespace RealTraits {

    /*******************************************************************************/
    /// @name General real value traits
    /// @{

    /**
     * @brief Common traits for all types used as real values
     *
     * @tparam T_Type  The type of the real value.
     */
    template<typename T_Type, typename = void>
    struct TraitsImplementation {
      public:

        using Type = CODI_DD(T_Type, double);  ///< See TraitsImplementation.

        using Real = Type;         ///< Inner type of the real value.
        using PassiveReal = Type;  ///< The original computation type that was used in the application.

        static int constexpr MaxDerivativeOrder = 0;  ///< CoDiPack derivative order of the type.

        /// Get the basic primal value of the type.
        static CODI_INLINE PassiveReal const& getPassiveValue(Type const& v) {
          return v;
        }
    };

    /**
     * @brief Function for checking if all values of the type are finite.
     *
     * @tparam T_Type  The type of the real value.
     */
    template<typename T_Type, typename = void>
    struct IsTotalFinite {
      public:

        using Type = CODI_DD(T_Type, double);  ///< See IsTotalFinite.

        /// Checks if the values are all finite.
        static CODI_INLINE bool isTotalFinite(Type const& v) {
          using std::isfinite;
          return isfinite(v);
        }
    };

    /**
     * @brief Function for checking if the value of the type is completely zero.
     *
     * @tparam T_Type  The type of the real value.
     */
    template<typename T_Type, typename = void>
    struct IsTotalZero {
      public:

        using Type = CODI_DD(T_Type, double);  ///< See IsTotalZero.

        /// Checks if the values are completely zero.
        static CODI_INLINE bool isTotalZero(Type const& v) {
          return Type() == v;
        }
    };

    /// \copydoc codi::RealTraits::TraitsImplementation::Real
    template<typename Type>
    using Real = typename TraitsImplementation<Type>::Real;

    /// \copydoc codi::RealTraits::TraitsImplementation::PassiveReal
    template<typename Type>
    using PassiveReal = typename TraitsImplementation<Type>::PassiveReal;

    /// \copydoc codi::RealTraits::TraitsImplementation::MaxDerivativeOrder
    template<typename Type>
    CODI_INLINE size_t constexpr MaxDerivativeOrder() {
      return TraitsImplementation<Type>::MaxDerivativeOrder;
    }

    /// \copydoc codi::RealTraits::TraitsImplementation::getPassiveValue()
    template<typename Type>
    CODI_INLINE PassiveReal<Type> const& getPassiveValue(Type const& v) {
      return TraitsImplementation<Type>::getPassiveValue(v);
    }

    /// \copydoc codi::RealTraits::IsTotalFinite
    template<typename Type>
    CODI_INLINE bool isTotalFinite(Type const& v) {
      return IsTotalFinite<Type>::isTotalFinite(v);
    }

    /// \copydoc codi::RealTraits::IsTotalZero
    template<typename Type>
    CODI_INLINE bool isTotalZero(Type const& v) {
      return IsTotalZero<Type>::isTotalZero(v);
    }

    /// @}
    /*******************************************************************************/
    /// @name Traits for generalized data extraction
    /// @{

    /**
     * @brief Data handling methods for aggregated types that contain CoDiPack active types.
     *
     * An aggregated type is for example std::complex<codi::RealReverse>, which contains two CoDiPack values. The
     * accessor methods in this class make the corresponding aggregates of primal values and identifiers available,
     * that is, std::complex<double> for the primal values and std::complex<int> for the identifiers.
     *
     * @tparam T_Type  Any type that contains a CoDiPack type.
     */
    template<typename T_Type, typename = void>
    struct DataExtraction {
      public:
        CODI_STATIC_ASSERT(false && std::is_void<T_Type>::value,
                           "Instantiation of unspecialized RealTraits::DataExtraction.");

        using Type = CODI_DD(T_Type, ActiveTypeProxy);  ///< See DataExtraction.

        using Real = typename Type::Real;              ///< Type of a corresponding aggregate of primal values.
        using Identifier = typename Type::Identifier;  ///< Type of a corresponding aggregate of identifiers.

        /// Extract an aggregate of primal values from an aggregate of active types.
        CODI_INLINE static Real getValue(Type const& v);

        /// Extract an aggregate of identifiers from an aggregate of active types.
        CODI_INLINE static Identifier getIdentifier(Type const& v);

        /// Set the primal values of an aggregate of active types.
        CODI_INLINE static void setValue(Type& v, Real const& value);

        /// Set the identifiers of a type of aggregated active types.
        CODI_INLINE static void setIdentifier(Type& v, Identifier const& identifier);
    };

    /**
     * @brief Tape registration methods for aggregated types that contain CoDiPack active types.
     *
     * An aggregated type is for example std::complex<codi::RealReverse>, which contains two CoDiPack values. The
     * methods in this class take each of these values into account. For `registerInput`, for example, both the real and
     * imaginary part of the complex type are registered.
     *
     * @tparam T_Type  Any type that contains a CoDiPack type.
     */
    template<typename T_Type, typename = void>
    struct TapeRegistration {
      public:
        CODI_STATIC_ASSERT(false && std::is_void<T_Type>::value,
                           "Instantiation of unspecialized RealTraits::TapeRegistration.");

        using Type = CODI_DD(T_Type, CODI_ANY);  ///< See TapeRegistration.

        using Real = typename DataExtraction<Type>::Real;  ///< See DataExtraction::Real.

        /// Register all active types of an aggregated type as tape inputs.
        CODI_INLINE static void registerInput(Type& v);

        /// Register all active types of an aggregated type as tape outputs.
        CODI_INLINE static void registerOutput(Type& v);

        /// Register all active types of an aggregated type as external function outputs.
        CODI_INLINE static Real registerExternalFunctionOutput(Type& v);
    };

    /**
     * @brief Methods that access inner values of aggregated types that contain CoDiPack active types.
     *
     * An aggregated type is for example std::complex<codi::RealReverse>, which contains two CoDiPack values. The
     * methods in this class access each of these values. The real part is the element 0 and the imaginary part is the
     * element 1.
     *
     * @tparam T_Type  Any type that contains a CoDiPack type. (E.g. codi::RealReverse, std::complex<codi::RealReverse>)
     */
    template<typename T_Type, typename = void>
    struct AggregatedTypeTraits {
      public:
        using Type = CODI_DD(T_Type, double);  ///< See AggregatedTypeTraits.
        using InnerType = CODI_ANY;            ///< Inner type of the aggregated type.
        using Real = CODI_ANY;  ///< Real version of the aggregated type without the active CoDiPack types.

        static int constexpr Elements = 0;  ///< Number of elements of the aggregated type.

        /// Array construction of the aggregated type. It is defined as
        /// \f$ w = T(v_0, v_1, ..., v_{N-1}) \f$ where \f$ N \f$  is the number of elements.
        static Type arrayConstructor(InnerType const* v) {
          CODI_UNUSED(v);
          CODI_STATIC_ASSERT(false && std::is_void<Type>::value,
                             "Instantiation of unspecialized AggregatedTypeTraits.");

          return Type{};
        }

        /// Adjoint implementation of element-wise construction, that is, T is our aggregated type and the construction
        /// is defined as w = T(v_0, v_1, ..., v_{N-1}) where N is the number of elements. Then this function needs to
        /// implement the adjoint formulation of this construction, which is defined as
        /// \f$ \bar v_{element} = dT/dv_{element}^T * \bar w \f$ .
        template<size_t element>
        static InnerType adjointOfConstructor(Type const& w, Type const& w_b) {
          CODI_UNUSED(w, w_b);
          CODI_STATIC_ASSERT(false && std::is_void<Type>::value,
                             "Instantiation of unspecialized AggregatedTypeTraits.");

          return InnerType{};
        }

        /// Implementation of the array access, which is defined as \f$ v = w[element] \f$ where \f$ w \f$ is an
        /// aggregated type.
        template<size_t element>
        static InnerType arrayAccess(Type const& w) {
          CODI_UNUSED(w);
          CODI_STATIC_ASSERT(false && std::is_void<Type>::value,
                             "Instantiation of unspecialized AggregatedTypeTraits.");

          return InnerType{};
        }

        /// \copydoc AggregatedTypeTraits::arrayAccess()
        template<size_t element>
        CODI_INLINE static InnerType& arrayAccess(Type& v) {
          CODI_UNUSED(v);
          CODI_STATIC_ASSERT(false && std::is_void<Type>::value,
                             "Instantiation of unspecialized AggregatedTypeTraits.");

          static InnerType temp = {};
          return temp;
        }

        /// Implementation of the adjoint array access. See arrayAccess for the equation definition. The adjoint is
        /// defined as \f$ \bar w += dw[element]/w^T * \bar v \f$.
        template<size_t element>
        static Type adjointOfArrayAccess(Type const& w, InnerType const& v_b) {
          CODI_UNUSED(w, v_b);
          CODI_STATIC_ASSERT(false && std::is_void<Type>::value,
                             "Instantiation of unspecialized AggregatedTypeTraits.");

          return Type{};
        }

        /// \copydoc codi::ComputeOperation::getMathRep
        static CODI_INLINE std::string getMathRep() {
          CODI_STATIC_ASSERT(false && std::is_void<Type>::value,
                             "Instantiation of unspecialized AggregatedTypeTraits.");

          return "";
        }
    };

    /// Base implementation of AggregatedTypeTraits for aggregated types that can be defined as an array.
    template<typename T_Type, typename T_InnerType, typename T_Real, int T_Elements>
    struct ArrayAggregatedTypeTraitsBase {
      public:
        static_assert(
            sizeof(T_Type) == T_Elements * sizeof(T_InnerType),
            "ArrayAggregatedTypeTraitsBase is designed for aggregated types that resemble arrays of their inner type."
            " The sizes in the given instantiation do not match this use case.");

        using Type = CODI_DD(T_Type, CODI_ANY);            ///< See AggregatedTypeTraits.
        using InnerType = CODI_DD(T_InnerType, CODI_ANY);  ///< See AggregatedTypeTraits.
        using Real = CODI_DD(T_Real, CODI_ANY);            ///< See AggregatedTypeTraits.

        static int constexpr Elements = CODI_DD(T_Elements, 1);  ///< See AggregatedTypeTraits.

        /// \copydoc codi::RealTraits::AggregatedTypeTraits::arrayConstructor()
        CODI_INLINE static Type arrayConstructor(InnerType const* v) {
          Type w{};

          InnerType* wArray = reinterpret_cast<InnerType*>(&w);

          static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
            wArray[i.value] = v[i.value];
          });

          return w;
        }

        /// \copydoc codi::RealTraits::AggregatedTypeTraits::adjointOfConstructor()
        template<size_t element>
        CODI_INLINE static InnerType adjointOfConstructor(Type const& w, Type const& w_b) {
          CODI_UNUSED(w);

          // We assume that dT/dv is the identity and is not holomorphic. That means no one is changed to a minus one on
          // a transpose.
          InnerType const* w_bArray = reinterpret_cast<InnerType const*>(&w_b);
          return ComputationTraits::transpose(w_bArray[element]);
        }

        /// \copydoc codi::RealTraits::AggregatedTypeTraits::arrayAccess()
        template<size_t element>
        CODI_INLINE static InnerType& arrayAccess(Type& w) {
          InnerType* wArray = reinterpret_cast<InnerType*>(&w);

          return wArray[element];
        }

        /// \copydoc codi::RealTraits::AggregatedTypeTraits::arrayAccess()
        template<size_t element>
        CODI_INLINE static InnerType const& arrayAccess(Type const& w) {
          InnerType const* wArray = reinterpret_cast<InnerType const*>(&w);

          return wArray[element];
        }

        /// \copydoc codi::RealTraits::AggregatedTypeTraits::adjointOfArrayAccess()
        template<size_t element>
        CODI_INLINE static Type adjointOfArrayAccess(Type const& w, InnerType const& v_b) {
          CODI_UNUSED(w);

          // We assume that dw[element]/w is the identity and is not holomorphic. That means no one is changed to a
          // minus one on a transpose.
          Type w_b{};
          InnerType* w_bArray = reinterpret_cast<InnerType*>(&w_b);
          static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
            w_bArray[i.value] = InnerType();
          });
          w_bArray[element] = v_b;

          return w_b;
        }
    };

    /// \copydoc codi::RealTraits::DataExtraction::getValue()
    template<typename Type>
    CODI_INLINE typename DataExtraction<Type>::Real getValue(Type const& v) {
      return DataExtraction<Type>::getValue(v);
    }

    /// \copydoc codi::RealTraits::DataExtraction::getIdentifier()
    template<typename Type>
    CODI_INLINE typename DataExtraction<Type>::Identifier getIdentifier(Type const& v) {
      return DataExtraction<Type>::getIdentifier(v);
    }

    /// \copydoc codi::RealTraits::DataExtraction::setValue()
    template<typename Type>
    CODI_INLINE void setValue(Type& v, typename DataExtraction<Type>::Real const& value) {
      return DataExtraction<Type>::setValue(v, value);
    }

    /// \copydoc codi::RealTraits::TapeRegistration::registerInput()
    template<typename Type>
    CODI_INLINE void registerInput(Type& v) {
      return TapeRegistration<Type>::registerInput(v);
    }

    /// \copydoc codi::RealTraits::TapeRegistration::registerOutput()
    template<typename Type>
    CODI_INLINE void registerOutput(Type& v) {
      return TapeRegistration<Type>::registerOutput(v);
    }

    /// \copydoc codi::RealTraits::TapeRegistration::registerExternalFunctionOutput()
    template<typename Type>
    CODI_INLINE typename DataExtraction<Type>::Real registerExternalFunctionOutput(Type& v) {
      return TapeRegistration<Type>::registerExternalFunctionOutput(v);
    }

    /// @}
    /*******************************************************************************/
    /// @name Detection of specific real value types
    /// @{

    /// Enable if helper when a type has been specialized for the AggregatedTypeTraits.
    template<typename Type>
    using EnableIfAggregatedTypeTraitsIsSpecialized =
        typename std::enable_if<(AggregatedTypeTraits<Type>::Elements != 0) &
                                (!ExpressionTraits::IsLhsExpression<Type>::value)>::type;

    /// Enable if helper for AggregatedActiveType.
    template<typename Expr>
    using EnableIfAggregatedActiveType =
        typename enable_if_base_of<AggregatedActiveTypeBase<typename Expr::Real, typename Expr::InnerActiveType,
                                                            typename Expr::Impl, Expr::isStatic>,
                                   Expr>::type;

    /// If the real type is not handled by CoDiPack
    template<typename Type>
    using IsPassiveReal = std::is_same<Type, PassiveReal<Type>>;

    /// Value entry of IsPassiveReal
    template<typename Expr>
    bool constexpr isPassiveReal = IsPassiveReal<Expr>::value;

    /// Negated enable if wrapper for IsPassiveReal
    template<typename Type>
    using EnableIfNotPassiveReal = typename std::enable_if<!IsPassiveReal<Type>::value>::type;

    /// Enable if wrapper for IsPassiveReal
    template<typename Type>
    using EnableIfPassiveReal = typename std::enable_if<IsPassiveReal<Type>::value>::type;

#ifndef DOXYGEN_DISABLE

    /// @}
    /*******************************************************************************/
    /// @name Various specializations.
    /// @{

    /// Specialization of DataExtraction for floating point types.
    template<typename T_Type>
    struct DataExtraction<T_Type, typename std::enable_if<std::is_floating_point<T_Type>::value>::type> {
      public:
        using Type = CODI_DD(T_Type, double);  ///< See DataExtraction.

        using Real = double;     ///< See DataExtraction::Real.
        using Identifier = int;  ///< See DataExtraction::Identifier.

        /// \copydoc DataExtraction::getValue()
        CODI_INLINE static Real getValue(Type const& v) {
          return v;
        }

        /// \copydoc DataExtraction::getIdentifier()
        CODI_INLINE static Identifier getIdentifier(Type const& v) {
          CODI_UNUSED(v);

          return 0;
        }

        /// \copydoc DataExtraction::setValue()
        CODI_INLINE static void setValue(Type& v, Real const& value) {
          v = value;
        }

        /// \copydoc DataExtraction::setIdentifier()
        CODI_INLINE static void setIdentifier(Type& v, Identifier const& identifier) {
          CODI_UNUSED(v, identifier);
        }
    };

    /// Specialization of DataExtraction for aggregated expression types.
    template<typename T_Type>
    struct DataExtraction<T_Type, EnableIfAggregatedTypeTraitsIsSpecialized<T_Type>> {
      public:

        using Type = CODI_DD(T_Type, double);

        using TypeTraits = AggregatedTypeTraits<Type>;
        using InnerType = typename TypeTraits::InnerType;
        static int constexpr Elements = TypeTraits::Elements;

        using InnerDataExtraction = DataExtraction<InnerType>;

        using Real = typename TypeTraits::Real;
        using Identifier = std::array<CODI_DD(typename InnerType::Identifier, int), Elements>;

        using RealTypeTraits = AggregatedTypeTraits<Real>;

        /// \copydoc DataExtraction::getValue()
        CODI_INLINE static Real getValue(Type const& v) {
          Real real{};

          static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
            RealTypeTraits::template arrayAccess<i.value>(real) =
                InnerDataExtraction::getValue(TypeTraits::template arrayAccess<i.value>(v));
          });

          return real;
        }

        /// \copydoc DataExtraction::getIdentifier()
        CODI_INLINE static Identifier getIdentifier(Type const& v) {
          Identifier res;

          static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
            res[i.value] = InnerDataExtraction::getIdentifier(TypeTraits::template arrayAccess<i.value>(v));
          });

          return res;
        }

        /// \copydoc DataExtraction::setValue()
        CODI_INLINE static void setValue(Type& v, Real const& value) {
          static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
            InnerDataExtraction::setValue(TypeTraits::template arrayAccess<i.value>(v),
                                          RealTypeTraits::template arrayAccess<i.value>(value));
          });
        }

        /// \copydoc DataExtraction::setIdentifier()
        CODI_INLINE static void setIdentifier(Type& v, Identifier const& identifier) {
          static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
            InnerDataExtraction::setIdentifier(TypeTraits::template arrayAccess<i.value>(v), identifier[i.value]);
          });
        }
    };

    /// Specialization of TapeRegistration for aggregated expression types.
    template<typename T_Type>
    struct TapeRegistration<T_Type, EnableIfAggregatedTypeTraitsIsSpecialized<T_Type>> {
      public:

        using Type = CODI_DD(T_Type, CODI_ANY);

        using TypeTraits = AggregatedTypeTraits<Type>;
        using InnerType = typename TypeTraits::InnerType;
        static int constexpr Elements = TypeTraits::Elements;

        using Real = typename TypeTraits::Real;

        using RealTypeTraits = RealTraits::AggregatedTypeTraits<Real>;
        using InnerTraits = TapeRegistration<InnerType>;

        /// \copydoc TapeRegistration::registerInput()
        CODI_INLINE static void registerInput(Type& v) {
          static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
            InnerTraits::registerInput(TypeTraits::template arrayAccess<i.value>(v));
          });
        }

        /// \copydoc TapeRegistration::registerOutput()
        CODI_INLINE static void registerOutput(Type& v) {
          static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
            InnerTraits::registerOutput(TypeTraits::template arrayAccess<i.value>(v));
          });
        }

        /// \copydoc TapeRegistration::registerExternalFunctionOutput()
        CODI_INLINE static Real registerExternalFunctionOutput(Type& v) {
          Real res{};

          static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
            RealTypeTraits::template arrayAccess<i.value>(res) =
                InnerTraits::registerExternalFunctionOutput(TypeTraits::template arrayAccess<i.value>(v));
          });

          return res;
        }
    };
    template<typename T_Type>
    struct AggregatedTypeTraits<T_Type, typename std::enable_if<std::is_floating_point<T_Type>::value>::type>
        : ArrayAggregatedTypeTraitsBase<T_Type, T_Type, T_Type, 1> {};

    template<typename T_Type>
    struct AggregatedTypeTraits<T_Type, typename std::enable_if<std::is_integral<T_Type>::value>::type>
        : ArrayAggregatedTypeTraitsBase<T_Type, T_Type, T_Type, 1> {};

    template<>
    struct AggregatedTypeTraits<int*> : ArrayAggregatedTypeTraitsBase<int*, int*, int*, 1> {};

#endif

    /// @}

  }
}
