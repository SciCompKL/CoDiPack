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

#include <cmath>
#include <complex>

#include "../config.h"
#include "../misc/macros.hpp"
#include "expressionTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Traits for values that can be used as real values e.g. double, float, codi::RealReverse etc..
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

        using Type = CODI_DD(T_Type, double);  ///< See TraitsImplementation

        using Real = Type;         ///< Inner type of the real value.
        using PassiveReal = Type;  ///< The original computation type, that was used in the application.

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

        using Type = CODI_DD(T_Type, double);  ///< See IsTotalFinite

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

        using Type = CODI_DD(T_Type, double);  ///< See IsTotalZero

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
     * accessor methods in this class access each of these value. For `getValue`, for example, a complex type of the
     * CoDiPack type's inner value is generated.
     *
     * @tparam T_Type  Any type that contains a CoDiPack type.
     */
    template<typename T_Type, typename = void>
    struct DataExtraction {
      public:
        CODI_STATIC_ASSERT(false && std::is_void<T_Type>::value,
                           "Instantiation of unspecialized RealTraits::DataExtraction.");

        using Type = CODI_DD(T_Type, ActiveTypeProxy);  ///< See DataExtraction.

        using Real = typename Type::Real;  ///< Type of primal values extracted from the type with AD values.
        using Identifier =
            typename Type::Identifier;  ///< Type of identifier values extracted from the type with AD values.

        /// Extract the primal values from a type of aggregated active types.
        CODI_INLINE static Real getValue(Type const& v);

        /// Extract the identifiers from a type of aggregated active types.
        CODI_INLINE static Identifier getIdentifier(Type const& v);

        /// Set the primal values of a type of aggregated active types.
        CODI_INLINE static void setValue(Type& v, Real const& value);
    };

    /**
     * @brief Tape registration methods for aggregated types that contain CoDiPack active types.
     *
     * An aggregated type is for example std::complex<codi::RealReverse>, which contains two CoDiPack values. The
     * methods in this class access each of these values in order to register the active types. For `registerInput`, the
     * real and imaginary part of the complex type are registered.
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

        /// Register all active types of a aggregated type as tape input.
        CODI_INLINE static void registerInput(Type& v);

        /// Register all active types of a aggregated type as tape output.
        CODI_INLINE static void registerOutput(Type& v);

        /// Register all active types of a aggregated type as external function outputs.
        CODI_INLINE static Real registerExternalFunctionOutput(Type& v);
    };

    /// \copydoc codi::RealTraits::DataExtraction::getValue()
    template<typename Type>
    typename DataExtraction<Type>::Real getValue(Type const& v) {
      return DataExtraction<Type>::getValue(v);
    }

    /// \copydoc codi::RealTraits::DataExtraction::getIdentifier()
    template<typename Type>
    typename DataExtraction<Type>::Identifier getIdentifier(Type const& v) {
      return DataExtraction<Type>::getIdentifier(v);
    }

    /// \copydoc codi::RealTraits::DataExtraction::setValue()
    template<typename Type>
    void setValue(Type& v, typename DataExtraction<Type>::Real const& value) {
      return DataExtraction<Type>::setValue(v, value);
    }

    /// \copydoc codi::RealTraits::TapeRegistration::registerInput()
    template<typename Type>
    void registerInput(Type& v) {
      return TapeRegistration<Type>::registerInput(v);
    }

    /// \copydoc codi::RealTraits::TapeRegistration::registerOutput()
    template<typename Type>
    void registerOutput(Type& v) {
      return TapeRegistration<Type>::registerOutput(v);
    }

    /// \copydoc codi::RealTraits::TapeRegistration::registerExternalFunctionOutput()
    template<typename Type>
    typename DataExtraction<Type>::Identifier registerExternalFunctionOutput(Type& v) {
      return TapeRegistration<Type>::registerExternalFunctionOutput(v);
    }

#ifndef DOXYGEN_DISABLE

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
    };

    /// Specialization of DataExtraction for complex types.
    template<typename T_InnerType>
    struct DataExtraction<std::complex<T_InnerType>> {
      public:

        using InnerType = CODI_DD(T_InnerType, ActiveTypeProxy);
        using Type = std::complex<T_InnerType>;  ///< See DataExtraction.

        using InnerExtraction = DataExtraction<InnerType>;

        using Real = std::complex<typename InnerExtraction::Real>;              ///< See DataExtraction::Real.
        using Identifier = std::complex<typename InnerExtraction::Identifier>;  ///< See DataExtraction::Identifier.

        /// \copydoc DataExtraction::getValue()
        CODI_INLINE static Real getValue(Type const& v) {
          InnerType const* vArray = reinterpret_cast<InnerType const*>(&v);

          return Real(InnerExtraction::getValue(vArray[0]), InnerExtraction::getValue(vArray[1]));
        }

        /// \copydoc DataExtraction::getIdentifier()
        CODI_INLINE static Identifier getIdentifier(Type const& v) {
          InnerType const* vArray = reinterpret_cast<InnerType const*>(&v);

          return Real(InnerExtraction::getIdentifier(vArray[0]), InnerExtraction::getIdentifier(vArray[1]));
        }

        /// \copydoc DataExtraction::setValue()
        CODI_INLINE static void setValue(Type& v, Real const& value) {
          InnerType* vArray = reinterpret_cast<InnerType*>(&v);

          InnerExtraction::setValue(vArray[0], std::real(value));
          InnerExtraction::setValue(vArray[1], std::imag(value));
        }
    };

    /// Specialization of TapeRegistration for complex types.
    template<typename T_InnerType>
    struct TapeRegistration<std::complex<T_InnerType>> {
      public:

        using InnerType = CODI_DD(T_InnerType, CODI_ANY);
        using Type = std::complex<T_InnerType>;  ///< See DataExtraction.

        using InnerRegistration = TapeRegistration<InnerType>;

        using Real = typename DataExtraction<Type>::Real;  ///< See DataExtraction::Real.

        /// \copydoc TapeRegistration::registerInput()
        CODI_INLINE static void registerInput(Type& v) {
          InnerType* vArray = reinterpret_cast<InnerType*>(&v);

          InnerRegistration::registerInput(vArray[0]);
          InnerRegistration::registerInput(vArray[1]);
        }

        /// \copydoc TapeRegistration::registerOutput()
        CODI_INLINE static void registerOutput(Type& v) {
          InnerType* vArray = reinterpret_cast<InnerType*>(&v);

          InnerRegistration::registerOutput(vArray[0]);
          InnerRegistration::registerOutput(vArray[1]);
        }

        /// \copydoc TapeRegistration::registerExternalFunctionOutput()
        CODI_INLINE static Real registerExternalFunctionOutput(Type& v) {
          InnerType* vArray = reinterpret_cast<InnerType*>(&v);

          return Real(InnerRegistration::registerExternalFunctionOutput(vArray[0]),
                      InnerRegistration::registerExternalFunctionOutput(vArray[1]));
        }
    };
#endif

    /// @}
    /*******************************************************************************/
    /// @name Detection of specific real value types
    /// @{

    /// If the real type is not handled by CoDiPack
    template<typename Type>
    using IsPassiveReal = std::is_same<Type, PassiveReal<Type>>;

#if CODI_IS_CPP14
    /// Value entry of IsPassiveReal
    template<typename Expr>
    bool constexpr isPassiveReal = IsPassiveReal<Expr>::value;
#endif

    /// Negated enable if wrapper for IsPassiveReal
    template<typename Type>
    using EnableIfNotPassiveReal = typename std::enable_if<!IsPassiveReal<Type>::value>::type;

    /// Enable if wrapper for IsPassiveReal
    template<typename Type>
    using EnableIfPassiveReal = typename std::enable_if<IsPassiveReal<Type>::value>::type;

    /// @}

  }
}
