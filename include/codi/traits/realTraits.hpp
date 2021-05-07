#pragma once

#include <cmath>

#include "../aux/macros.hpp"
#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  namespace RealTraits {

    /*******************************************************************************/
    /// @name General real value traits
    /// @{

    /**
     * @brief Common traits for all types used as real values
     *
     * @tparam _Type  The type of the real value.
     */
    template<typename _Type, typename = void>
    struct TraitsImplementation {
      public:

        using Type = CODI_DD(_Type, double);  ///< See TraitsImplementation

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
     * @tparam _Type  The type of the real value.
     */
    template<typename _Type, typename = void>
    struct IsTotalFinite {
      public:

        using Type = CODI_DD(_Type, double);  ///< See IsTotalFinite

        /// Checks if the values are all finite.
        static CODI_INLINE bool isTotalFinite(Type const& v) {
          using std::isfinite;
          return isfinite(v);
        }
    };

    /**
     * @brief Function for checking if the value of the type is completely zero.
     *
     * @tparam _Type  The type of the real value.
     */
    template<typename _Type, typename = void>
    struct IsTotalZero {
      public:

        using Type = CODI_DD(_Type, double);  ///< See IsTotalZero

        /// Checks if the values are completely zero.
        static CODI_INLINE bool isTotalZero(Type const& v) {
          return Type() == v;
        }
    };

    /// \copydoc codi::TraitsImplementation::Real
    template<typename Type>
    using Real = typename TraitsImplementation<Type>::Real;

    /// \copydoc codi::TraitsImplementation::PassiveReal
    template<typename Type>
    using PassiveReal = typename TraitsImplementation<Type>::PassiveReal;

    /// \copydoc codi::TraitsImplementation::MaxDerivativeOrder
    template<typename Type>
    CODI_INLINE size_t constexpr MaxDerivativeOrder() {
      return TraitsImplementation<Type>::MaxDerivativeOrder;
    }

    /// \copydoc codi::TraitsImplementation::getPassiveValue()
    template<typename Type>
    CODI_INLINE PassiveReal<Type> const& getPassiveValue(Type const& v) {
      return TraitsImplementation<Type>::getPassiveValue(v);
    }

    /// \copydoc codi::IsTotalFinite
    template<typename Type>
    CODI_INLINE bool isTotalFinite(Type const& v) {
      return IsTotalFinite<Type>::isTotalFinite(v);
    }

    /// \copydoc codi::IsTotalZero
    template<typename Type>
    CODI_INLINE bool isTotalZero(Type const& v) {
      return IsTotalZero<Type>::isTotalZero(v);
    }

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
