#pragma once

#include <iostream>

#include "../aux/macros.hpp"
#include "../config.h"
#include "../traits/expressionTraits.hpp"
#include "../traits/realTraits.hpp"
#include "logic/nodeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Base class for all CoDiPack expressions.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * This interface resembles an rvalue in C++.
   *
   * @tparam T_Real  Original primal value of the statement/expression.
   * @tparam T_Impl  Class implementing this interface.
   */
  template<typename T_Real, typename T_Impl>
  struct ExpressionInterface : public NodeInterface<T_Impl> {
    public:

      using Real = CODI_DD(T_Real, double);               ///< See ExpressionInterface.
      using Impl = CODI_DD(T_Impl, ExpressionInterface);  ///< See ExpressionInterface.

      using ActiveResult = CODI_UNDEFINED;  ///< Type into which the expression can be converted. Usually also the type
                                            ///< from which is constructed.

      /// Cast to the implementation.
      CODI_INLINE Impl const& cast() const {
        return static_cast<Impl const&>(*this);
      }

#if CODI_ImplicitConversion
      /// Implicit cast for CoDiPack expressions.
      CODI_INLINE operator const Real() const {
        Warning::implicitCast<Config::ImplicitConversionWarning>();

        return cast().getValue();
      }
#endif

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      using StoreAs = ExpressionInterface;  ///< Defines how this expression is stored in an expression tree.

      /// Compute the primal value that is usually evaluated by the statement/expression.
      CODI_INLINE Real const getValue() const;

      /// Get the Jacobian with respect to the given argument.
      ///
      /// This is just the local Jacobian and not the one for the whole expression tree.
      template<size_t argNumber>
      CODI_INLINE Real getJacobian() const;

      /// @}

    private:
      ExpressionInterface& operator=(ExpressionInterface const&) = delete;
  };

#ifndef DOXYGEN_DISABLE
  template<typename T_Type>
  struct RealTraits::TraitsImplementation<T_Type, ExpressionTraits::EnableIfExpression<T_Type>> {
    public:

      using Type = CODI_DD(T_Type, CODI_T(ExpressionInterface<double, T_Type>));
      using Real = typename Type::Real;

      using PassiveReal = RealTraits::PassiveReal<Real>;

      static int constexpr MaxDerivativeOrder = 1 + RealTraits::MaxDerivativeOrder<Real>();

      static CODI_INLINE PassiveReal const& getPassiveValue(Type const& v) {
        return RealTraits::getPassiveValue(v.getValue());
      }
  };
#endif

  /// Write the primal value to the stream.
  template<typename Expr>
  ExpressionTraits::EnableIfExpression<Expr, std::ostream>& operator<<(std::ostream& out, Expr const& v) {
    out << v.getValue();

    return out;
  }
}
