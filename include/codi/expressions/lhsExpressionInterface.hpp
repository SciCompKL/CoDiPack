#pragma once

#include <iostream>

#include "../aux/macros.hpp"
#include "../config.h"
#include "../tapes/interfaces/gradientAccessTapeInterface.hpp"
#include "../tapes/interfaces/internalStatementRecordingTapeInterface.hpp"
#include "../traits/expressionTraits.hpp"
#include "../traits/realTraits.hpp"
#include "expressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Base class for all CoDiPack lvalue expression.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * This interface resembles an lvalue in C++.
   *
   * @tparam _Real  Original primal value of the statement/expression.
   * @tparam _Gradient  Gradient values computed by the tape implementation.
   * @tparam _Tape  The tape that manages the lvalues of the expression.
   * @tparam _Impl  Class implementing this interface.
   */
  template<typename _Real, typename _Gradient, typename _Tape, typename _Impl>
  struct LhsExpressionInterface : public ExpressionInterface<_Real, _Impl> {
    public:

      using Real = CODI_DD(_Real, double);        ///< See LhsExpressionInterface
      using Gradient = CODI_DD(_Gradient, Real);  ///< See LhsExpressionInterface
      using Tape =
          CODI_DD(_Tape, CODI_T(CODI_UNION<InternalStatementRecordingTapeInterface<int>,
                                           GradientAccessTapeInterface<double, int>>));  ///< See LhsExpressionInterface
      using Impl = CODI_DD(_Impl, LhsExpressionInterface);                               ///< See LhsExpressionInterface

      using Identifier = typename Tape::Identifier;       ///< See GradientAccessTapeInterface
      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      Real const& value() const;  ///< Get a constant reference to the lvalue represented by the expression.
      Real& value();              ///< Get a reference to the lvalue represented by the expression.

      Identifier const& getIdentifier() const;  ///< Get a constant reference to the identifier of the tape for this
                                                ///< expression.
      Identifier& getIdentifier();  ///< Get a constant reference to the identifier of the tape for this expression.

      static Tape& getGlobalTape();  ///< Get a reference to the tape which manages this expression.

      /// @}
      /*******************************************************************************/
      /// @name General implementation
      /// @{

      /// Cast to the implementation.
      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);

      }
      using ExpressionInterface<Real, Impl>::cast;

      /// Get the gradient of this lvalue from the tape
      CODI_INLINE Gradient& gradient() {
        return Impl::getGlobalTape().gradient(cast().getIdentifier());
      }

      /// Get the gradient of this lvalue from the tape
      CODI_INLINE Gradient const& gradient() const {
        return const_cast<Tape const&>(Impl::getGlobalTape()).gradient(cast().getIdentifier());
      }

      /// Get the gradient of this lvalue from the tape
      CODI_INLINE Gradient getGradient() const {
        return cast().gradient();
      }

      /// Set the gradient of this lvalue in the tape
      CODI_INLINE void setGradient(Gradient const& g) {
        cast().gradient() = g;
      }

      /// Get the primal value of this lvalue
      CODI_INLINE Real const& getValue() const {
        return cast().value();
      }

      /// Set the primal value of this lvalue
      CODI_INLINE void setValue(Real const& v) {
        cast().value() = v;
      }

      /// Assignment operator for passive values. Calls store on the InternalStatementRecordingTapeInterface.
      CODI_INLINE Impl& operator=(PassiveReal const& rhs) {
        Impl::getGlobalTape().store(cast(), rhs);
        return cast();
      }

      /// Assignment operator for expressions. Calls store on the InternalStatementRecordingTapeInterface.
      template<typename Rhs>
      CODI_INLINE Impl& operator=(ExpressionInterface<Real, Rhs> const& rhs) {
        Impl::getGlobalTape().store(cast(), rhs.cast());
        return cast();
      }

      /// Assignment operator for lhs expressions. Calls store on the InternalStatementRecordingTapeInterface.
      CODI_INLINE Impl& operator=(LhsExpressionInterface const& rhs) {
        Impl::getGlobalTape().store(cast(), rhs);
        return cast();
      }

      /// @}
      /*******************************************************************************/
      /// @name Implementation of NodeInterface
      /// @{

      static bool constexpr EndPoint = true;  ///< \copydoc codi::NodeInterface::EndPoint

      /// \copydoc codi::NodeInterface::forEachLink
      template<typename Logic, typename... Args>
      CODI_INLINE void forEachLink(TraversalLogic<Logic>& logic, Args&&... args) const {
        CODI_UNUSED(logic, args...);
      }

      /// \copydoc codi::NodeInterface::forEachLinkConstExpr
      template<typename Logic, typename... Args>
      CODI_INLINE static typename Logic::ResultType constexpr forEachLinkConstExpr(Args&&... CODI_UNUSED_ARG(args)) {
        return Logic::NeutralElement;
      }

    protected:

      /// Helper function to initialize the primal value and the identifier by the tape.
      ///
      /// Only needs to be called if the members are freshly created.
      CODI_INLINE void init() {
        Impl::getGlobalTape().initIdentifier(cast().value(), cast().getIdentifier());
      }

      /// Helper function to deconstruct the primal value and the identifier by the tape.
      ///
      /// Only needs to be called if the members will get destroyed.
      CODI_INLINE void destroy() {
        Impl::getGlobalTape().destroyIdentifier(cast().value(), cast().getIdentifier());
      }

      /// @}
  };

  /// Write the primal value to the stream.
  template<typename Real, typename Gradient, typename Tape, typename Impl>
  std::ostream& operator<<(std::ostream& out, LhsExpressionInterface<Real, Gradient, Tape, Impl> const& v) {
    return out << v.cast().value();
  }

#ifndef DOXYGEN_DISABLE
  template<typename _Type>
  struct RealTraits::TraitsImplementation<_Type, ExpressionTraits::EnableIfLhsExpression<_Type>> {
    public:

      using Type = CODI_DD(
          _Type, CODI_T(LhsExpressionInterface<double, double, InternalStatementRecordingTapeInterface<CODI_ANY>, _Type>));
      using Real = typename Type::Real;

      using PassiveReal = RealTraits::PassiveReal<Real>;

      static int constexpr MaxDerivativeOrder = 1 + RealTraits::MaxDerivativeOrder<Real>();

      static CODI_INLINE PassiveReal const& getPassiveValue(Type const& v) {
        return RealTraits::getPassiveValue(v.getValue());
      }
  };
#endif
}