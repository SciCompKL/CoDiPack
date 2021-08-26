#pragma once

#include "../aux/macros.hpp"
#include "../config.h"
#include "../tapes/interfaces/fullTapeInterface.hpp"
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
   * @tparam T_Real  Original primal value of the statement/expression.
   * @tparam T_Gradient  Gradient values computed by the tape implementation.
   * @tparam T_Tape  The tape that manages the lvalues of the expression.
   *                Minimal interface: InternalStatementRecordingTapeInterface, GradientAccessTapeInterface
   * @tparam T_Impl  Class implementing this interface.
   */
  template<typename T_Real, typename T_Gradient, typename T_Tape, typename T_Impl>
  struct LhsExpressionInterface : public ExpressionInterface<T_Real, T_Impl> {
    public:

      using Real = CODI_DD(T_Real, double);        ///< See LhsExpressionInterface.
      using Gradient = CODI_DD(T_Gradient, Real);  ///< See LhsExpressionInterface.
      using Tape = CODI_DD(T_Tape,
                           CODI_T(FullTapeInterface<double, double, int, CODI_ANY>));  ///< See LhsExpressionInterface.
      using Impl = CODI_DD(T_Impl, LhsExpressionInterface);                             ///< See LhsExpressionInterface.

      using Identifier = typename Tape::Identifier;       ///< See GradientAccessTapeInterface.
      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.

      LhsExpressionInterface() = default;  ///< Constructor
      LhsExpressionInterface(LhsExpressionInterface const& other) = default; ///< Constructor

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      Real const& value() const;  ///< Get a constant reference to the lvalue represented by the expression.
      Real& value();              ///< Get a reference to the lvalue represented by the expression.

      Identifier const& getIdentifier() const;  ///< Get a constant reference to the identifier of the tape for this
                                                ///< expression. See also @ref IdentifierManagement
      Identifier& getIdentifier();  ///< Get a reference to the identifier of the tape for this expression. See also
                                    ///< @ref IdentifierManagement

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

      /// Get the gradient of this lvalue from the tape.
      CODI_INLINE Gradient& gradient() {
        return Impl::getGlobalTape().gradient(cast().getIdentifier());
      }

      /// Get the gradient of this lvalue from the tape.
      CODI_INLINE Gradient const& gradient() const {
        return const_cast<Tape const&>(Impl::getGlobalTape()).gradient(cast().getIdentifier());
      }

      /// Get the gradient of this lvalue from the tape.
      CODI_INLINE Gradient getGradient() const {
        return cast().gradient();
      }

      /// Set the gradient of this lvalue in the tape.
      CODI_INLINE void setGradient(Gradient const& g) {
        cast().gradient() = g;
      }

      /// Get the primal value of this lvalue.
      CODI_INLINE Real const& getValue() const {
        return cast().value();
      }

      /// Set the primal value of this lvalue.
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
      /// To be called in constructors of the implementing class.
      CODI_INLINE void init() {
        Impl::getGlobalTape().initIdentifier(cast().value(), cast().getIdentifier());
      }

      /// Helper function to deconstruct the primal value and the identifier by the tape.
      ///
      /// To be called in the destructor of the implementing class.
      CODI_INLINE void destroy() {
        Impl::getGlobalTape().destroyIdentifier(cast().value(), cast().getIdentifier());
      }

      /// @}
  };

  /// Read the primal value from a stream
  template<typename Expr>
  ExpressionTraits::EnableIfLhsExpression<Expr, std::istream>& operator>>(std::istream& stream, Expr& v){
    typename Expr::Real temp;

    stream >> temp;
    v.setValue(temp);

    return stream;
  }


#ifndef DOXYGEN_DISABLE

  /// Specialization of RealTraits::DataExtraction for CoDiPack types.
  template<typename T_Type>
  struct RealTraits::DataExtraction<T_Type, ExpressionTraits::EnableIfLhsExpression<T_Type>> {
    public:
      using Type = CODI_DD(T_Type,
                           CODI_T(LhsExpressionInterface<double, int, CODI_ANY, CODI_ANY>));  ///< See DataExtraction.

      using Real = typename Type::Real;              ///< See DataExtraction::Real.
      using Identifier = typename Type::Identifier;  ///< See DataExtraction::Identifier.

      /// \copydoc DataExtraction::getValue()
      CODI_INLINE static Real getValue(Type const& v) {
        return v.getValue();
      }

      /// \copydoc DataExtraction::getIdentifier()
      CODI_INLINE static Identifier getIdentifier(Type const& v) {
        return v.getIdentifier();
      }

      /// \copydoc DataExtraction::setValue()
      CODI_INLINE static void setValue(Type& v, Real const& value) {
        v.setValue(value);
      }
  };

  /// Specialization of RealTraits::DataRegistration for CoDiPack types.
  template<typename T_Type>
  struct RealTraits::TapeRegistration<T_Type, ExpressionTraits::EnableIfLhsExpression<T_Type>> {
      using Type = CODI_DD(T_Type,
                           CODI_T(LhsExpressionInterface<double, int, CODI_ANY, CODI_ANY>));  ///< See DataRegistration.

      using Real = typename DataExtraction<Type>::Real;  ///< See DataExtraction::Real.

      /// \copydoc DataRegistration::registerInput()
      CODI_INLINE static void registerInput(Type& v) {
        Type::getGlobalTape().registerInput(v);
      }

      /// \copydoc DataRegistration::registerOutput()
      CODI_INLINE static void registerOutput(Type& v) {
        Type::getGlobalTape().registerOutput(v);
      }

      /// \copydoc DataRegistration::registerExternalFunctionOutput()
      CODI_INLINE static Real registerExternalFunctionOutput(Type& v) {
        return Type::getGlobalTape().registerExternalFunctionOutput(v);
      }
  };
#endif
}
