#pragma once

#include "../aux/macros.hpp"
#include "../config.h"
#include "../tapes/interfaces/fullTapeInterface.hpp"
#include "../traits/realTraits.hpp"
#include "assignmentOperators.hpp"
#include "incrementOperators.hpp"
#include "lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Represents a concrete lvalue in the CoDiPack expression tree.
   *
   * The class uses members for the storing the value and the identifier.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam _Tape  The tape that manages all expressions created with this type.
   */
  template<typename _Tape>
  struct ActiveType : public LhsExpressionInterface<typename _Tape::Real, typename _Tape::Gradient, _Tape, ActiveType<_Tape> >,
                      public AssignmentOperators<_Tape, ActiveType<_Tape>>,
                      public IncrementOperators<_Tape, ActiveType<_Tape>> {
    public:

      /// See ActiveType.
      /// For reverse AD the tape needs to implement the ReverseTapeInterface
      /// For forward AD the 'tape' needs to implement the InternalStatementRecordingInterface and GradientAccessTapeInterface
      using Tape = CODI_DD(_Tape, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));

      using Real = typename Tape::Real;              ///< See LhsExpressionInterface
      using PassiveReal = RealTraits::PassiveReal<Real>;     ///< Basic computation type
      using Identifier = typename Tape::Identifier;  ///< See LhsExpressionInterface
      using Gradient = typename Tape::Gradient;      ///< See LhsExpressionInterface

      using Base = LhsExpressionInterface<Real, Gradient, Tape, ActiveType>;  ///< Base class abbreviation

    private:

      Real primalValue;
      Identifier identifier;

      static Tape globalTape;

    public:

      /// Constructor
      CODI_INLINE ActiveType() : primalValue(), identifier() {
        Base::init();
      }

      /// Constructor
      CODI_INLINE ActiveType(ActiveType<Tape> const& v) : primalValue(), identifier() {
        Base::init();
        this->getGlobalTape().store(*this, v);
      }

      /// Constructor
      CODI_INLINE ActiveType(PassiveReal const& value) : primalValue(value), identifier() {
        Base::init();
      }

      /// Constructor
      template<class Rhs>
      CODI_INLINE ActiveType(ExpressionInterface<Real, Rhs> const& rhs) : primalValue(), identifier() {
        Base::init();
        this->getGlobalTape().store(*this, rhs.cast());
      }

      /// Destructor
      CODI_INLINE ~ActiveType() {
        Base::destroy();
      }

      /// See LhsExpressionInterface::operator =(ExpressionInterface const&)
      CODI_INLINE ActiveType<Tape>& operator=(ActiveType<Tape> const& v) {
        static_cast<LhsExpressionInterface<Real, Gradient, Tape, ActiveType>&>(*this) = v;
        return *this;
      }
      using LhsExpressionInterface<Real, Gradient, Tape, ActiveType>::operator=;

      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      using StoreAs = ActiveType const&; ///< \copydoc codi::ExpressionInterface::StoreAs

      /// @}
      /*******************************************************************************/
      /// @name Implementation of LhsExpressionInterface
      /// @{

      /// \copydoc codi::LhsExpressionInterface::getIdentifier()
      CODI_INLINE Identifier& getIdentifier() {
        return identifier;
      }

      /// \copydoc codi::LhsExpressionInterface::getIdentifier() const
      CODI_INLINE Identifier const& getIdentifier() const {
        return identifier;
      }

      /// \copydoc codi::LhsExpressionInterface::value()
      CODI_INLINE Real& value() {
        return primalValue;
      }

      /// \copydoc codi::LhsExpressionInterface::value() const
      CODI_INLINE Real const& value() const {
        return primalValue;
      }

      /// \copydoc codi::LhsExpressionInterface::getGlobalTape()
      static CODI_INLINE Tape& getGlobalTape() {
        return globalTape;
      }

      /// @}
  };

  template<typename Tape>
  Tape ActiveType<Tape>::globalTape = Tape();
}
