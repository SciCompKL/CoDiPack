#pragma once

#include "../aux/macros.hpp"
#include "../config.h"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../expressions/logic/helpers/jacobianComputationLogic.hpp"
#include "../traits/expressionTraits.hpp"
#include "../traits/realTraits.hpp"
#include "../traits/tapeTraits.hpp"
#include "interfaces/gradientAccessTapeInterface.hpp"
#include "interfaces/internalStatementRecordingInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Implementation of a forward AD mode through the internal expression interfaces.
   *
   * For an explanation of the forward AD mode please see the Section \ref sec_forwardAD "forward AD".
   *
   * The store method implementation performs a reverse AD sweep on the expression itself. The result is then added to
   * the tangent data of the left hand side type.
   *
   * The identifier data in the LhsExpressionInterface implementations is used by this class to store the tangent data
   * for each value.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * @tparam _Real        The computation type of a tape usually defined by ActiveType::Real.
   * @tparam _Gradient    The gradient type of a tape usually defined by ActiveType::Gradient.

   */
  template<typename _Real, typename _Gradient>
  struct ForwardEvaluation : public InternalStatementRecordingInterface<_Gradient>,
                             public GradientAccessTapeInterface<_Gradient, _Gradient> {
    public:

      using Real = CODI_DD(_Real, double);  ///< See ForwardEvaluation
      using Gradient = CODI_DD(_Gradient, double); ///< See ForwardEvaluation

      using PassiveReal = RealTraits::PassiveReal<Real>; ///< Basic computation type
      using Identifier = Gradient;  ///< Same as the gradient type. Tangent data is stored in the active types.

      /*******************************************************************************/
      /// @name Implementation of InternalStatementRecordingInterface
      /// @{

      static bool constexpr AllowJacobianOptimization = true;  ///< See InternalStatementRecordingInterface

      /// \copydoc codi::InternalStatementRecordingInterface::initIdentifier()
      template<typename Real>
      void initIdentifier(Real& value, Identifier& identifier) {
        CODI_UNUSED(value);
        identifier = Identifier();
      }

      /// \copydoc codi::InternalStatementRecordingInterface::destroyIdentifier()
      template<typename Real>
      void destroyIdentifier(Real& value, Identifier& identifier) {
        CODI_UNUSED(value, identifier);
      }

      /// @}

    private:

      struct LocalReverseLogic : public JacobianComputationLogic<Real, LocalReverseLogic> {
          template<typename Node>
          CODI_INLINE void handleJacobianOnActive(Node const& node, Real jacobian, Gradient& lhsGradient) {
            CODI_ENABLE_CHECK(Config::IgnoreInvalidJacobies, RealTraits::isTotalFinite(jacobian)) {
              lhsGradient += node.gradient() * jacobian;
            }
          }
      };

    public:

      /// @{

      /// \copydoc codi::InternalStatementRecordingInterface::store()
      template<typename Lhs, typename Rhs>
      void store(LhsExpressionInterface<Real, Gradient, ForwardEvaluation, Lhs>& lhs,
                 ExpressionInterface<Real, Rhs> const& rhs) {

        LocalReverseLogic reversal;

        Gradient newGradient = Gradient();
        reversal.eval(rhs.cast(), 1.0, newGradient);

        lhs.cast().value() = rhs.cast().getValue();
        lhs.cast().gradient() = newGradient;
      }

      /// \copydoc codi::InternalStatementRecordingInterface::store() <br>
      /// Optimization for copy statements.
      template<typename Lhs, typename Rhs>
      void store(LhsExpressionInterface<Real, Gradient, ForwardEvaluation, Lhs>& lhs,
                 LhsExpressionInterface<Real, Gradient, ForwardEvaluation, Rhs> const& rhs) {

        lhs.cast().value() = rhs.cast().getValue();
        lhs.cast().gradient() = rhs.cast().getGradient();
      }

      /// \copydoc codi::InternalStatementRecordingInterface::store() <br>
      /// Specialization for passive assignments.
      template<typename Lhs>
      void store(LhsExpressionInterface<Real, Gradient, ForwardEvaluation, Lhs>& lhs, PassiveReal const& rhs) {
        lhs.cast().value() = rhs;
        lhs.cast().gradient() = Gradient();
      }

      /// @}
      /*******************************************************************************/
      /// @name Implementation of GradientAccessTapeInterface
      /// @{

      /// \copydoc codi::GradientAccessTapeInterface::setGradient()
      void setGradient(Identifier& identifier, Gradient const& gradient) {
        identifier = gradient;
      }

      /// \copydoc codi::GradientAccessTapeInterface::getGradient()
      Gradient const& getGradient(Identifier const& identifier) const {
        return identifier;
      }

      /// \copydoc codi::GradientAccessTapeInterface::gradient(Identifier const&)
      Gradient& gradient(Identifier& identifier) {
        return identifier;
      }

      /// \copydoc codi::GradientAccessTapeInterface::gradient(Identifier const&) const
      Gradient const& gradient(Identifier const& identifier) const {
        return identifier;
      }

      /// @}

    private:

      void setGradient(Identifier const& identifier, Gradient const& gradient) {
        CODI_UNUSED(identifier, gradient);
      }

      Gradient& gradient(Identifier const& identifier) {
        CODI_UNUSED(identifier);
        static Gradient temp = Gradient();
        return temp;
      }
  };


  /// \copydoc codi::RealTraits::IsTotalFinite <br>
  /// Value and gradient are tested if they are finite.
  template<typename _Type>
  struct RealTraits::IsTotalFinite<_Type, TapeTraits::EnableIfForwardTape<typename _Type::Tape>> {
    public:

      using Type = CODI_DD(
                      _Type,
                      TEMPLATE(LhsExpressionInterface<double, double, InternalExpressionTapeInterface<ANY>, _Type>)
                    ); ///< See RealTraits::IsTotalFinite

      /// \copydoc codi::RealTraits::IsTotalFinite::isTotalFinite()
      static CODI_INLINE bool isTotalFinite(Type const& v) {
        using std::isfinite;
        return RealTraits::isTotalFinite(v.getValue()) && RealTraits::isTotalFinite(v.getGradient());
      }
  };

  /// \copydoc codi::RealTraits::IsTotalZero <br>
  /// Value and gradient are tested if they are zero.
  template<typename _Type>
  struct RealTraits::IsTotalZero<_Type, TapeTraits::EnableIfForwardTape<typename _Type::Tape>> {
    public:

      using Type = CODI_DD(
                      _Type,
                      TEMPLATE(LhsExpressionInterface<double, double, InternalExpressionTapeInterface<ANY>, _Type>)
                    ); ///< See RealTraits::IsTotalZero
      using Real = typename Type::Real; ///< See codi::LhsExpressionInterface::Real

      /// \copydoc codi::RealTraits::IsTotalFinite::isTotalZero()
      static CODI_INLINE bool isTotalZero(Type const& v) {
        return Real() == v.getValue() && typename Type::Gradient() == v.getGradient();
      }
  };
}

