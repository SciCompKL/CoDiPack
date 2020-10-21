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

  template<typename _Real, typename _Gradient>
  struct ForwardEvaluation : public InternalStatementRecordingInterface<_Gradient>,
                             public GradientAccessTapeInterface<_Gradient, _Gradient> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);
      using Gradient = CODI_DECLARE_DEFAULT(_Gradient, double);

      using PassiveReal = PassiveRealType<Real>;
      using Identifier = Gradient;

      /*******************************************************************************
       * Section: Implementation of InternalExpressionTapeInterface
       *
       * Description: TODO
       *
       */

      static bool constexpr AllowJacobianOptimization = true;

      template<typename Real>
      void initIdentifier(Real& value, Identifier& identifier) {
        CODI_UNUSED(value);
        identifier = Identifier();
      }

      template<typename Real>
      void destroyIdentifier(Real& value, Identifier& identifier) {
        CODI_UNUSED(value, identifier);
      }

      struct LocalReverseLogic : public JacobianComputationLogic<Real, LocalReverseLogic> {
          template<typename Node>
          CODI_INLINE void handleJacobianOnActive(Node const& node, Real jacobian, Gradient& lhsGradient) {
            CODI_ENABLE_CHECK(Config::IgnoreInvalidJacobies, isTotalFinite(jacobian)) {
              lhsGradient += node.gradient() * jacobian;
            }
          }
      };

      template<typename Lhs, typename Rhs>
      void store(LhsExpressionInterface<Real, Gradient, ForwardEvaluation, Lhs>& lhs,
                 ExpressionInterface<Real, Rhs> const& rhs) {

        LocalReverseLogic reversal;

        Gradient newGradient = Gradient();
        reversal.eval(rhs.cast(), 1.0, newGradient);

        lhs.cast().value() = rhs.cast().getValue();
        lhs.cast().gradient() = newGradient;
      }

      template<typename Lhs, typename Rhs>
      void store(LhsExpressionInterface<Real, Gradient, ForwardEvaluation, Lhs>& lhs,
                 LhsExpressionInterface<Real, Gradient, ForwardEvaluation, Rhs> const& rhs) {

        lhs.cast().value() = rhs.cast().getValue();
        lhs.cast().gradient() = rhs.cast().getGradient();
      }

      template<typename Lhs>
      void store(LhsExpressionInterface<Real, Gradient, ForwardEvaluation, Lhs>& lhs, PassiveReal const& rhs) {
        lhs.cast().value() = rhs;
        lhs.cast().gradient() = Gradient();
      }

      /*******************************************************************************
       * Section: Implementation of GradientAccessTapeInterface
       *
       * Description: TODO
       *
       */

      void setGradient(Identifier& identifier, Gradient const& gradient) {
        identifier = gradient;
      }

      Gradient const& getGradient(Identifier const& identifier) const {
        return identifier;
      }

      Gradient& gradient(Identifier& identifier) {
        return identifier;
      }

      Gradient const& gradient(Identifier const& identifier) const {
        return identifier;
      }

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

  template<typename _Type>
  struct IsTotalFinite<_Type, enableIfForwardTape<typename _Type::Tape>> {
    public:

      using Type = CODI_DECLARE_DEFAULT(
                      _Type,
                      TEMPLATE(LhsExpressionInterface<double, double, InternalExpressionTapeInterface<ANY>, _Type>)
                    );

      static CODI_INLINE bool isTotalFinite(Type const& v) {
        using std::isfinite;
        return isfinite(v.getValue()) && isfinite(v.getGradient());
      }
  };

  template<typename _Type>
  struct IsTotalZero<_Type, enableIfForwardTape<typename _Type::Tape>> {
    public:

      using Type = CODI_DECLARE_DEFAULT(
                      _Type,
                      TEMPLATE(LhsExpressionInterface<double, double, InternalExpressionTapeInterface<ANY>, _Type>)
                    );
      using Real = typename RealTraits<Type>::Real;

      static CODI_INLINE bool isTotalZero(Type const& v) {
        return Real() == v.getValue() && typename Type::Gradient() == v.getGradient();
      }
  };
}

