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

#include "../config.h"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../expressions/logic/helpers/jacobianComputationLogic.hpp"
#include "../misc/macros.hpp"
#include "../traits/expressionTraits.hpp"
#include "../traits/realTraits.hpp"
#include "../traits/tapeTraits.hpp"
#include "interfaces/gradientAccessTapeInterface.hpp"
#include "interfaces/internalStatementRecordingTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Implementation of a tape-free forward AD mode through the internal expression interfaces.
   *
   * For an explanation of the forward AD mode please see the Section \ref sec_forwardAD "forward AD".
   *
   * The store method implementation performs a reverse AD sweep on the expression itself. The result is then added to
   * the tangent data of the left hand side type.
   *
   * The identifier data in the LhsExpressionInterface implementations is used by this class to store the tangent data
   * for each value.
   *
   * This way, a tape-free foward mode is implemented in a manner that is consistent with CoDiPack's taping interface,
   * even if no tape is actually recorded.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * @tparam T_Real        The computation type of a tape usually defined by ActiveType::Real.
   * @tparam T_Gradient    The gradient type of a tape usually defined by ActiveType::Gradient.

   */
  template<typename T_Real, typename T_Gradient>
  struct ForwardEvaluation : public InternalStatementRecordingTapeInterface<T_Gradient>,
                             public GradientAccessTapeInterface<T_Gradient, T_Gradient> {
    public:

      using Real = CODI_DD(T_Real, double);          ///< See ForwardEvaluation.
      using Gradient = CODI_DD(T_Gradient, double);  ///< See ForwardEvaluation.

      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.
      using Identifier = Gradient;  ///< Same as the gradient type. Tangent data is stored in the active types.

      /*******************************************************************************/
      /// @name Implementation of InternalStatementRecordingTapeInterface
      /// @{

      static bool constexpr AllowJacobianOptimization = true;  ///< See InternalStatementRecordingTapeInterface

      /// \copydoc codi::InternalStatementRecordingTapeInterface::initIdentifier()
      template<typename Real>
      void initIdentifier(Real& value, Identifier& identifier) {
        CODI_UNUSED(value);
        identifier = Identifier();
      }

      /// \copydoc codi::InternalStatementRecordingTapeInterface::destroyIdentifier()
      template<typename Real>
      void destroyIdentifier(Real& value, Identifier& identifier) {
        CODI_UNUSED(value, identifier);
      }

      /// @}

    private:

      struct LocalReverseLogic : public JacobianComputationLogic<LocalReverseLogic> {
        public:
          template<typename Node>
          CODI_INLINE void handleJacobianOnActive(Node const& node, Real jacobian, Gradient& lhsGradient) {
            if (CODI_ENABLE_CHECK(Config::IgnoreInvalidJacobians, RealTraits::isTotalFinite(jacobian))) {
              lhsGradient += node.gradient() * jacobian;
            }
          }
      };

    public:

      /// @{

      /// \copydoc codi::InternalStatementRecordingTapeInterface::store()
      template<typename Lhs, typename Rhs>
      void store(LhsExpressionInterface<Real, Gradient, ForwardEvaluation, Lhs>& lhs,
                 ExpressionInterface<Real, Rhs> const& rhs) {
        LocalReverseLogic reversal;

        Gradient newGradient = Gradient();
        reversal.eval(rhs.cast(), 1.0, newGradient);

        lhs.cast().value() = rhs.cast().getValue();
        lhs.cast().gradient() = newGradient;
      }

      /// \copydoc codi::InternalStatementRecordingTapeInterface::store() <br>
      /// Optimization for copy statements.
      template<typename Lhs, typename Rhs>
      void store(LhsExpressionInterface<Real, Gradient, ForwardEvaluation, Lhs>& lhs,
                 LhsExpressionInterface<Real, Gradient, ForwardEvaluation, Rhs> const& rhs) {
        lhs.cast().value() = rhs.cast().getValue();
        lhs.cast().gradient() = rhs.cast().getGradient();
      }

      /// \copydoc codi::InternalStatementRecordingTapeInterface::store() <br>
      /// Specialization for passive assignments.
      template<typename Lhs>
      void store(LhsExpressionInterface<Real, Gradient, ForwardEvaluation, Lhs>& lhs, Real const& rhs) {
        lhs.cast().value() = rhs;
        lhs.cast().gradient() = Gradient();
      }

      /// @}
      /*******************************************************************************/
      /// @name Implementation of GradientAccessTapeInterface
      /// @{

      /// \copydoc codi::GradientAccessTapeInterface::setGradient()
      /// <br> Implementation: Automatic adjoints management has no effect. The forward mode does not maintain internal
      /// adjoints.
      void setGradient(Identifier& identifier, Gradient const& gradient,
                       AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        CODI_UNUSED(adjointsManagement);
        identifier = gradient;
      }

      /// \copydoc codi::GradientAccessTapeInterface::getGradient()
      /// <br> Implementation: Automatic adjoints management has no effect. The forward mode does not maintain internal
      /// adjoints.
      Gradient const& getGradient(Identifier const& identifier,
                                  AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) const {
        CODI_UNUSED(adjointsManagement);
        return identifier;
      }

      /// \copydoc codi::GradientAccessTapeInterface::gradient(Identifier const&, AdjointsManagement)
      /// <br> Implementation: Automatic adjoints management has no effect. The forward mode does not maintain internal
      /// adjoints.
      Gradient& gradient(Identifier& identifier,
                         AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        CODI_UNUSED(adjointsManagement);
        return identifier;
      }

      /// \copydoc codi::GradientAccessTapeInterface::gradient(Identifier const&, AdjointsManagement) const
      /// <br> Implementation: Automatic adjoints management has no effect. The forward mode does not maintain internal
      /// adjoints.
      Gradient const& gradient(Identifier const& identifier,
                               AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) const {
        CODI_UNUSED(adjointsManagement);
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
  /// Tests whether both value and gradient are finite.
  template<typename T_Type>
  struct RealTraits::IsTotalFinite<T_Type, TapeTraits::EnableIfForwardTape<typename T_Type::Tape>> {
    public:

      using Type = CODI_DD(
          T_Type, CODI_T(LhsExpressionInterface<double, double, InternalStatementRecordingTapeInterface<CODI_ANY>,
                                                T_Type>));  ///< See RealTraits::IsTotalFinite.

      /// \copydoc codi::RealTraits::IsTotalFinite::isTotalFinite()
      static CODI_INLINE bool isTotalFinite(Type const& v) {
        using std::isfinite;
        return RealTraits::isTotalFinite(v.getValue()) && RealTraits::isTotalFinite(v.getGradient());
      }
  };

  /// \copydoc codi::RealTraits::IsTotalZero <br>
  /// Tests whether both value and gradient are zero.
  template<typename T_Type>
  struct RealTraits::IsTotalZero<T_Type, TapeTraits::EnableIfForwardTape<typename T_Type::Tape>> {
    public:

      using Type = CODI_DD(
          T_Type, CODI_T(LhsExpressionInterface<double, double, InternalStatementRecordingTapeInterface<CODI_ANY>,
                                                T_Type>));  ///< See RealTraits::IsTotalZero.
      using Real = typename Type::Real;                     ///< See
                                                            ///< codi::LhsExpressionInterface::Real.

      /// \copydoc codi::RealTraits::IsTotalFinite::isTotalZero()
      static CODI_INLINE bool isTotalZero(Type const& v) {
        return Real() == v.getValue() && typename Type::Gradient() == v.getGradient();
      }
  };
}
