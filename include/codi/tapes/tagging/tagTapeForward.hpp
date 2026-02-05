/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), RPTU University Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, RPTU University Kaiserslautern-Landau)
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
 *  - SciComp, RPTU University Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../misc/macros.hpp"
#include "../../traits/realTraits.hpp"
#include "../interfaces/gradientAccessTapeInterface.hpp"
#include "../interfaces/internalStatementRecordingTapeInterface.hpp"
#include "tagData.hpp"
#include "tagTapeBase.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Tape for tagging variables and find errors in the AD workflow.
   *
   * Mimics a CoDiPack forward evaluation.
   *
   * See TagTapeBase for detailed information and functionality.
   *
   * @tparam T_Real  The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_tag   The type of the tag, usually int.
   */
  template<typename T_Real, typename T_Tag>
  struct TagTapeForward : public InternalStatementRecordingTapeInterface<T_Tag>,
                          public GradientAccessTapeInterface<T_Tag, T_Tag>,
                          public TagTapeBase<T_Real, T_Tag, TagData<T_Tag>, TagTapeForward<T_Real, T_Tag>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See TagTapeForward.
      using Tag = CODI_DD(T_Tag, double);    ///< See TagTapeForward.

      /// Required definition for event system.
      struct TapeTypes {
          /// Required definition for event system.
          struct IndexManager {
              /// Required definition for event system.
              using Index = int;
          };
      };

      using Gradient = TagData<Tag>;  ///< See TapeTypesInterface.
      using Identifier = Gradient;    ///< Same as the gradient type. Tangent data is stored in the active types.
      using ActiveTypeTapeData = TagData<Tag>;  ///< See TapeTypesInterface.
      using Position = EmptyPosition;           ///< See TapeTypesInterface.

      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.

      using Base = TagTapeBase<T_Real, T_Tag, TagData<T_Tag>, TagTapeForward>;  ///< Base class abbreviation.

    private:

      /// Temporary gradient.
      Gradient tempGradient;

    public:

      /// Constructor.
      TagTapeForward() : Base(), tempGradient() {}

      using Base::AllowJacobianOptimization;
      using Base::destroyTapeData;
      using Base::getIdentifier;
      using Base::initTapeData;
      using Base::store;

      /*******************************************************************************/
      /// @name Implementation of GradientAccessTapeInterface
      /// @{

      /// Do nothing.
      void setGradient(Identifier& identifier, Gradient const& gradient) {
        CODI_UNUSED(identifier, gradient);
      }

      /// Do nothing.
      Gradient const& getGradient(Identifier const& identifier) const {
        CODI_UNUSED(identifier);

        return tempGradient;
      }

      /// Do nothing.
      Gradient& gradient(Identifier& identifier) {
        CODI_UNUSED(identifier);

        return tempGradient;
      }

      /// Do nothing.
      Gradient const& gradient(Identifier const& identifier) const {
        CODI_UNUSED(identifier);

        return tempGradient;
      }

      /// Do not allow.
      void setGradient(Identifier const& identifier, Gradient const& gradient) = delete;

      /// Do not allow.
      Gradient& gradient(Identifier const& identifier) = delete;

      /// @}
  };
}
