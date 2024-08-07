/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2024 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../data/position.hpp"
#include "../misc/internalAdjointsInterface.hpp"
#include "../misc/tapeParameters.hpp"
#include "forwardEvaluationTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Allows user defined vectors for the forward and adjoint evaluation, and for clearing adjoints.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * The two additional evaluate methods allow for the evaluation of the tape with a custom adjoint vector, and the
   * additional clearing method allows clearing the custom adjoint vector according to the recorded tape.
   *
   * The adjoint vector type (template parameter AdjointVector in the member functions) must be a accessible with
   * operator[]. Suitable choices are pointers, e.g., Adjoint*, or classes with overloaded operator[] like
   * std::vector<Adjoint>.
   *
   * codi::AdjointVectorTraits::GradientImplementation must be specialized for AdjointVector.
   * The type of the vector entries deduced from these traits (gradient type) must support the following operators:
   *  - operator =
   *  - operator *(Tape::Real, Adjoint) (Scalar multiplication from the left)
   *  - operator +=
   * The gradient type must also specialize #codi::GradientTraits::TraitsImplementation.
   *
   * Here is an example for an evaluation with a custom adjoint vector
   * (documentation/examples/customAdjointVectorEvaluationTapeInterface.cpp):
   * \snippet examples/customAdjointVectorEvaluationTapeInterface.cpp Custom vector
   *
   * @tparam T_Position  Global tape position, usually chosen as Tape::Position.
   */
  template<typename T_Position>
  struct CustomAdjointVectorEvaluationTapeInterface : public virtual ForwardEvaluationTapeInterface<T_Position> {
    public:

      using Position = CODI_DD(T_Position, EmptyPosition);  ///< See CustomAdjointVectorEvaluationTapeInterface.

      /*******************************************************************************/
      /// @name Interface definition

      /**
       * \copybrief codi::PositionalEvaluationTapeInterface::evaluate
       *
       * Tape evaluation with a custom adjoint vector.
       *
       * @tparam AdjointVector  See CustomAdjointVectorEvaluationTapeInterface documentation.
       */
      template<typename AdjointVector>
      void evaluate(Position const& start, Position const& end, AdjointVector&& data);

      // clang-format off
      /**
       * \copybrief codi::ForwardEvaluationTapeInterface::evaluateForward(T_Position const&, T_Position const&, AdjointsManagement)
       *
       * Tape evaluation with a custom adjoint vector.
       *
       * @tparam AdjointVector  See CustomAdjointVectorEvaluationTapeInterface documentation.
       */
      // clang-format on
      template<typename AdjointVector>
      void evaluateForward(Position const& start, Position const& end, AdjointVector&& data);

      /**
       * \copybrief codi::ReverseTapeInterface::clearAdjoints
       *
       * Clear custom adjoint vector according to a tape recording.
       *
       * @tparam AdjointVector  See CustomAdjointVectorEvaluationTapeInterface documentation.
       */
      template<typename AdjointVector>
      void clearCustomAdjoints(Position const& start, Position const& end, AdjointVector&& data);

      /**
       * @brief Obtain a representation of the tape's internal adjoint vector that can be used as custom adjoints.
       *
       * To avoid that functionality has to be implemented both for custom, external and internal adjoints, this method
       * provides access to the internal adjoints so that they can be used as if they were custom adjoints.
       *
       * Warning: If you use this method, proceed with care. If internal adjoints are modified due to side effect of
       * other methods, the object returned here might become invalid, or, conversely, modifications of the returned
       * object other than reading/writing adjoints might interfere with the tape's management of internal adjoints.
       *
       * @tparam InternalAdjoints  Placeholder for the implementation-dependent return type.
       */
      CODI_DD(CODI_T(template<typename InternalAdjoints>), )
      CODI_DD(InternalAdjoints, CODI_T(InternalAdjointsInterface<double, int, CODI_DEFAULT_TAPE>))
      getInternalAdjoints();
  };
}
