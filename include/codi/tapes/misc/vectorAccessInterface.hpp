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

#include <cstddef>

#include "../../config.h"
#include "../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Unified access to the adjoint vector and primal vector in a tape evaluation.
   *
   * The interface abstracts the vector access such that custom vectors modes ( \ref
   * Example_11_External_function_user_data) can be handled in a generalized way for external functions. All definitions
   * in this interface are based on the primal evaluation type of the tape. This means that also all vector
   * definitions need to be evaluated with this type.
   *
   * In general, this interface allows to evaluate the \ref sec_forwardAD "forward" and \ref sec_reverseAD "reverse" AD
   * equations. All mathematical symbols in this documentation refer to the linked equations.
   *
   * All identifiers in this interface are tape identifiers and can be obtained with #codi::ActiveType::getIdentifier.
   *
   * The interface provides different access types for the user which can be separated into five categories (all
   * functions listed in their typical order of use):
   *
   *  - Indirect adjoint access:
   *    - setLhsAdjoint(): Identify the lhs variable \f$ w \f$ of the forward statement. Create an internal copy of
   *      \f$ \bar w \f$ for updateAdjointWithLhs(), and set \f$ \bar w \f$ to zero.
   *    - updateAdjointWithLhs(): Use the lhs defined by setLhsAdjoint() and perform the corrsponding update on
   *      \f$ \bar u \f$.
   *
   *  - Indirect tangent access:
   *    - updateTangentWithLhs(): Update an internal value with the \f$ \dot u \f$ value.
   *    - setLhsTangent(): Set \f$ \dot w \f$ to the internal value.
   *
   *  - Direct adjoint vector access: The arrays need to have the size of getVectorSize()
   *    - getAdjointVec(): Get the adjoint vector at the specified location.
   *    - resetAdjointVec(): Reset the adjoint vector at the specified location to zero.
   *    - updateAdjointVec(): Update the adjoint vector at the specified location with the provided components.
   *
   *  - Direct adjoint component access:
   *    - Same as the 'direct adjoint vector access' but all functions just work on one component.
   *    - Same function without the 'Vec' suffix.
   *
   *  - Primal access: (Optional)
   *    - Only available if 'hasPrimals()' is true
   *    - setPrimal(): Set the primal value
   *    - getPrimal(): Get the primal value
   *    - This access is required for primal values tapes, which need to update or revert primal values during the
   *      tape evaluation.
   *
   * @tparam T_Real        The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_Identifier  The adjoint/tangent identification of a tape, usually chosen as ActiveType::Identifier.
   */
  template<typename T_Real, typename T_Identifier>
  struct VectorAccessInterface {
    public:

      using Real = CODI_DD(T_Real, double);           ///< See VectorAccessInterface.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See VectorAccessInterface.

      virtual ~VectorAccessInterface() {}  ///< Destructor

      /*******************************************************************************/
      /// @name Misc

      virtual size_t getVectorSize() const = 0;          ///< Vector size in the current tape evaluation.
      virtual bool isLhsZero() = 0;                      ///< True if the adjoint set with setLhsAdjoint is zero.
      virtual VectorAccessInterface* clone() const = 0;  ///< Obtain a heap-allocated copy of the vector access inter-
                                                         ///< face. The user is responsible for deleting the pointer.

      /*******************************************************************************/
      /// @name Indirect adjoint access

      virtual void setLhsAdjoint(Identifier const& index) = 0;  ///< Identify the lhs variable \f$ w \f$. Create an
                                                                ///< internal copy of \f$ \bar w \f$ and set
                                                                ///< \f$ \bar w \f$ to zero.
      virtual void updateAdjointWithLhs(Identifier const& index,
                                        Real const& jacobian) = 0;  ///< Perform \f$ \bar u_{\text{index}} \aeq
                                                                    /// \text{jacobian} * \bar w \f$

      /*******************************************************************************/
      /// @name Indirect tangent access

      virtual void setLhsTangent(Identifier const& index) = 0;  ///< Perform \f$ \dot w = \text{internalMem} \f$.
                                                                ///< Internal memory is reset afterwards.
      virtual void updateTangentWithLhs(Identifier const& index,
                                        Real const& jacobian) = 0;  ///< Perform \f$ \text{internalMem} \aeq jacobian *
                                                                    ///< \dot
                                                                    /// u_{\text{index}} \f$.

      /*******************************************************************************/
      /// @name Direct adjoint access

      virtual void resetAdjoint(Identifier const& index, size_t dim) = 0;  ///< Set the adjoint component to zero.
      virtual void resetAdjointVec(Identifier const& index) = 0;           ///< Set the adjoint entry to zero.

      virtual Real getAdjoint(Identifier const& index, size_t dim) = 0;          ///< Get the adjoint component.
      virtual void getAdjointVec(Identifier const& index, Real* const vec) = 0;  ///< Get the adjoint entry.
      virtual Real const* getAdjointVec(Identifier const& index) = 0;            ///< Get the adjoint entry.

      virtual void updateAdjoint(Identifier const& index, size_t dim,
                                 Real const& adjoint) = 0;  ///< Update the adjoint component.
      virtual void updateAdjointVec(Identifier const& index, Real const* const vec) = 0;  ///< Update the adjoint entry.

      /*******************************************************************************/
      /// @name Primal access

      virtual void setPrimal(Identifier const& index, Real const& primal) = 0;  ///< Set the primal value.
      virtual Real getPrimal(Identifier const& index) = 0;                      ///< Get the primal value.

      virtual bool hasPrimals() = 0;  ///< True if the tape/vector interface has primal values.
  };
}
