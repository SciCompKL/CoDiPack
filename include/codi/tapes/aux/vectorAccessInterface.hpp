#pragma once

#include <cstddef>

#include "../../aux/macros.hpp"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Unified access to the adjoint vector and primal vector in a tape evaluation.
   *
   * The interface abstracts the vector access such that custom vectors modes ( \ref Example_11_External_function_user_data)
   * can be handled in a generalized way for external functions. All definitions in this interface are based upon the
   * primal evaluation type of the tape. This means that also all vector definitions need to be evaluated with this
   * type.
   *
   * In general this interface allows to evaluate the \ref sec_forwardAD "forward" and \ref sec_reverseAD "reverse" AD
   * equations. All mathematical symbols in this documentation refer to the linked equations.
   *
   * All identifiers in this interface are tape identifiers and can be obtained with #codi::ActiveType::getIdentifier.
   *
   * The interface provides different access types for the user which can be separated into five categories (all function
   * listed in their typical use order):
   *
   *  - Indirect adjoint access:
   *    - setLhsAdjoint(): Define \f$ \bar w \f$ and set it to zero. Internal copy is kept for updateAdjointWithLhs()
   *    - updateAdjointWithLhs(): Use the lhs defined by setLhsAdjoint() and update \f$ \bar u \f$.
   *
   *  - Indirect tangent access:
   *    - updateTangentWithLhs(): Update an internal value with the \f$ \dot u \f$ value.
   *    - setLhsTangent(): Set \f$ \dot w \f$ to the internal value.
   *
   *  - Direct adjoint vector access: The arrays need to have the size of getVectorSize()
   *    - getAdjointVec(): Get the adjoint vector at the specified location.
   *    - resetAdjointVec(): Reset the adjoint vector at the specified location to zero.
   *    - updateAdjointVec(): Update the adjoint vector ad the specified location with the provided vector.
   *
   *  - Direct adjoint component access:
   *    - Same as the 'direct adjoint vector access' but all function just work on one component.
   *    - Same function without the 'Vec' suffix.
   *
   *  - Primal access: (Optional)
   *    - Only available if 'hasPrimals()' is true
   *    - setPrimal(): Set the primal value
   *    - getPrimal(): Get the primal value
   *    - This access is required for primal values tapes, that need to update or revert primal values during the
   *      tape evaluation.
   *
   * @tparam _Real        The computation type of a tape usually defined by ActiveType::Real.
   * @tparam _Identifier  The adjoint/tangent identification of a tape usually defined by ActiveType::Identifier.
   */
  template<typename _Real, typename _Identifier>
  struct VectorAccessInterface {

      using Real = CODI_DD(_Real, double); ///< See VectorAccessInterface
      using Identifier = CODI_DD(_Identifier, int); ///< See VectorAccessInterface

      virtual ~VectorAccessInterface() {} ///< Destructor

      /*******************************************************************************/
      /// @name Misc

      virtual size_t getVectorSize() const = 0;  ///< Vector size in the current tape evaluation
      virtual bool isLhsZero() = 0;  ///< If adjoint set with setLhsAdjoint is zero.

      /*******************************************************************************/
      /// @name Indirect adjoint access

      virtual void setLhsAdjoint(Identifier const& index) = 0; ///< Set \f$ \bar w \f$ and set it to zero. The lhs is copied internally.
      virtual void updateAdjointWithLhs(Identifier const& index, Real const& jacobi) = 0; ///< Perform \f$ \bar u_{\text{index}} \aeq \text{jacobi} * \bar w \f$

      /*******************************************************************************/
      /// @name Indirect tangent access

      virtual void setLhsTangent(Identifier const& index) = 0; ///< Perform \f$ \dot w = \text{internalMem} \f$. Internal memory is reset afterwards.
      virtual void updateTangentWithLhs(Identifier const& index, Real const& jacobi) = 0; ///< Perform \f$ \text{internalMem} \aeq jacobi * \dot u_{\text{index}} \f$.

      /*******************************************************************************/
      /// @name Direct adjoint access

      virtual void resetAdjoint(Identifier const& index, size_t dim) = 0; ///< Set the adjoint component to zero
      virtual void resetAdjointVec(Identifier const& index) = 0; ///< Set the adjoint entry to zero

      virtual Real getAdjoint(Identifier const& index, size_t dim) = 0; ///< Get the adjoint component
      virtual void getAdjointVec(Identifier const& index, Real* const vec) = 0; ///< Get the adjoint entry

      virtual void updateAdjoint(Identifier const& index, size_t dim, Real const& adjoint) = 0; ///< Update the adjoint component
      virtual void updateAdjointVec(Identifier const& index, Real const* const vec) = 0; ///< Update the adjoint entry

      /*******************************************************************************/
      /// @name Primal access

      virtual void setPrimal(Identifier const& index, Real const& primal) = 0; ///< Set the primal value
      virtual Real getPrimal(Identifier const& index) = 0;  ///< Get the primal value

      virtual bool hasPrimals() = 0;  ///< If the tape/vector interface has primal values.
  };
}

