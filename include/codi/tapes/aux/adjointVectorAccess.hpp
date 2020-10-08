#pragma once

#include <cstddef>

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../tools/data/direction.hpp"
#include "../../traits/realTraits.hpp"
#include "vectorAccessInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Implementation of VectorAccessInterface for adjoint vectors.
   *
   * The adjoint vector is used as is, it needs to have the correct size. No bounds checking is performed
   *
   * @tparam _Real        The computation type of a tape usually defined by ActiveType::Real.
   * @tparam _Identifier  The adjoint/tangent identification of a tape usually defined by ActiveType::Identifier.
   * @tparam _Gradient    The gradient type of a tape usually defined by ActiveType::Gradient.
   */
  template<typename _Real, typename _Identifier, typename _Gradient>
  struct AdjointVectorAccess : public VectorAccessInterface<_Real, _Identifier>  {

      using Real = CODI_DECLARE_DEFAULT(_Real, double); ///< See AdjointVectorAccess
      using Identifier = CODI_DECLARE_DEFAULT(_Identifier, int);  ///< See AdjointVectorAccess
      using Gradient = CODI_DECLARE_DEFAULT(_Gradient, double); ///< See AdjointVectorAccess

    private:

      Gradient* adjointVector; ///< Pointer to the gradient vector

      Gradient lhs; ///< Temporary storage for indirect adjoint or tangent updates.

    public:


      /// Size of adjointVector needs to big enough to. No bounds checking is performed.
      AdjointVectorAccess(Gradient* adjointVector) :
        adjointVector(adjointVector),
        lhs() {}

      /*******************************************************************************/
      /// @name Misc

      /// \copydoc VectorAccessInterface::getVectorSize
      size_t getVectorSize() const {
        return GradientTraits::dim<Gradient>();
      }

      /// \copydoc VectorAccessInterface::isLhsZero
      bool isLhsZero() {
        return isTotalZero(lhs);
      }

      /*******************************************************************************/
      /// @name Indirect adjoint access

      /// \copydoc VectorAccessInterface::setLhsAdjoint
      void setLhsAdjoint(Identifier const& index) {
        lhs = adjointVector[index];
        adjointVector[index] = Gradient();
      }

      /// \copydoc VectorAccessInterface::updateAdjointWithLhs
      void updateAdjointWithLhs(Identifier const& index, Real const& jacobi) {
        adjointVector[index] += jacobi * lhs;
      }

      /*******************************************************************************/
      /// @name Indirect tangent access


      /// \copydoc VectorAccessInterface::setLhsTangent
      void setLhsTangent(Identifier const& index) {
        adjointVector[index] = lhs;
        lhs = Gradient();
      }

      /// \copydoc VectorAccessInterface::updateTangentWithLhs
      void updateTangentWithLhs(Identifier const& index, Real const& jacobi) {
        lhs +=  jacobi * adjointVector[index];
      }

      /*******************************************************************************/
      /// @name Direct adjoint access


      /// \copydoc VectorAccessInterface::resetAdjoint
      void resetAdjoint(Identifier const& index, size_t dim) {
        GradientTraits::at(adjointVector[index], dim) = typename GradientTraits::Real<Gradient>();
      }

      /// \copydoc VectorAccessInterface::resetAdjointVec
      void resetAdjointVec(Identifier const& index) {
        adjointVector[index] = Gradient();
      }

      /// \copydoc VectorAccessInterface::getAdjoint
      Real getAdjoint(Identifier const& index, size_t dim) {
        CODI_UNUSED(dim);

        return (Real) GradientTraits::at(adjointVector[index], dim);
      }

      /// \copydoc VectorAccessInterface::getAdjointVec
      void getAdjointVec(Identifier const& index, Real* const vec) {
        for (size_t i = 0; i < getVectorSize(); ++i) {
          vec[i] = (Real) GradientTraits::at(adjointVector[index], i);
        }
      }

      /// \copydoc VectorAccessInterface::updateAdjoint
      void updateAdjoint(Identifier const& index, size_t dim, Real const& adjoint) {
        GradientTraits::at(adjointVector[index], dim) += adjoint;
      }

      /// \copydoc VectorAccessInterface::updateAdjointVec
      void updateAdjointVec(Identifier const& index, Real const* const vec) {
        for (size_t i = 0; i < getVectorSize(); ++i) {
          GradientTraits::at(adjointVector[index], i) += vec[i];
        }
      }

      /*******************************************************************************/
      /// @name Primal access

      /// \copydoc VectorAccessInterface::setPrimal <br><br>
      /// Implementation: Not implemented, empty function
      void setPrimal(Identifier const& index, Real const& primal) {
        CODI_UNUSED(index, primal);
      }

      /// \copydoc VectorAccessInterface::getPrimal <br><br>
      /// Implementation: Not implemented, returns zero
      Real getPrimal(Identifier const& index) {
        CODI_UNUSED(index);

        return Real();
      }

      /// \copydoc VectorAccessInterface::setPrimal <br><br>
      /// Implementation: Always returns false
      bool hasPrimals() {
        return false;
      }
  };
}

