#pragma once

#include <cstddef>

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../tools/data/direction.hpp"
#include "../../traits/realTraits.hpp"
#include "vectorAccessInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Identifier, typename _Gradient>
  struct AdjointVectorAccess : public VectorAccessInterface<_Real, _Identifier>  {

      using Real = CODI_DECLARE_DEFAULT(_Real, double);
      using Identifier = CODI_DECLARE_DEFAULT(_Identifier, int);
      using Gradient = CODI_DECLARE_DEFAULT(_Gradient, double);

      Gradient* adjointVector;

      Gradient lhs;

      AdjointVectorAccess(Gradient* adjointVector) :
        adjointVector(adjointVector),
        lhs() {}

      size_t getVectorSize() const {
        return GradientTraits::dim<Gradient>();
      }

      void resetAdjoint(Identifier const& index, size_t dim) {
        GradientTraits::at(adjointVector[index], dim) = typename GradientTraits::Real<Gradient>();
      }

      void resetAdjointVec(Identifier const& index) {
        adjointVector[index] = Gradient();
      }

      Real getAdjoint(Identifier const& index, size_t dim) {
        CODI_UNUSED(dim);

        return (Real) GradientTraits::at(adjointVector[index], dim);
      }

      void getAdjointVec(Identifier const& index, Real* const vec) {
        for (size_t i = 0; i < getVectorSize(); ++i) {
          vec[i] = (Real) GradientTraits::at(adjointVector[index], i);
        }
      }

      void updateAdjoint(Identifier const& index, size_t dim, Real const& adjoint) {
        GradientTraits::at(adjointVector[index], dim) += adjoint;
      }

      void updateAdjointVec(Identifier const& index, Real const* const vec) {
        for (size_t i = 0; i < getVectorSize(); ++i) {
          GradientTraits::at(adjointVector[index], i) += vec[i];
        }
      }

      void setLhsAdjoint(Identifier const& index) {
        lhs = adjointVector[index];
        adjointVector[index] = Gradient();
      }

      void updateAdjointWithLhs(Identifier const& index, Real const& jacobi) {
        adjointVector[index] += jacobi * lhs;
      }

      void setLhsTangent(Identifier const& index) {
        adjointVector[index] = lhs;
        lhs = Gradient();
      }

      bool isLhsZero() {
        return isTotalZero(lhs);
      }

      void updateTangentWithLhs(Identifier const& index, Real const& jacobi) {
        lhs += jacobi * adjointVector[index];
      }

      void setPrimal(Identifier const& index, Real const& primal) {
        CODI_UNUSED(index, primal);
      }

      Real getPrimal(Identifier const& index) {
        CODI_UNUSED(index);

        return Real();
      }

      bool hasPrimals() {
        return false;
      }
  };
}

