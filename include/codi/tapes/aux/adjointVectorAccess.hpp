#pragma once

#include <cstddef>

#include "../../aux/macros.h"
#include "../../config.h"
#include "../../traits/realTraits.hpp"
#include "vectorAccessInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Identifier, typename _Gradient>
  struct AdjointVectorAccess : public VectorAccessInterface<_Real, _Identifier>  {

      using Real = DECLARE_DEFAULT(_Real, double);
      using Identifier = DECLARE_DEFAULT(_Identifier, int);
      using Gradient = DECLARE_DEFAULT(_Gradient, double);

      Gradient* adjointVector;

      Gradient lhs;


      AdjointVectorAccess(Gradient* adjointVector) :
        adjointVector(adjointVector),
        lhs() {}

      size_t getVectorSize() const {
        return 1;
      }

      void resetAdjoint(Identifier const& index, size_t dim) {
        CODI_UNUSED(dim);

        adjointVector[index] = Gradient();
      }

      void resetAdjointVec(Identifier const& index) {
        adjointVector[index] = Gradient();
      }

      Real getAdjoint(Identifier const& index, size_t dim) {
        CODI_UNUSED(dim);

        return (Real) adjointVector[index];
      }

      void getAdjointVec(Identifier const& index, Real* const vec) {
        *vec = (Real)adjointVector[index];
      }

      void updateAdjoint(Identifier const& index, size_t dim, Real const& adjoint) {
        CODI_UNUSED(dim);

        adjointVector[index] += adjoint;
      }

      void updateAdjointVec(Identifier const& index, Real const* const vec) {
        adjointVector[index] += *vec;
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
        lhs +=  jacobi * adjointVector[index];
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

