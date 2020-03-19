#pragma once

#include <cstddef>

#include "../../aux/macros.h"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Identifier>
  struct VectorAccessInterface {

      using Real = DECLARE_DEFAULT(_Real, double);
      using Identifier = DECLARE_DEFAULT(_Identifier, int);

      virtual ~VectorAccessInterface() {}

      virtual size_t getVectorSize() const = 0;

      virtual void resetAdjoint(Identifier const& index, size_t dim) = 0;
      virtual void resetAdjointVec(Identifier const& index) = 0;

      virtual Real getAdjoint(Identifier const& index, size_t dim) = 0;
      virtual void getAdjointVec(Identifier const& index, Real* const vec) = 0;

      virtual void updateAdjoint(Identifier const& index, size_t dim, Real const& adjoint) = 0;
      virtual void updateAdjointVec(Identifier const& index, Real const* const vec) = 0;

      virtual void setLhsAdjoint(Identifier const& index) = 0;
      virtual void updateAdjointWithLhs(Identifier const& index, Real const& jacobi) = 0;

      virtual void setLhsTangent(Identifier const& index) = 0;
      virtual void updateTangentWithLhs(Identifier const& index, Real const& jacobi) = 0;

      virtual bool isLhsZero() = 0;

      virtual void setPrimal(Identifier const& index, Real const& primal) = 0;
      virtual Real getPrimal(Identifier const& index) = 0;

      virtual bool hasPrimals() = 0;
  };
}

