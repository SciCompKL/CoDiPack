#pragma once

#include <cstddef>

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "adjointVectorAccess.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Identifier, typename _Gradient>
  struct PrimalAdjointVectorAccess : public AdjointVectorAccess<_Real, _Identifier, _Gradient>  {

      using Real = CODI_DECLARE_DEFAULT(_Real, double);
      using Identifier = CODI_DECLARE_DEFAULT(_Identifier, int);
      using Gradient = CODI_DECLARE_DEFAULT(_Gradient, double);

      using Base = AdjointVectorAccess<Real, Identifier, Gradient>;

      Real* primalVector;

      PrimalAdjointVectorAccess(Gradient* adjointVector, Real* primalVector) :
        Base(adjointVector),
        primalVector(primalVector) {}

      void setPrimal(Identifier const& index, Real const& primal) {
        primalVector[index] = primal;
      }

      Real getPrimal(Identifier const& index) {
        return primalVector[index];
      }

      bool hasPrimals() {
        return true;
      }
  };
}
