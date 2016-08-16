/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#pragma once

#include "typeTraits.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  template<typename AdjointData, typename Real, typename IndexType>
  class ExpressionHandle {
    private:
      typedef typename TypeTraits<Real>::PassiveReal PassiveReal;
    public:
      typedef void (*StatementFuncPointer)(const Real& seed, const IndexType* indices, const PassiveReal* passiveValues, const Real* primalValues, Real* adjointValues);

      const StatementFuncPointer adjointFunc;
      const size_t maxActiveVariables;
      const size_t maxPassiveVariables;

      ExpressionHandle(const StatementFuncPointer adjointFunc, const size_t maxActiveVariables, const size_t maxPassiveVariables) :
        adjointFunc(adjointFunc),
        maxActiveVariables(maxActiveVariables),
        maxPassiveVariables(maxPassiveVariables) {}
  };

  template<typename AdjointData, typename Real, typename IndexType, typename Expr>
  class ExpressionHandleStore {
    private:
      static const ExpressionHandle<AdjointData, Real, IndexType> handle;
    public:
      static const ExpressionHandle<AdjointData, Real, IndexType>* getHandle() {
        return &handle;
      }
  };

  template<typename AdjointData, typename Real, typename IndexType, typename Expr>
  const ExpressionHandle<AdjointData, Real, IndexType> ExpressionHandleStore<AdjointData, Real, IndexType, Expr>::handle(Expr::template evalAdjointOffset<IndexType, 0, 0>, ExpressionTraits<Expr>::maxActiveVariables, ExpressionTraits<Expr>::maxPassiveVariables);
}
