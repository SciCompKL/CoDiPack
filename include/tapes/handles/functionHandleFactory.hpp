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

#include "../../configure.h"
#include "../../typeTraits.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  template<typename Real, typename IndexType, typename GradientValue=Real>
  struct FunctionHandleFactory {

    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    typedef void (*Handle)(const GradientValue& adj, const StatementInt& passiveActives, size_t& indexPos, IndexType* &indices, size_t& constantPos, typename TypeTraits<Real>::PassiveReal* &constants, Real* primalVector, GradientValue* adjoints);

    template<typename Expr, typename Tape>
    static CODI_INLINE Handle createHandle() {

      return &curryEvaluateHandle<Expr, Tape>;
    }

    template<typename Expr, typename Tape>
    static CODI_INLINE void curryEvaluateHandle(const GradientValue& adj, const StatementInt& passiveActives, size_t& indexPos, IndexType* &indices, size_t& constantPos, PassiveReal* &constants, Real* primalVector, GradientValue* adjoints) {
      Tape::evaluateHandle(Expr::template evalAdjoint<IndexType, 0, 0>, ExpressionTraits<Expr>::maxActiveVariables, ExpressionTraits<Expr>::maxConstantVariables, adj, passiveActives, indexPos, indices, constantPos, constants, primalVector, adjoints);
    }

    template<typename Tape>
    static CODI_INLINE void callHandle(Handle handle, const GradientValue& adj, const StatementInt& passiveActives, size_t& indexPos, IndexType* &indices, size_t& constantPos, PassiveReal* &constants, Real* primalVector, GradientValue* adjoints) {

      handle(adj, passiveActives, indexPos, indices, constantPos, constants, primalVector, adjoints);
    }

  };
}
