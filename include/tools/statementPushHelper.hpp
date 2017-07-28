/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2017 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */


#pragma once

#include "../configure.h"
#include "../exceptions.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  template<typename CoDiType>
  struct StatementPushHelper {

      typedef typename CoDiType::Real Real;
      typedef typename CoDiType::GradientData GradientData;

      typedef typename CoDiType::TapeType Tape;

      GradientData indexVector[MaxStatementIntSize];
      Real jacobiVector[MaxStatementIntSize];
      size_t vectorPos;

      void startPushStatement() {
        vectorPos = 0;
      }

      void pushArgument(const CoDiType& arg, const Real& jacobi) {
        Tape& tape = CoDiType::getGlobalTape();

        if(MaxStatementIntSize == vectorPos ) {
          CODI_EXCEPTION("Adding more than %zu arguments to a statement.", MaxStatementIntSize);
        }

        ENABLE_CHECK (OptTapeActivity, tape.isActive()) {
          ENABLE_CHECK(OptCheckZeroIndex, 0 != arg.getGradientData()) {
            ENABLE_CHECK(OptIgnoreInvalidJacobies, codi::isfinite(jacobi)) {
              ENABLE_CHECK(OptJacobiIsZero, !isTotalZero(jacobi)) {

                indexVector[vectorPos] = arg.getGradientData();
                jacobiVector[vectorPos] = jacobi;
                vectorPos += 1;
              }
            }
          }
        }
      }

      void endPushStatement(CoDiType& lhs, const Real& primal) {
        Tape& tape = CoDiType::getGlobalTape();

        ENABLE_CHECK (OptTapeActivity, tape.isActive()) {
          if(0 != vectorPos) {
            tape.storeManual(primal, lhs.getGradientData(), vectorPos);

            for(size_t i = 0; i < vectorPos; ++i) {
              tape.pushJacobiManual(jacobiVector[i], 0.0, indexVector[i]);
            }
          }
        }

        lhs.value() = primal;
      }

      template<typename ArgIter, typename JacobiIter>
      void pushStatement(CoDiType& lhs, const Real& primal, const ArgIter startArg, const ArgIter endArg, const JacobiIter startJac) {

        startPushStatement();

        JacobiIter jacPos = startJac;
        ArgIter argPos = startArg;
        while(argPos != endArg) {
          pushArgument(*argPos, *jacPos);

          ++jacPos;
          ++argPos;
        }

        endPushStatement(lhs, primal);

      }

      template<typename ArgVector, typename JacobiVector>
      void pushStatement(CoDiType& lhs, const Real& primal, const ArgVector& arguments, const JacobiVector& jacobies, const size_t size) {

        startPushStatement();

        for(size_t i = 0; i < size; ++i) {
          pushArgument(arguments[i], jacobies[i]);
        }

        endPushStatement(lhs, primal);
      }
  };

  template<typename CoDiType>
  struct ForwardStatementPushHelper {

      typedef typename CoDiType::Real Real;
      typedef typename CoDiType::GradientValue GradientValue;

      GradientValue lhsTangent;

      void startPushStatement() {
        lhsTangent = GradientValue();
      }

      void pushArgument(const CoDiType& arg, const Real& jacobi) {
        ENABLE_CHECK(OptIgnoreInvalidJacobies, codi::isfinite(jacobi)) {
          lhsTangent += jacobi * arg.getGradient();
        }
      }

      void endPushStatement(CoDiType& lhs, const Real& primal) {
        lhs.gradient() = lhsTangent;
        lhs.value() = primal;
      }

      template<typename ArgIter, typename JacobiIter>
      void pushStatement(CoDiType& lhs, const Real& primal, const ArgIter startArg, const ArgIter endArg, const JacobiIter startJac) {

        startPushStatement();

        JacobiIter jacPos = startJac;
        ArgIter argPos = startArg;
        while(argPos != endArg) {
          pushArgument(*argPos, *jacPos);

          ++jacPos;
          ++argPos;
        }

        endPushStatement(lhs, primal);
      }

      template<typename ArgVector, typename JacobiVector>
      void pushStatement(CoDiType& lhs, const Real& primal, const ArgVector& arguments, const JacobiVector& jacobies, const size_t size) {

        startPushStatement();

        for(size_t i = 0; i < size; ++i) {
          pushArgument(arguments[i], jacobies[i]);
        }

        endPushStatement(lhs, primal);
      }
  };
}
