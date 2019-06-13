/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
#include "../gradientTraits.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  template<typename CoDiType, bool ZeroChecks = true>
  struct Algorithms {

    private:
      using Tape = typename CoDiType::TapeType;
      using Position = typename Tape::Position;
      using Real = typename CoDiType::Real;
      using GradientData = typename CoDiType::GradientData;
      using GradientValue = typename CoDiType::GradientValue;

      using GT = GradientValueTraits<GradientValue>;

    public:
      enum class EvaluationType {
          Forward,
          Reverse
      };

      static CODI_INLINE EvaluationType getEvaluationChoice(size_t const inputs, size_t const outputs) {
        if(inputs < outputs) {
          return EvaluationType::Forward;
        } else {
          return EvaluationType::Reverse;
        }
      }

      template<typename Jac, bool usePreAccEval = true>
      static CODI_INLINE void computeJacobian(
          Tape& tape, Position const& start, Position const& end,
          GradientData const * input, size_t const inputSize,
          GradientData const * output, size_t const outputSize,
          Jac& jac)
      {
        constexpr size_t gradDim = GT::getVectorSize();

        EvaluationType evalType = getEvaluationChoice(inputSize, outputSize);
        if(EvaluationType::Forward == evalType) {

          for(size_t j = 0; j < inputSize; j += gradDim) {
            seedGradient(tape, j, input, inputSize);

            if(usePreAccEval) {
              tape.evaluateForwardPreacc(start, end);
            } else {
              tape.evaluateForward(start, end);
            }

            for(size_t i = 0; i < outputSize; i += 1) {
              for(size_t curDim = 0; curDim < gradDim && j + curDim < inputSize; curDim += 1) {
                jac(i,j + curDim) = GT::at(tape.getGradient(output[i]), curDim);
              }
            }

            unseedGradient(tape, j, input, inputSize);
          }
        } else if(EvaluationType::Reverse == evalType) {

          for(size_t i = 0; i < outputSize; i += gradDim) {
            seedGradient(tape, i, output, outputSize);

            if(usePreAccEval) {
              tape.evaluatePreacc(end, start);
            } else {
              tape.evaluate(end, start);
            }

            for(size_t j = 0; j < inputSize; j += 1) {
              for(size_t curDim = 0; curDim < gradDim && i + curDim < outputSize; curDim += 1) {
                jac(i + curDim,j) = GT::at(tape.getGradient(input[j]), curDim);
                GT::at(tape.gradient(input[j]), curDim) = typename GT::Data();
              }
            }

            unseedGradient(tape, i, output, outputSize);

            if(!ZeroAdjointReverse) {
              tape.clearAdjoints(end, start);
            }
          }
        } else {
          CODI_EXCEPTION("Evaluation mode not implemented for preaccumulation. Mode is: %d.", (int)evalType);
        }
      }

      template<typename Jac>
      static CODI_INLINE void computeJacobian(
          Position const& start, Position const& end,
          GradientData const * input, size_t const inputSize,
          GradientData const * output, size_t const outputSize,
          Jac& jac)
      {
        computeJacobian(CoDiType::getGlobalTape(), start, end, input, inputSize, output, outputSize, jac);
      }

    private:
      static CODI_INLINE void seedGradient(Tape& tape, size_t const pos, GradientData const* values, const size_t size) {
        constexpr size_t gradDim = GT::getVectorSize();

        for(size_t curDim = 0; curDim < gradDim && pos + curDim < size; curDim += 1) {
          ENABLE_CHECK(ZeroChecks, tape.isActive(values[pos + curDim])) {
            GT::at(tape.gradient(values[pos + curDim]), curDim) = typename GT::Data(1.0);
          }
        }
      }

      static CODI_INLINE void unseedGradient(Tape& tape, size_t const pos, GradientData const* values, const size_t size) {
        constexpr size_t gradDim = GT::getVectorSize();

        for(size_t curDim = 0; curDim < gradDim && pos + curDim < size; curDim += 1) {
          ENABLE_CHECK(ZeroChecks, tape.isActive(values[pos + curDim])) {
            GT::at(tape.gradient(values[pos + curDim]), curDim) = typename GT::Data();
          }
        }
      }
  };


}
