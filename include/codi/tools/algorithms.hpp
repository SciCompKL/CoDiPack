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

      template<typename Hes>
      static CODI_INLINE void computeHessianPrimalValueTape(
          Tape& tape, Position const& start, Position const& end,
          GradientData const * input, size_t const inputSize,
          GradientData const * output, size_t const outputSize,
          Hes& hes)
      {
        EvaluationType evalType = getEvaluationChoice(inputSize, outputSize);
        if(EvaluationType::Forward == evalType) {
          computeHessianPrimalValueTapeForward(tape, start, end, input, inputSize, output, outputSize, hes);
        } else if(EvaluationType::Reverse == evalType) {
          computeHessianPrimalValueTapeReverse(tape, start, end, input, inputSize, output, outputSize, hes);
        } else {
          CODI_EXCEPTION("Evaluation mode not implemented for preaccumulation. Mode is: %d.", (int)evalType);
        }
      }

      template<typename Hes>
      static CODI_INLINE void computeHessianPrimalValueTapeForward(
          Tape& tape, Position const& start, Position const& end,
          GradientData const * input, size_t const inputSize,
          GradientData const * output, size_t const outputSize,
          Hes& hes)
      {
        using GT1st = GT;
        constexpr size_t gradDim1st = GT1st::getVectorSize();
        using GT2nd = GradientValueTraits<typename Real::GradientValue>;
        constexpr size_t gradDim2nd = GT2nd::getVectorSize();

        // Assume tape that was just recorded
        tape.revertPrimals(start);

        for(size_t j = 0; j < inputSize; j += gradDim2nd) {
          for(size_t curDim = 0; curDim < gradDim2nd && j + curDim < inputSize; curDim += 1) {
            GT2nd::at(tape.primalValue(input[j + curDim]).gradient(), curDim) = typename GT2nd::Data(1.0);
          }

          // The k = j init is no prolbem, it will evaluated slightly more elements around the diagonal
          for(size_t k = j; k < inputSize; k += gradDim1st) {
            seedGradient(tape, k, input, inputSize);

            tape.evaluateForward(start, end);

            for(size_t i = 0; i < outputSize; i += 1) {
              for(size_t vecPos1st = 0; vecPos1st < gradDim1st && k + vecPos1st < inputSize; vecPos1st += 1) {
                for(size_t vecPos2nd = 0; vecPos2nd < gradDim2nd && j + vecPos2nd < inputSize; vecPos2nd += 1) {
                  hes(i, j + vecPos2nd, k + vecPos1st) = GT2nd::at(GT1st::at(tape.getGradient(output[i]), vecPos1st).gradient(), vecPos2nd);
                  hes(i, k + vecPos1st, j + vecPos2nd) = hes(i, j + vecPos2nd, k + vecPos1st); // symmetry
                }
              }
            }

            unseedGradient(tape, k, input, inputSize);
          }

          for(size_t curDim = 0; curDim < gradDim2nd && j + curDim < inputSize; curDim += 1) {
            GT2nd::at(tape.primalValue(input[j + curDim]).gradient(), curDim) = typename GT2nd::Data();
          }
        }
      }

      template<typename Hes>
      static CODI_INLINE void computeHessianPrimalValueTapeReverse(
          Tape& tape, Position const& start, Position const& end,
          GradientData const * input, size_t const inputSize,
          GradientData const * output, size_t const outputSize,
          Hes& hes)
      {
        using GT1st = GT;
        constexpr size_t gradDim1st = GT1st::getVectorSize();
        using GT2nd = GradientValueTraits<typename Real::GradientValue>;
        constexpr size_t gradDim2nd = GT2nd::getVectorSize();

        // Assume tape that was just recorded
        tape.revertPrimals(start);

        for(size_t j = 0; j < inputSize; j += gradDim2nd) {
          for(size_t curDim = 0; curDim < gradDim2nd && j + curDim < inputSize; curDim += 1) {
            GT2nd::at(tape.primalValue(input[j + curDim]).gradient(), curDim) = typename GT2nd::Data(1.0);
          }

          // propagate the new derivative information
          tape.evaluatePrimal(start, end);

          for(size_t i = 0; i < outputSize; i += gradDim1st) {
            seedGradient(tape, i, output, outputSize);

            // propaget the derivatives backward for second order derivatives
            tape.evaluatePreacc(end, start);

            for(size_t k = 0; k < inputSize; k += 1) {
              for(size_t vecPos1st = 0; vecPos1st < gradDim1st && i + vecPos1st < outputSize; vecPos1st += 1) {
                for(size_t vecPos2nd = 0; vecPos2nd < gradDim2nd && j + vecPos2nd < inputSize; vecPos2nd += 1) {
                  hes(i + vecPos1st, j + vecPos2nd, k) = GT2nd::at(GT1st::at(tape.gradient(input[k]), vecPos1st).gradient(), vecPos2nd);
                }
              }

              tape.gradient(input[k]) = GradientValue();
            }

            unseedGradient(tape, i, output, outputSize);

            if(!ZeroAdjointReverse) {
              tape.clearAdjoints(end, start);
            }
          }

          for(size_t curDim = 0; curDim < gradDim2nd && j + curDim < inputSize; curDim += 1) {
            GT2nd::at(tape.primalValue(input[j + curDim]).gradient(), curDim) = typename GT2nd::Data();
          }

          if(j + gradDim2nd < inputSize) {
            tape.revertPrimals(start);
          }
        }
      }

      template<typename Func, typename VecIn, typename VecOut, typename Hes>
      static CODI_INLINE void computeHessian(Func func, VecIn& input, VecOut& output, Hes& hes) {
        EvaluationType evalType = getEvaluationChoice(input.size(), output.size());
        if(EvaluationType::Forward == evalType) {
          computeHessianForward(func, input, output, hes);
        } else if(EvaluationType::Reverse == evalType) {
          computeHessianReverse(func, input, output, hes);
        } else {
          CODI_EXCEPTION("Evaluation mode not implemented for preaccumulation. Mode is: %d.", (int)evalType);
        }
      }

      template<typename Func, typename VecIn, typename VecOut, typename Hes>
      static CODI_INLINE void computeHessianForward(Func func, VecIn& input, VecOut& output, Hes& hes) {
        using GT1st = GT;
        constexpr size_t gradDim1st = GT1st::getVectorSize();
        using GT2nd = GradientValueTraits<typename Real::GradientValue>;
        constexpr size_t gradDim2nd = GT2nd::getVectorSize();

        Tape& tape = CoDiType::getGlobalTape();

        for(size_t j = 0; j < input.size(); j += gradDim2nd) {
          for(size_t curDim = 0; curDim < gradDim2nd && j + curDim < input.size(); curDim += 1) {
            GT2nd::at(input[j + curDim].value().gradient(), curDim)= typename GT2nd::Data(1.0);
          }

          // propagate the new derivative information
          recordTape(func, input, output);

          // The k = j init is no prolbem, it will evaluated slightly more elements around the diagonal
          for(size_t k = j; k < input.size(); k += gradDim1st) {
            setGradientOnValue1stOrder(tape, k, input.data(), input.size(), typename GT::Data(1.0));

            // propaget the derivatives backward for second order derivatives
            tape.evaluateForwardPreacc(tape.getZeroPosition(), tape.getPosition());

            for(size_t i = 0; i < output.size(); i += 1) {
              for(size_t vecPos1st = 0; vecPos1st < gradDim1st && k + vecPos1st < input.size(); vecPos1st += 1) {
                for(size_t vecPos2nd = 0; vecPos2nd < gradDim2nd && j + vecPos2nd < input.size(); vecPos2nd += 1) {
                  hes(i, j + vecPos2nd, k + vecPos1st) = GT2nd::at(GT1st::at(tape.getGradient(output[i].getGradientData()), vecPos1st).gradient(), vecPos2nd);
                  hes(i, k + vecPos1st, j + vecPos2nd) = hes(i, j + vecPos2nd, k + vecPos1st); // symmetry
                }
              }
            }

            setGradientOnValue1stOrder(tape, k, input.data(), input.size(), typename GT::Data());
          }

          for(size_t curDim = 0; curDim < gradDim2nd && j + curDim < input.size(); curDim += 1) {
            GT2nd::at(input[j + curDim].value().gradient(), curDim)= typename GT2nd::Data();
          }

          tape.reset();
        }
      }

      template<typename Func, typename VecIn, typename VecOut, typename Hes>
      static CODI_INLINE void computeHessianReverse(Func func, VecIn& input, VecOut& output, Hes& hes) {
        using GT1st = GT;
        constexpr size_t gradDim1st = GT1st::getVectorSize();
        using GT2nd = GradientValueTraits<typename Real::GradientValue>;
        constexpr size_t gradDim2nd = GT2nd::getVectorSize();

        Tape& tape = CoDiType::getGlobalTape();

        for(size_t j = 0; j < input.size(); j += gradDim2nd) {
          for(size_t curDim = 0; curDim < gradDim2nd && j + curDim < input.size(); curDim += 1) {
            GT2nd::at(input[j + curDim].value().gradient(), curDim)= typename GT2nd::Data(1.0);
          }

          // propagate the new derivative information
          recordTape(func, input, output);

          for(size_t i = 0; i < output.size(); i += gradDim1st) {
            setGradientOnValue1stOrder(tape, i, output.data(), output.size(), typename GT::Data(1.0));

            // propaget the derivatives backward for second order derivatives
            tape.evaluatePreacc(tape.getPosition(), tape.getZeroPosition());

            for(size_t k = 0; k < input.size(); k += 1) {
              for(size_t vecPos1st = 0; vecPos1st < gradDim1st && i + vecPos1st < output.size(); vecPos1st += 1) {
                for(size_t vecPos2nd = 0; vecPos2nd < gradDim2nd && j + vecPos2nd < input.size(); vecPos2nd += 1) {
                  hes(i + vecPos1st, j + vecPos2nd, k) = GT2nd::at(GT1st::at(tape.gradient(input[k].getGradientData()), vecPos1st).gradient(), vecPos2nd);
                }
              }

              tape.gradient(input[k].getGradientData()) = GradientValue();
            }

            setGradientOnValue1stOrder(tape, i, output.data(), output.size(), typename GT::Data());

            if(!ZeroAdjointReverse) {
              tape.clearAdjoints(tape.getPosition(), tape.getZeroPosition());
            }
          }

          for(size_t curDim = 0; curDim < gradDim2nd && j + curDim < input.size(); curDim += 1) {
            GT2nd::at(input[j + curDim].value().gradient(), curDim)= typename GT2nd::Data();
          }

          tape.reset();
        }
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

      template<typename T>
      static CODI_INLINE void setGradientOnValue1stOrder(Tape& tape, size_t const pos, CoDiType* values, const size_t size, T value) {
        constexpr size_t gradDim = GT::getVectorSize();

        for(size_t curDim = 0; curDim < gradDim && pos + curDim < size; curDim += 1) {
          ENABLE_CHECK(ZeroChecks, tape.isActive(values[pos + curDim].getGradientData())) {
            GT::at(tape.gradient(values[pos + curDim].getGradientData()), curDim) = value;
          }
        }
      }


      template<typename Func, typename VecIn, typename VecOut>
      static CODI_INLINE void recordTape(Func func, VecIn& input, VecOut& output) {

        Tape& tape = CoDiType::getGlobalTape();
        tape.setActive();
        for(size_t curIn = 0; curIn < input.size(); curIn += 1) {
          tape.registerInput(input[curIn]);
        }

        func(input, output);

        for(size_t curOut = 0; curOut < output.size(); curOut += 1) {
          tape.registerOutput(output[curOut]);
        }
        tape.setPassive();
      }
  };


}
