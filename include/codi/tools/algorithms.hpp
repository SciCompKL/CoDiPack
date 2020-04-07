/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2020 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *     Max Sagebaum
 *     Tim Albring
 *     Johannes Blühdorn
 */


#pragma once

#include "../configure.h"
#include "../exceptions.hpp"
#include "../gradientTraits.hpp"
#include "data/jacobian.hpp"
#include "data/staticDummy.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Basic algorithms for tape evaluation in CoDiPack
   *
   * This class provides algorithms for:
   *  - Jacobian assembly
   *  - Hessian assembly
   *
   * All algorithms try to make the best choice for the evaluation mode, either forward or reverse. Which mode is
   * selected can be get in advance with the method Algorithms::getEvaluationChoice.
   *
   * \section FunctionDefinition Mathematical definitions of variables
   *
   * The function is defined by:
   *  \f[ y = f(x) \f]
   * where \f$ x \in \R^n \f$ is the input vector with the size \f$ n \f$ and \f$y \in \R^m \f$ is the output vector
   * with the size \f$ m \f$
   *
   * The Jacobian of \f$ f \f$ is defined by:
   *  \f[ J = \frac{\d f}{\d x} \in \R^{m \times n} \f]
   * The number of rows (\f$ m \f$) represents the number of output variables and the number of columns (\f$ n \f$)
   * represents the number of input variable. The derivative for the i-th output with respect to the j-th input is
   * represented by \f$ J_{i,j} \f$.
   *
   * The Jacobian of \f$ f \f$ is defined by:
   *  \f[ H = \frac{\d^2 f}{\d^2 x} \in \R^{m \times n \times n} \f]
   * The first dimension (\f$ m \f$) represents the number of output variables, the second and third dimension (\f$ n \f$) represents the
   * first and second derivative with respect to the input variables.
   * The second derivative for the i-th output with respect to the j-th and k-th input is
   * represented by \f$ H_{i,j,k} \f$.
   *
   * @tparam CoDiType      An ActiveReal type that has a tape which implements the ReverseTapeInterface
   * @tparam ActiveChecks  If activity checks for the seeding of gradient data should be performed. [Default: true]
   */
  template<typename CoDiType, bool ActiveChecks = true>
  struct Algorithms {

      using Tape = typename CoDiType::TapeType; /**< Tape type of the active real*/
      using Position = typename Tape::Position; /**< Position type of the tape */
      using Real = typename CoDiType::Real; /**< Inner type of the active real */
      using GradientData = typename CoDiType::GradientData; /**< Gradient data (identifier) of the tape */
      using GradientValue = typename CoDiType::GradientValue; /**< Gradient value type. */

      using GT = GradientValueTraits<GradientValue>; /**< Traits for the gradient value type. */

      /**
       * @brief Enumeration for the possible evaluation choices of a tape.
       */
      enum class EvaluationType {
          Forward,
          Reverse
      };

      /**
       * @brief Which path the algorithm will select for the evaluation.
       *
       * If the number of inputs are smaller than the outputs a forward evaluation will be chosen. Otherwise a
       * reverse evaluation.
       *
       * @param[in]  inputs  Number of inputs (n)
       * @param[in] outputs  Number of outputs (m)
       * @return
       */
      static CODI_INLINE EvaluationType getEvaluationChoice(size_t const inputs, size_t const outputs) {
        if(inputs <= outputs) {
          return EvaluationType::Forward;
        } else {
          return EvaluationType::Reverse;
        }
      }

      /**
       * @brief Compute the Jacobian with multiple tape sweeps.
       *
       * It has to hold start < end
       *
       * The algorithm expects that no gradient data has been seeded with non zero values.
       * After the return the algorithm ensures that all gradient data have zero values.
       *
       * @param[in,out]   tape  The tape with the recorded information for the derivative computation.
       * @param[in]      start  Starting position of the derivative information on the tape.
       * @param[in]        end  Ending position of the derivative information on the tape.
       * @param[in]      input  The gradient information for the input variables of f.
       * @param[in]  inputSize  The size of the input array.
       * @param[in]     output  The gradient information for the output variables of f.
       * @param[in] outputSize  The size of the output array.
       * @param[in,out]    jac  The Jacobian where the computed derivative information is stored.
       *
       * @tparam           Jac  Type of the Jacobian matrix. This type needs to implement the interface of the
       *                          Jacobian class.
       * @tparam usePreAccEval  If the specialized preaccumulation evaluation routines should be used.
       *                          For short tapes the specialized functions should be used. For large tapes the regular
       *                          functions are more performant.
       */
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
            setGradientOnGradientData(tape, j, input, inputSize, typename GT::Data(1.0));

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

            setGradientOnGradientData(tape, j, input, inputSize, typename GT::Data());
          }
        } else if(EvaluationType::Reverse == evalType) {

          for(size_t i = 0; i < outputSize; i += gradDim) {
            setGradientOnGradientData(tape, i, output, outputSize, typename GT::Data(1.0));

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

            setGradientOnGradientData(tape, i, output, outputSize, typename GT::Data());

            if(!ZeroAdjointReverse) {
              tape.clearAdjoints(end, start);
            }
          }
        } else {
          CODI_EXCEPTION("Evaluation mode not implemented for preaccumulation. Mode is: %d.", (int)evalType);
        }
      }

      /**
       * @brief Compute the Jacobian with multiple tape sweeps of the global tape.
       *
       * This method calls the more generalized version computeJacobian(Tape& tape, Position const& start, Position const& end, GradientData const * input, size_t const inputSize, GradientData const * output, size_t const outputSize, Jac& jac)
       * It has to hold start < end
       *
       * @param[in]      start  Starting position of the derivative information on the tape.
       * @param[in]        end  Ending position of the derivative information on the tape.
       * @param[in]      input  The gradient information for the input variables of f.
       * @param[in]  inputSize  The size of the input array.
       * @param[in]     output  The gradient information for the output variables of f.
       * @param[in] outputSize  The size of the output array.
       * @param[in,out]    jac  The Jacobian where the computed derivative information is stored.
       *
       * @tparam           Jac  Type of the Jacobian matrix. This type needs to implement the interface of the
       *                          Jacobian class.
       */
      template<typename Jac>
      static CODI_INLINE void computeJacobian(
          Position const& start, Position const& end,
          GradientData const * input, size_t const inputSize,
          GradientData const * output, size_t const outputSize,
          Jac& jac)
      {
        computeJacobian(CoDiType::getGlobalTape(), start, end, input, inputSize, output, outputSize, jac);
      }

      /**
       * @brief Compute the Hessian with multiple tape sweeps.
       *
       * This algorithm is only available if the tape implements the ReversePrimalValueTapeInterface. Since it performs
       * repeated primal evaluations to change original seeding of the tape. It also requires that the tape can compute
       * second order derivatives via a nested first order forward type.
       *
       * The algorithm expects that no gradient data has been seeded with non zero values.
       * It also expects that the current tape state was just recorded. This is, the primal values in the tape represent
       * the output values of f.
       * After the return the algorithm ensures that all gradient data have zero values.
       *
       * It has to hold start < end.
       *
       * @param[in,out]   tape  The tape with the recorded information for the derivative computation.
       * @param[in]      start  Starting position of the derivative information on the tape.
       * @param[in]        end  Ending position of the derivative information on the tape.
       * @param[in]      input  The gradient information for the input variables of f.
       * @param[in]  inputSize  The size of the input array.
       * @param[in]     output  The gradient information for the output variables of f.
       * @param[in] outputSize  The size of the output array.
       * @param[in,out]    hes  The Hessian where the computed derivative information is stored.
       * @param[in,out]    jac  Optional Jacobian matrix. If set also the gradient values are extracted during the
       *                          Hessian evaluation.
       *
       * @tparam           Hes  Type of the Hessian matrix. This type needs to implement the interface of the
       *                          Hessian class.
       * @tparam           Jac  Type of the Jacobian matrix. The type needs to implement the interface of the
       *                          Jacobian class.
       */
      template<typename Hes, typename Jac = DummyJacobian>
      static CODI_INLINE void computeHessianPrimalValueTape(
          Tape& tape, Position const& start, Position const& end,
          GradientData const * input, size_t const inputSize,
          GradientData const * output, size_t const outputSize,
          Hes& hes, Jac& jac = StaticDummy<DummyJacobian>::dummy)
      {
        EvaluationType evalType = getEvaluationChoice(inputSize, outputSize);
        if(EvaluationType::Forward == evalType) {
          computeHessianPrimalValueTapeForward(tape, start, end, input, inputSize, output, outputSize, hes, jac);
        } else if(EvaluationType::Reverse == evalType) {
          computeHessianPrimalValueTapeReverse(tape, start, end, input, inputSize, output, outputSize, hes, jac);
        } else {
          CODI_EXCEPTION("Evaluation mode not implemented for preaccumulation. Mode is: %d.", (int)evalType);
        }
      }

      /**
       * @brief Forward version of the hessian computation.
       *
       * Two input variables are seeded with gradient information and then a forward evaluation is performed.
       * Before each evaluation the primal values of the tape are reverted to the start position.
       *
       * The algorithm exploits symmetry and will perform n * (n + 1) / 2 forward tape evaluation. Vector gradient
       * values for the first and second order derivatives will reduce the tape evaluations accordingly.
       *
       * \copydetails computeHessianPrimalValueTape
       */
      template<typename Hes, typename Jac = DummyJacobian>
      static CODI_INLINE void computeHessianPrimalValueTapeForward(
          Tape& tape, Position const& start, Position const& end,
          GradientData const * input, size_t const inputSize,
          GradientData const * output, size_t const outputSize,
          Hes& hes, Jac& jac = StaticDummy<DummyJacobian>::dummy)
      {
        using GT1st = GT;
        constexpr size_t gradDim1st = GT1st::getVectorSize();
        using GT2nd = GradientValueTraits<typename Real::GradientValue>;
        constexpr size_t gradDim2nd = GT2nd::getVectorSize();

        // Assume tape that was just recorded
        tape.revertPrimals(start);

        for(size_t j = 0; j < inputSize; j += gradDim2nd) {
          setGradient2ndOnGradientData(tape, j, input, inputSize, typename GT2nd::Data(1.0));

          // The k = j init is no problem, it will evaluated slightly more elements around the diagonal
          for(size_t k = j; k < inputSize; k += gradDim1st) {
            setGradientOnGradientData(tape, k, input, inputSize, typename GT::Data(1.0));

            tape.evaluateForward(start, end);

            for(size_t i = 0; i < outputSize; i += 1) {
              for(size_t vecPos1st = 0; vecPos1st < gradDim1st && k + vecPos1st < inputSize; vecPos1st += 1) {
                for(size_t vecPos2nd = 0; vecPos2nd < gradDim2nd && j + vecPos2nd < inputSize; vecPos2nd += 1) {
                  hes(i, j + vecPos2nd, k + vecPos1st) = GT2nd::at(GT1st::at(tape.getGradient(output[i]), vecPos1st).gradient(), vecPos2nd);
                  hes(i, k + vecPos1st, j + vecPos2nd) = hes(i, j + vecPos2nd, k + vecPos1st); // symmetry
                }
              }

              if(j == 0) {
                for(size_t vecPos = 0; vecPos < gradDim1st && k + vecPos < inputSize; vecPos += 1) {
                  jac(i, k + vecPos) = GT1st::at(tape.getGradient(output[i]), vecPos).value();
                }
              }
            }

            setGradientOnGradientData(tape, k, input, inputSize, typename GT::Data());
          }

          setGradient2ndOnGradientData(tape, j, input, inputSize, typename GT2nd::Data());
        }
      }

      /**
       * @brief Reverse version of the hessian computation.
       *
       * One input variable is seeded with gradient information and then a forward evaluation is performed.
       * Afterwards one output variable is seeded with gradient information and then a reverse evaluation is performed.
       *
       * The algorithm can not exploit symmetry and will perform n forward evaluation and n * m reverse evaluations.
       * Vector gradient values for the first and second order derivatives will reduce the tape evaluations accordingly.
       *
       * \copydetails computeHessianPrimalValueTape
       */
      template<typename Hes, typename Jac = DummyJacobian>
      static CODI_INLINE void computeHessianPrimalValueTapeReverse(
          Tape& tape, Position const& start, Position const& end,
          GradientData const * input, size_t const inputSize,
          GradientData const * output, size_t const outputSize,
          Hes& hes, Jac& jac = StaticDummy<DummyJacobian>::dummy)
      {
        using GT1st = GT;
        constexpr size_t gradDim1st = GT1st::getVectorSize();
        using GT2nd = GradientValueTraits<typename Real::GradientValue>;
        constexpr size_t gradDim2nd = GT2nd::getVectorSize();

        // Assume tape that was just recorded
        tape.revertPrimals(start);

        for(size_t j = 0; j < inputSize; j += gradDim2nd) {
          setGradient2ndOnGradientData(tape, j, input, inputSize, typename GT2nd::Data(1.0));

          // propagate the new derivative information
          tape.evaluatePrimal(start, end);

          for(size_t i = 0; i < outputSize; i += gradDim1st) {
            setGradientOnGradientData(tape, i, output, outputSize, typename GT::Data(1.0));

            // propagate the derivatives backward for second order derivatives
            tape.evaluatePreacc(end, start);

            for(size_t k = 0; k < inputSize; k += 1) {
              for(size_t vecPos1st = 0; vecPos1st < gradDim1st && i + vecPos1st < outputSize; vecPos1st += 1) {
                for(size_t vecPos2nd = 0; vecPos2nd < gradDim2nd && j + vecPos2nd < inputSize; vecPos2nd += 1) {
                  hes(i + vecPos1st, j + vecPos2nd, k) = GT2nd::at(GT1st::at(tape.gradient(input[k]), vecPos1st).gradient(), vecPos2nd);
                }
              }

              if(j == 0) {
                for(size_t vecPos1st = 0; vecPos1st < gradDim1st && i + vecPos1st < outputSize; vecPos1st += 1) {
                  jac(i + vecPos1st, k) = GT1st::at(tape.getGradient(input[k]), vecPos1st).value();
                }
              }

              tape.gradient(input[k]) = GradientValue();
            }

            setGradientOnGradientData(tape, i, output, outputSize, typename GT::Data());

            if(!ZeroAdjointReverse) {
              tape.clearAdjoints(end, start);
            }
          }

          setGradient2ndOnGradientData(tape, j, input, inputSize, typename GT2nd::Data());

          if(j + gradDim2nd < inputSize) {
            tape.revertPrimals(start);
          }
        }
      }

      /**
       * @brief Compute the Hessian with multiple tape recordings and sweeps.
       *
       * The algorithm will repeatedly evaluated the function func and record the evaluation on the global tape. It
       * expects that the tape is empty. The tape needs to be able to compute second order derivatives via a nested
       * first order forward type.
       *
       * After the return, the algorithm ensures that the tape is empty.
       *
       * It has to hold start < end.
       *
       * @param[in]           func  The function for the recording of the tape. It needs to be a function object that
       *                              will accept the call: func(input, output)
       * @param[in,out]      input  The input values for the function. The gradient seeding will be changed on these values.
       * @param[in,out]     output  The output values for the function. They will be overwritten in a call to func.
       * @param[in,out]    hes  The Hessian where the computed derivative information is stored.
       * @param[in,out]    jac  Optional Jacobian matrix. If set also the gradient values are extracted during the
       *                          Hessian evaluation.
       *
       * @tparam           Hes  Type of the Hessian matrix. This type needs to implement the interface of the
       *                          Hessian class.
       * @tparam           Jac  Type of the Jacobian matrix. The type needs to implement the interface of the
       *                          Jacobian class.
       */
      template<typename Func, typename VecIn, typename VecOut, typename Hes, typename Jac = DummyJacobian>
      static CODI_INLINE void computeHessian(Func func, VecIn& input, VecOut& output, Hes& hes, Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        EvaluationType evalType = getEvaluationChoice(input.size(), output.size());
        if(EvaluationType::Forward == evalType) {
          computeHessianForward(func, input, output, hes, jac);
        } else if(EvaluationType::Reverse == evalType) {
          computeHessianReverse(func, input, output, hes, jac);
        } else {
          CODI_EXCEPTION("Evaluation mode not implemented for preaccumulation. Mode is: %d.", (int)evalType);
        }
      }

      /**
       * @brief Forward version of the hessian computation with a function object.
       *
       * One input variable is seeded with gradient information and then a tape is recorded.
       * Afterwards a second input variable is seeded with gradient information and the tape is evaluated multiple times in
       * the forward mode.
       * Before each recording the global tape is reset.
       *
       * The algorithm will record n tapes and exploits symmetry for the reverse evaluation which will result in
       * n * (n + 1) / 2 forward tape evaluations. Vector gradient values for the first and second order derivatives
       * will reduce the tape evaluations accordingly.
       *
       * \copydetails computeHessian
       */
      template<typename Func, typename VecIn, typename VecOut, typename Hes, typename Jac = DummyJacobian>
      static CODI_INLINE void computeHessianForward(Func func, VecIn& input, VecOut& output, Hes& hes, Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        using GT1st = GT;
        constexpr size_t gradDim1st = GT1st::getVectorSize();
        using GT2nd = GradientValueTraits<typename Real::GradientValue>;
        constexpr size_t gradDim2nd = GT2nd::getVectorSize();

        Tape& tape = CoDiType::getGlobalTape();

        for(size_t j = 0; j < input.size(); j += gradDim2nd) {
          setGradient2ndOnCoDiValue(j, input.data(), input.size(), typename GT2nd::Data(1.0));

          // propagate the new derivative information
          recordTape(func, input, output);

          // The k = j init is no problem, it will evaluated slightly more elements around the diagonal
          for(size_t k = j; k < input.size(); k += gradDim1st) {
            setGradientOnCoDiValue(tape, k, input.data(), input.size(), typename GT::Data(1.0));

            // propagate the derivatives forward for second order derivatives
            tape.evaluateForwardPreacc(tape.getZeroPosition(), tape.getPosition());

            for(size_t i = 0; i < output.size(); i += 1) {
              for(size_t vecPos1st = 0; vecPos1st < gradDim1st && k + vecPos1st < input.size(); vecPos1st += 1) {
                for(size_t vecPos2nd = 0; vecPos2nd < gradDim2nd && j + vecPos2nd < input.size(); vecPos2nd += 1) {
                  hes(i, j + vecPos2nd, k + vecPos1st) = GT2nd::at(GT1st::at(tape.getGradient(output[i].getGradientData()), vecPos1st).gradient(), vecPos2nd);
                  hes(i, k + vecPos1st, j + vecPos2nd) = hes(i, j + vecPos2nd, k + vecPos1st); // symmetry
                }
              }

              if(j == 0) {
                for(size_t vecPos = 0; vecPos < gradDim1st && k + vecPos < input.size(); vecPos += 1) {
                  jac(i, k + vecPos) = GT1st::at(tape.getGradient(output[i].getGradientData()), vecPos).value();
                }
              }
            }

            setGradientOnCoDiValue(tape, k, input.data(), input.size(), typename GT::Data());
          }

          setGradient2ndOnCoDiValue(j, input.data(), input.size(), typename GT2nd::Data());

          tape.reset();
        }
      }

      /**
       * @brief Reverse version of the hessian computation with a function object.
       *
       * One input variable is seeded with gradient information and then a tape is recorded.
       * Afterwards an output variable is seeded with gradient information and the tape is evaluated once in the reverse
       * mode.
       * Before each recording the global tape is reset.
       *
       * The algorithm will record n tapes and perform m * n reverse tape evaluations.
       * Vector gradient values for the first and second order derivatives will reduce the tape evaluations accordingly.
       *
       * \copydetails computeHessian
       */
      template<typename Func, typename VecIn, typename VecOut, typename Hes, typename Jac = DummyJacobian>
      static CODI_INLINE void computeHessianReverse(Func func, VecIn& input, VecOut& output, Hes& hes, Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        using GT1st = GT;
        constexpr size_t gradDim1st = GT1st::getVectorSize();
        using GT2nd = GradientValueTraits<typename Real::GradientValue>;
        constexpr size_t gradDim2nd = GT2nd::getVectorSize();

        Tape& tape = CoDiType::getGlobalTape();

        for(size_t j = 0; j < input.size(); j += gradDim2nd) {
          setGradient2ndOnCoDiValue(j, input.data(), input.size(), typename GT2nd::Data(1.0));

          // propagate the new derivative information
          recordTape(func, input, output);

          for(size_t i = 0; i < output.size(); i += gradDim1st) {
            setGradientOnCoDiValue(tape, i, output.data(), output.size(), typename GT::Data(1.0));

            // propagate the derivatives backward for second order derivatives
            tape.evaluatePreacc(tape.getPosition(), tape.getZeroPosition());

            for(size_t k = 0; k < input.size(); k += 1) {
              for(size_t vecPos1st = 0; vecPos1st < gradDim1st && i + vecPos1st < output.size(); vecPos1st += 1) {
                for(size_t vecPos2nd = 0; vecPos2nd < gradDim2nd && j + vecPos2nd < input.size(); vecPos2nd += 1) {
                  hes(i + vecPos1st, j + vecPos2nd, k) = GT2nd::at(GT1st::at(tape.gradient(input[k].getGradientData()), vecPos1st).gradient(), vecPos2nd);
                }
              }

              if(j == 0) {
                for(size_t vecPos1st = 0; vecPos1st < gradDim1st && i + vecPos1st < output.size(); vecPos1st += 1) {
                  jac(i + vecPos1st, k) = GT1st::at(tape.getGradient(input[k].getGradientData()), vecPos1st).value();
                }
              }

              tape.gradient(input[k].getGradientData()) = GradientValue();
            }

            setGradientOnCoDiValue(tape, i, output.data(), output.size(), typename GT::Data());

            if(!ZeroAdjointReverse) {
              tape.clearAdjoints(tape.getPosition(), tape.getZeroPosition());
            }
          }

          setGradient2ndOnCoDiValue(j, input.data(), input.size(), typename GT2nd::Data());

          tape.reset();
        }
      }

    private:

      /**
       * @brief Set gradient values for a vector mode evaluation.
       *
       * @param[in,out] tape  The tape on which the seeding is set.
       * @param[in]      pos  The position in the value vector.
       * @param[in]   values  The value vector.
       * @param[in]     size  The size of the value vector.
       * @param[in]    value  The value which is set as gradient information.
       *
       * @tparam T  The type of the provided value.
       */
      template<typename T>
      static CODI_INLINE void setGradientOnGradientData(Tape& tape, size_t const pos, GradientData const* values, const size_t size, T value) {
        constexpr size_t gradDim = GT::getVectorSize();

        for(size_t curDim = 0; curDim < gradDim && pos + curDim < size; curDim += 1) {
          ENABLE_CHECK(ActiveChecks, tape.isActive(values[pos + curDim])) {
            GT::at(tape.gradient(values[pos + curDim]), curDim) = value;
          }
        }
      }

      /**
       * @brief Set gradient values for a vector mode evaluation on the second order gradient information.
       *
       * @param[in,out] tape  The tape on which the seeding is set.
       * @param[in]      pos  The position in the value vector.
       * @param[in]   values  The value vector.
       * @param[in]     size  The size of the value vector.
       * @param[in]    value  The value which is set as gradient information.
       *
       * @tparam T  The type of the provided value.
       */
      template<typename T>
      static CODI_INLINE void setGradient2ndOnGradientData(Tape& tape, size_t const pos, GradientData const* values, const size_t size, T value) {
        using GT2nd = GradientValueTraits<typename Real::GradientValue>;
        constexpr size_t gradDim2nd = GT2nd::getVectorSize();

        for(size_t curDim = 0; curDim < gradDim2nd && pos + curDim < size; curDim += 1) {
          GT2nd::at(tape.primalValue(values[pos + curDim]).gradient(), curDim) = value;
        }
      }

      /**
       * @brief Set gradient values for a vector mode evaluation.
       *
       * @param[in,out] tape  The tape on which the seeding is set.
       * @param[in]      pos  The position in the value vector.
       * @param[in]   values  The value vector.
       * @param[in]     size  The size of the value vector.
       * @param[in]    value  The value which is set as gradient information.
       *
       * @tparam T  The type of the provided value.
       */
      template<typename T>
      static CODI_INLINE void setGradientOnCoDiValue(Tape& tape, size_t const pos, CoDiType* values, const size_t size, T value) {
        constexpr size_t gradDim = GT::getVectorSize();

        for(size_t curDim = 0; curDim < gradDim && pos + curDim < size; curDim += 1) {
          ENABLE_CHECK(ActiveChecks, tape.isActive(values[pos + curDim].getGradientData())) {
            GT::at(tape.gradient(values[pos + curDim].getGradientData()), curDim) = value;
          }
        }
      }

      /**
       * @brief Set gradient values for a vector mode evaluation on the second order gradient information.
       *
       * @param[in]      pos  The position in the value vector.
       * @param[in]   values  The value vector.
       * @param[in]     size  The size of the value vector.
       * @param[in]    value  The value which is set as gradient information.
       *
       * @tparam T  The type of the provided value.
       */
      template<typename T>
      static CODI_INLINE void setGradient2ndOnCoDiValue(size_t const pos, CoDiType* values, const size_t size, T value) {
        using GT2nd = GradientValueTraits<typename Real::GradientValue>;
        constexpr size_t gradDim2nd = GT2nd::getVectorSize();

        for(size_t curDim = 0; curDim < gradDim2nd && pos + curDim < size; curDim += 1) {
          // No check required since this are forward types.
          GT2nd::at(values[pos + curDim].value().gradient(), curDim) = value;
        }
      }


      /**
       * @brief Record the tape from the evaluation of the function object.
       *
       * @param[in]       func  The function object for the evaluation. It needs to accept the call func(input, output).
       * @param[in,out]  input  The values are registered as inputs.
       * @param[in,out] output  The values are registered as outputs and overwritten by the call to the function object.
       *
       * @tparam   Func  Function object.
       * @tparam  VecIn  Type of the input vector. Needs to implement the std interface.
       * @tparam VecOut  Type of the output vector.  Needs to implement the std interface.
       */
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
