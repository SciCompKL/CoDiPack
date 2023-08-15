/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../config.h"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../misc/exceptions.hpp"
#include "../tapes/misc/tapeParameters.hpp"
#include "../traits/gradientTraits.hpp"
#include "data/dummy.hpp"
#include "data/jacobian.hpp"
#include "data/staticDummy.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Basic algorithms for tape evaluation in CoDiPack.
   *
   * This class provides algorithms for:
   *  - Jacobian assembly
   *  - Hessian assembly
   *
   * All algorithms try to make the best choice for the evaluation mode depending on the number of inputs and outputs,
   * either forward or reverse. Which mode is selected can be queried in advance with the method
   * Algorithms::getEvaluationChoice.
   *
   * See \ref sec_namingConventions for the naming conventions.
   *
   * Hessians have to implement the HessianInterface and Jacobians have to implement the JacobianInterface.
   *
   * @tparam T_Type  An ActiveReal type that has a tape which implements the ReverseTapeInterface
   * @tparam ActiveChecks  If activity checks for the seeding of gradient data should be performed. [Default: true]
   */
  template<typename T_Type, bool T_ActiveChecks = true>
  struct Algorithms {
    public:

      /// See Algorithms.
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);

      static bool constexpr ActiveChecks = T_ActiveChecks;  ///< See Algorithms.

      using Tape = typename Type::Tape;              ///< See LhsExpressionInterface.
      using Position = typename Tape::Position;      ///< See LhsExpressionInterface.
      using Real = typename Type::Real;              ///< See LhsExpressionInterface.
      using Identifier = typename Type::Identifier;  ///< See LhsExpressionInterface.
      using Gradient = typename Type::Gradient;      ///< See LhsExpressionInterface.

      using GT = GradientTraits::TraitsImplementation<Gradient>;  ///< Shortcut for traits of gradient.

      /// Evaluation modes for the derivative computation.
      enum class EvaluationType {
        Forward,
        Reverse
      };

      /// Returns the preferred evaluation mode depending on the number of inputs and outputs.
      /// If the number of inputs is smaller than the number of outputs, a forward evaluation will be chosen.
      /// Otherwise, a reverse evaluation is preferred.
      static CODI_INLINE EvaluationType getEvaluationChoice(size_t const inputs, size_t const outputs) {
        if (inputs <= outputs) {
          return EvaluationType::Forward;
        } else {
          return EvaluationType::Reverse;
        }
      }

      /**
       * @brief Compute the Jacobian with multiple tape sweeps.
       *
       * The following prerequisites must be met.
       * - It must hold start < end.
       * - The gradient for every identifier used in the tape section [start, end] is seeded with zero.
       * - Every input to the tape section [start, end] that is a computational dependency of the specified outputs is
       *   contained in the specified inputs.
       *
       * The function has the following behavior.
       * - After return, the gradient for every input identifier, output identifier and all identifiers used in the tape
       *   section [start, end] on a left hand side are zero.
       * - If an output identifier is specified multiple times
       *   - and the forward mode is used, every row of the Jacobian corresponding to that output identifier except for
       *     the last one is zero.
       *   - and the reverse mode is used, every row of the Jacobian corresponding to that output identifier is
       *     identical.
       * - If an input identifier is specified multiple times
       *   - and the forward mode is used, every column of the Jacobian corresponding to that input identifier is
       *     identical.
       *   - and the reverse mode is used, every column of the Jacobian corresponding to that input identifier except
       *     for the first one is zero.
       *
       * The following usage is recommended.
       * - There should be no duplicate identifiers among the specified inputs.
       * - There should be no duplicate identifiers among the specified outputs.
       *
       * The following usage is allowed.
       * - The specified input identifiers and output identifiers need not be disjoint.
       *
       * The algorithm conforms with the mechanism for mutual exclusion of adjoint vector usage and adjoint vector
       * reallocation and can therefore be applied in multithreaded taping.
       *
       * In the case of manual adjoints management, the caller is responsible for ensuring sufficient adjoint vector
       * size and for declaring usage of the adjoints, see codi::AdjointsManagement for details.
       *
       * #### Parameters
       * [in,out] __jac__  Has to implement JacobianInterface.
       */
      template<typename Jac, bool keepState = true>
      static CODI_INLINE void computeJacobian(Tape& tape, Position const& start, Position const& end,
                                              Identifier const* input, size_t const inputSize, Identifier const* output,
                                              size_t const outputSize, Jac& jac,
                                              AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        size_t constexpr gradDim = GT::dim;

        // internally, automatic management is implemented in an optimized way that uses manual management
        if (AdjointsManagement::Automatic == adjointsManagement) {
          tape.resizeAdjointVector();
          tape.beginUseAdjointVector();
        }

        EvaluationType evalType = getEvaluationChoice(inputSize, outputSize);
        if (EvaluationType::Forward == evalType) {
          for (size_t j = 0; j < inputSize; j += gradDim) {
            setGradientOnIdentifier(tape, j, input, inputSize, typename GT::Real(1.0), AdjointsManagement::Manual);

            if (keepState) {
              tape.evaluateForwardKeepState(start, end, AdjointsManagement::Manual);
            } else {
              tape.evaluateForward(start, end, AdjointsManagement::Manual);
            }

            for (size_t i = 0; i < outputSize; i += 1) {
              for (size_t curDim = 0; curDim < gradDim && j + curDim < inputSize; curDim += 1) {
                jac(outputSize - i - 1, j + curDim) =
                    GT::at(tape.getGradient(output[outputSize - i - 1], AdjointsManagement::Manual), curDim);
                if (Gradient() != output[i]) {
                  GT::at(tape.gradient(output[outputSize - i - 1], AdjointsManagement::Manual), curDim) =
                      typename GT::Real();
                }
              }
            }

            setGradientOnIdentifier(tape, j, input, inputSize, typename GT::Real(), AdjointsManagement::Manual);
          }

          tape.clearAdjoints(end, start, AdjointsManagement::Manual);

        } else if (EvaluationType::Reverse == evalType) {
          for (size_t i = 0; i < outputSize; i += gradDim) {
            setGradientOnIdentifier(tape, i, output, outputSize, typename GT::Real(1.0), AdjointsManagement::Manual);

            if (keepState) {
              tape.evaluateKeepState(end, start, AdjointsManagement::Manual);
            } else {
              tape.evaluate(end, start, AdjointsManagement::Manual);
            }

            for (size_t j = 0; j < inputSize; j += 1) {
              for (size_t curDim = 0; curDim < gradDim && i + curDim < outputSize; curDim += 1) {
                jac(i + curDim, j) = GT::at(tape.getGradient(input[j], AdjointsManagement::Manual), curDim);
                GT::at(tape.gradient(input[j], AdjointsManagement::Manual), curDim) = typename GT::Real();
              }
            }

            setGradientOnIdentifier(tape, i, output, outputSize, typename GT::Real(), AdjointsManagement::Manual);

            if (!Config::ReversalZeroesAdjoints) {
              tape.clearAdjoints(end, start, AdjointsManagement::Manual);
            }
          }
        } else {
          CODI_EXCEPTION("Evaluation mode not implemented. Mode is: %d.", (int)evalType);
        }

        if (AdjointsManagement::Automatic == adjointsManagement) {
          tape.endUseAdjointVector();
        }
      }

      // clang-format off
      /// \copybrief computeJacobian(Tape&, Position const&, Position const&, Identifier const*, size_t const, Identifier const*, size_t const, Jac& jac, AdjointsManagement)
      /// \n This method uses the global tape for the Jacobian evaluation.
      /// \copydetails computeJacobian(Tape&, Position const&, Position const&, Identifier const*, size_t const, Identifier const*, size_t const, Jac& jac, AdjointsManagement)
      // clang-format on
      template<typename Jac>
      static CODI_INLINE void computeJacobian(Position const& start, Position const& end, Identifier const* input,
                                              size_t const inputSize, Identifier const* output, size_t const outputSize,
                                              Jac& jac,
                                              AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        computeJacobian(Type::getTape(), start, end, input, inputSize, output, outputSize, jac, adjointsManagement);
      }

      /**
       * @brief Compute the Hessian with multiple tape sweeps.
       *
       * This algorithm is only available if the tape implements the PrimalEvaluationTapeInterface. It performs repeated
       * primal evaluations to change the original seeding of the tape. It also requires that the tape can compute
       * second order derivatives via a nested first order forward type.
       *
       * The algorithm expects that no gradients have been seeded with nonzero values.
       * It also expects that the current tape state was just recorded, that is, the primal values in the tape represent
       * the output values of f.
       * At return, it is ensured that all gradients have zero values.
       *
       * It has to hold start < end.
       *
       * #### Parameters
       * [in,out] __hes__  Has to implement HessianInterface. \n
       * [in,out] __jac__  Optional: Jacobian values are also extracted. Has to implement JacobianInterface.
       */
      template<typename Hes, typename Jac = DummyJacobian>
      static CODI_INLINE void computeHessianPrimalValueTape(Tape& tape, Position const& start, Position const& end,
                                                            Identifier const* input, size_t const inputSize,
                                                            Identifier const* output, size_t const outputSize, Hes& hes,
                                                            Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        EvaluationType evalType = getEvaluationChoice(inputSize, outputSize);
        if (EvaluationType::Forward == evalType) {
          computeHessianPrimalValueTapeForward(tape, start, end, input, inputSize, output, outputSize, hes, jac);
        } else if (EvaluationType::Reverse == evalType) {
          computeHessianPrimalValueTapeReverse(tape, start, end, input, inputSize, output, outputSize, hes, jac);
        } else {
          CODI_EXCEPTION("Evaluation mode not implemented. Mode is: %d.", (int)evalType);
        }
      }

      /**
       * @brief Forward version of the Hessian computation.
       *
       * Two input variables are seeded with gradient information and then a forward evaluation is performed.
       * Before each evaluation, the primal values of the tape are reverted to the start position.
       *
       * The algorithm exploits symmetry and will perform n * (n + 1) / 2 forward tape evaluation. Vector gradient
       * values for the first and second order derivatives will reduce the tape evaluations accordingly.
       *
       * \copydetails computeHessianPrimalValueTape
       */
      template<typename Hes, typename Jac = DummyJacobian>
      static CODI_INLINE void computeHessianPrimalValueTapeForward(Tape& tape, Position const& start,
                                                                   Position const& end, Identifier const* input,
                                                                   size_t const inputSize, Identifier const* output,
                                                                   size_t const outputSize, Hes& hes,
                                                                   Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        using GT1st = GT;
        size_t constexpr gradDim1st = GT1st::dim;
        using GT2nd = GradientTraits::TraitsImplementation<CODI_DD(typename Real::Gradient, double)>;
        size_t constexpr gradDim2nd = GT2nd::dim;

        // Assume that the tape was just recorded.
        tape.revertPrimals(start);

        for (size_t j = 0; j < inputSize; j += gradDim2nd) {
          setGradient2ndOnIdentifier(tape, j, input, inputSize, typename GT2nd::Real(1.0));

          // The k = j init is no problem, it will evaluate slightly more elements around the diagonal.
          for (size_t k = j; k < inputSize; k += gradDim1st) {
            setGradientOnIdentifier(tape, k, input, inputSize, typename GT1st::Real(1.0));

            tape.evaluateForward(start, end);

            for (size_t i = 0; i < outputSize; i += 1) {
              for (size_t vecPos1st = 0; vecPos1st < gradDim1st && k + vecPos1st < inputSize; vecPos1st += 1) {
                for (size_t vecPos2nd = 0; vecPos2nd < gradDim2nd && j + vecPos2nd < inputSize; vecPos2nd += 1) {
                  hes(i, j + vecPos2nd, k + vecPos1st) =
                      GT2nd::at(GT1st::at(tape.getGradient(output[i]), vecPos1st).gradient(), vecPos2nd);
                  hes(i, k + vecPos1st, j + vecPos2nd) = hes(i, j + vecPos2nd, k + vecPos1st);  // Symmetry
                }
              }

              if (j == 0) {
                for (size_t vecPos = 0; vecPos < gradDim1st && k + vecPos < inputSize; vecPos += 1) {
                  jac(i, k + vecPos) = GT1st::at(tape.getGradient(output[i]), vecPos).value();
                }
              }
            }

            setGradientOnIdentifier(tape, k, input, inputSize, typename GT1st::Real());
          }

          setGradient2ndOnIdentifier(tape, j, input, inputSize, typename GT2nd::Real());
        }
      }

      /**
       * @brief Reverse version of the Hessian computation.
       *
       * One input variable is seeded with gradient information and then a forward evaluation is performed.
       * Afterwards, one output variable is seeded with gradient information and then a reverse evaluation is performed.
       *
       * The algorithm cannot exploit symmetry and will perform n forward evaluation and n * m reverse evaluations.
       * Vector gradient values for the first and second order derivatives will reduce the tape evaluations accordingly.
       *
       * \copydetails computeHessianPrimalValueTape
       */
      template<typename Hes, typename Jac = DummyJacobian>
      static CODI_INLINE void computeHessianPrimalValueTapeReverse(Tape& tape, Position const& start,
                                                                   Position const& end, Identifier const* input,
                                                                   size_t const inputSize, Identifier const* output,
                                                                   size_t const outputSize, Hes& hes,
                                                                   Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        using GT1st = GT;
        size_t constexpr gradDim1st = GT1st::dim;
        using GT2nd = GradientTraits::TraitsImplementation<CODI_DD(typename Real::Gradient, double)>;
        size_t constexpr gradDim2nd = GT2nd::dim;

        // Assume that the tape was just recorded.
        tape.revertPrimals(start);

        for (size_t j = 0; j < inputSize; j += gradDim2nd) {
          setGradient2ndOnIdentifier(tape, j, input, inputSize, typename GT2nd::Real(1.0));

          // Propagate the new derivative information.
          tape.evaluatePrimal(start, end);

          for (size_t i = 0; i < outputSize; i += gradDim1st) {
            setGradientOnIdentifier(tape, i, output, outputSize, typename GT1st::Real(1.0));

            // Propagate the derivatives backward for second order derivatives.
            tape.evaluateKeepState(end, start);

            for (size_t k = 0; k < inputSize; k += 1) {
              for (size_t vecPos1st = 0; vecPos1st < gradDim1st && i + vecPos1st < outputSize; vecPos1st += 1) {
                for (size_t vecPos2nd = 0; vecPos2nd < gradDim2nd && j + vecPos2nd < inputSize; vecPos2nd += 1) {
                  hes(i + vecPos1st, j + vecPos2nd, k) =
                      GT2nd::at(GT1st::at(tape.gradient(input[k]), vecPos1st).gradient(), vecPos2nd);
                }
              }

              if (j == 0) {
                for (size_t vecPos1st = 0; vecPos1st < gradDim1st && i + vecPos1st < outputSize; vecPos1st += 1) {
                  jac(i + vecPos1st, k) = GT1st::at(tape.getGradient(input[k]), vecPos1st).value();
                }
              }

              tape.gradient(input[k]) = Gradient();
            }

            setGradientOnIdentifier(tape, i, output, outputSize, typename GT1st::Real());

            if (!Config::ReversalZeroesAdjoints) {
              tape.clearAdjoints(end, start);
            }
          }

          setGradient2ndOnIdentifier(tape, j, input, inputSize, typename GT2nd::Real());

          if (j + gradDim2nd < inputSize) {
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
       * #### Parameters
       * [in]     __func__  The function for the recording of the tape. It needs to be a function object that
       *                         will accept the call: func(input, output) \n
       * [in,out] __hes__  Has to implement HessianInterface. \n
       * [in,out] __jac__  Optional: Jacobian values are also extracted. Has to implement JacobianInterface.
       */
      template<typename Func, typename VecIn, typename VecOut, typename Hes, typename Jac = DummyJacobian>
      static CODI_INLINE void computeHessian(Func func, VecIn& input, VecOut& output, Hes& hes,
                                             Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        EvaluationType evalType = getEvaluationChoice(input.size(), output.size());
        if (EvaluationType::Forward == evalType) {
          computeHessianForward(func, input, output, hes, jac);
        } else if (EvaluationType::Reverse == evalType) {
          computeHessianReverse(func, input, output, hes, jac);
        } else {
          CODI_EXCEPTION("Evaluation mode not implemented. Mode is: %d.", (int)evalType);
        }
      }

      /**
       * @brief Forward version of the Hessian computation with a function object.
       *
       * One input variable is seeded with gradient information and then a tape is recorded.
       * Afterwards a second input variable is seeded with gradient information and the tape is evaluated multiple times
       * in the forward mode. Before each recording the global tape is reset.
       *
       * The algorithm will record n tapes and exploits symmetry for the forward evaluation which will result in
       * n * (n + 1) / 2 forward tape evaluations. Vector gradient values for the first and second order derivatives
       * will reduce the tape evaluations accordingly.
       *
       * \copydetails computeHessian
       */
      template<typename Func, typename VecIn, typename VecOut, typename Hes, typename Jac = DummyJacobian>
      static CODI_INLINE void computeHessianForward(Func func, VecIn& input, VecOut& output, Hes& hes,
                                                    Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        using GT1st = GT;
        size_t constexpr gradDim1st = GT1st::dim;
        using GT2nd = GradientTraits::TraitsImplementation<CODI_DD(typename Real::Gradient, double)>;
        size_t constexpr gradDim2nd = GT2nd::dim;

        Tape& tape = Type::getTape();

        for (size_t j = 0; j < input.size(); j += gradDim2nd) {
          setGradient2ndOnCoDiValue(j, input.data(), input.size(), typename GT2nd::Real(1.0));

          // Propagate the new derivative information.
          recordTape(func, input, output);

          // The k = j init is no problem, it will evaluate slightly more elements around the diagonal.
          for (size_t k = j; k < input.size(); k += gradDim1st) {
            setGradientOnCoDiValue(tape, k, input.data(), input.size(), typename GT1st::Real(1.0));

            // Propagate the derivatives forward for second order derivatives.
            tape.evaluateForwardKeepState(tape.getZeroPosition(), tape.getPosition());

            for (size_t i = 0; i < output.size(); i += 1) {
              for (size_t vecPos1st = 0; vecPos1st < gradDim1st && k + vecPos1st < input.size(); vecPos1st += 1) {
                for (size_t vecPos2nd = 0; vecPos2nd < gradDim2nd && j + vecPos2nd < input.size(); vecPos2nd += 1) {
                  hes(i, j + vecPos2nd, k + vecPos1st) = GT2nd::at(
                      GT1st::at(tape.getGradient(output[i].getIdentifier()), vecPos1st).gradient(), vecPos2nd);
                  hes(i, k + vecPos1st, j + vecPos2nd) = hes(i, j + vecPos2nd, k + vecPos1st);  // Symmetry
                }
              }

              if (j == 0) {
                for (size_t vecPos = 0; vecPos < gradDim1st && k + vecPos < input.size(); vecPos += 1) {
                  jac(i, k + vecPos) = GT1st::at(tape.getGradient(output[i].getIdentifier()), vecPos).value();
                }
              }
            }

            setGradientOnCoDiValue(tape, k, input.data(), input.size(), typename GT1st::Real());
          }

          setGradient2ndOnCoDiValue(j, input.data(), input.size(), typename GT2nd::Real());

          tape.reset();
        }
      }

      /**
       * @brief Reverse version of the Hessian computation with a function object.
       *
       * One input variable is seeded with gradient information and then a tape is recorded.
       * Afterwards an output variable is seeded with gradient information and the tape is evaluated once in the reverse
       * mode.
       * Before each recording, the global tape is reset.
       *
       * The algorithm will record n tapes and perform m * n reverse tape evaluations.
       * Vector gradient values for the first and second order derivatives will reduce the tape evaluations accordingly.
       *
       * \copydetails computeHessian
       */
      template<typename Func, typename VecIn, typename VecOut, typename Hes, typename Jac = DummyJacobian>
      static CODI_INLINE void computeHessianReverse(Func func, VecIn& input, VecOut& output, Hes& hes,
                                                    Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        using GT1st = GT;
        size_t constexpr gradDim1st = GT1st::dim;
        using GT2nd = GradientTraits::TraitsImplementation<CODI_DD(typename Real::Gradient, double)>;
        size_t constexpr gradDim2nd = GT2nd::dim;

        Tape& tape = Type::getTape();

        for (size_t j = 0; j < input.size(); j += gradDim2nd) {
          setGradient2ndOnCoDiValue(j, input.data(), input.size(), typename GT2nd::Real(1.0));

          // Propagate the new derivative information.
          recordTape(func, input, output);

          for (size_t i = 0; i < output.size(); i += gradDim1st) {
            setGradientOnCoDiValue(tape, i, output.data(), output.size(), typename GT1st::Real(1.0));

            // Propagate the derivatives backward for second order derivatives.
            tape.evaluateKeepState(tape.getPosition(), tape.getZeroPosition());

            for (size_t k = 0; k < input.size(); k += 1) {
              for (size_t vecPos1st = 0; vecPos1st < gradDim1st && i + vecPos1st < output.size(); vecPos1st += 1) {
                for (size_t vecPos2nd = 0; vecPos2nd < gradDim2nd && j + vecPos2nd < input.size(); vecPos2nd += 1) {
                  hes(i + vecPos1st, j + vecPos2nd, k) =
                      GT2nd::at(GT1st::at(tape.gradient(input[k].getIdentifier()), vecPos1st).gradient(), vecPos2nd);
                }
              }

              if (j == 0) {
                for (size_t vecPos1st = 0; vecPos1st < gradDim1st && i + vecPos1st < output.size(); vecPos1st += 1) {
                  jac(i + vecPos1st, k) = GT1st::at(tape.getGradient(input[k].getIdentifier()), vecPos1st).value();
                }
              }

              tape.gradient(input[k].getIdentifier()) = Gradient();
            }

            setGradientOnCoDiValue(tape, i, output.data(), output.size(), typename GT1st::Real());

            if (!Config::ReversalZeroesAdjoints) {
              tape.clearAdjoints(tape.getPosition(), tape.getZeroPosition());
            }
          }

          setGradient2ndOnCoDiValue(j, input.data(), input.size(), typename GT2nd::Real());

          tape.reset();
        }
      }

    private:

      /**
       * @brief Sets the gradient for vector modes. Seeds the next GT::dim dimensions.
       *
       * Does not perform bounds checking for the gradient access.
       * Declares usage of the adjoint vector, see DataManagementTapeInterface.
       */
      template<typename T>
      static CODI_INLINE void setGradientOnIdentifier(
          Tape& tape, size_t const pos, Identifier const* identifiers, size_t const size, T value,
          AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        size_t constexpr gradDim = GT::dim;

        for (size_t curDim = 0; curDim < gradDim && pos + curDim < size; curDim += 1) {
          if (CODI_ENABLE_CHECK(ActiveChecks, 0 != identifiers[pos + curDim])) {
            GT::at(tape.gradient(identifiers[pos + curDim], adjointsManagement), curDim) = value;
          }
        }
      }

      /// Sets the gradient for 2nd order vector modes. Seeds the next GT2nd:dim dimensions.
      template<typename T>
      static CODI_INLINE void setGradient2ndOnIdentifier(Tape& tape, size_t const pos, Identifier const* identifiers,
                                                         size_t const size, T value) {
        using GT2nd = GradientTraits::TraitsImplementation<CODI_DD(typename Real::Gradient, double)>;
        size_t constexpr gradDim2nd = GT2nd::dim;

        for (size_t curDim = 0; curDim < gradDim2nd && pos + curDim < size; curDim += 1) {
          // No activity check on the identifier required since forward types are used.
          GT2nd::at(tape.primal(identifiers[pos + curDim]).gradient(), curDim) = value;
        }
      }

      /**
       * @brief Sets the gradient for 1st order vector modes. Seeds the next GT:dim dimensions.
       *
       * Does not perform bounds checking for the gradient access.
       * Declares usage of the adjoint vector, see DataManagementTapeInterface.
       */
      template<typename T>
      static CODI_INLINE void setGradientOnCoDiValue(Tape& tape, size_t const pos, Type* identifiers, size_t const size,
                                                     T value) {
        size_t constexpr gradDim = GT::dim;

        for (size_t curDim = 0; curDim < gradDim && pos + curDim < size; curDim += 1) {
          if (CODI_ENABLE_CHECK(ActiveChecks, 0 != identifiers[pos + curDim].getIdentifier())) {
            GT::at(tape.gradient(identifiers[pos + curDim].getIdentifier()), curDim) = value;
          }
        }
      }

      /// Sets the gradient for 2nd order vector modes. Seeds the next GT2nd:dim dimensions.
      template<typename T>
      static CODI_INLINE void setGradient2ndOnCoDiValue(size_t const pos, Type* identifiers, size_t const size,
                                                        T value) {
        using GT2nd = GradientTraits::TraitsImplementation<CODI_DD(typename Real::Gradient, double)>;
        size_t constexpr gradDim2nd = GT2nd::dim;

        for (size_t curDim = 0; curDim < gradDim2nd && pos + curDim < size; curDim += 1) {
          // No activity check on the identifier required since forward types are used.
          GT2nd::at(identifiers[pos + curDim].value().gradient(), curDim) = value;
        }
      }

      /// Record an evalaution of the function.
      template<typename Func, typename VecIn, typename VecOut>
      static CODI_INLINE void recordTape(Func func, VecIn& input, VecOut& output) {
        Tape& tape = Type::getTape();
        tape.setActive();
        for (size_t curIn = 0; curIn < input.size(); curIn += 1) {
          tape.registerInput(input[curIn]);
        }

        func(input, output);

        for (size_t curOut = 0; curOut < output.size(); curOut += 1) {
          tape.registerOutput(output[curOut]);
        }
        tape.setPassive();
      }
  };

}
