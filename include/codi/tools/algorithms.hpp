#pragma once

#include "../aux/exceptions.hpp"
#include "../config.h"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../traits/gradientTraits.hpp"
#include "data/dummy.hpp"
#include "data/staticDummy.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Type, bool _ActiveChecks = true>
  struct Algorithms {
    public:

      using Type = CODI_DECLARE_DEFAULT(_Type,
                                        CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

      static bool constexpr ActiveChecks = _ActiveChecks;

      using Tape = typename Type::Tape;
      using Position = typename Tape::Position;
      using Real = typename Type::Real;
      using Identifier = typename Type::Identifier;
      using Gradient = typename Type::Gradient;

      using GT = GradientTraits::TraitsImplementation<Gradient>;

      enum class EvaluationType
      {
        Forward,
        Reverse
      };

      static CODI_INLINE EvaluationType getEvaluationChoice(size_t const inputs, size_t const outputs) {
        if (inputs <= outputs) {
          return EvaluationType::Forward;
        } else {
          return EvaluationType::Reverse;
        }
      }

      template<typename Jac, bool keepState = true>
      static CODI_INLINE void computeJacobian(Tape& tape, Position const& start, Position const& end,
                                              Identifier const* input, size_t const inputSize, Identifier const* output,
                                              size_t const outputSize, Jac& jac) {
        size_t constexpr gradDim = GT::dim;

        EvaluationType evalType = getEvaluationChoice(inputSize, outputSize);
        if (EvaluationType::Forward == evalType) {
          for (size_t j = 0; j < inputSize; j += gradDim) {
            setGradientOnIdentifier(tape, j, input, inputSize, typename GT::Real(1.0));

            if (keepState) {
              tape.evaluateForwardKeepState(start, end);
            } else {
              tape.evaluateForward(start, end);
            }

            for (size_t i = 0; i < outputSize; i += 1) {
              for (size_t curDim = 0; curDim < gradDim && j + curDim < inputSize; curDim += 1) {
                jac(i, j + curDim) = GT::at(tape.getGradient(output[i]), curDim);
              }
            }

            setGradientOnIdentifier(tape, j, input, inputSize, typename GT::Real());
          }

          tape.clearAdjoints(end, start);

        } else if (EvaluationType::Reverse == evalType) {
          for (size_t i = 0; i < outputSize; i += gradDim) {
            setGradientOnIdentifier(tape, i, output, outputSize, typename GT::Real(1.0));

            if (keepState) {
              tape.evaluateKeepState(end, start);
            } else {
              tape.evaluate(end, start);
            }

            for (size_t j = 0; j < inputSize; j += 1) {
              for (size_t curDim = 0; curDim < gradDim && i + curDim < outputSize; curDim += 1) {
                jac(i + curDim, j) = GT::at(tape.getGradient(input[j]), curDim);
                GT::at(tape.gradient(input[j]), curDim) = typename GT::Real();
              }
            }

            setGradientOnIdentifier(tape, i, output, outputSize, typename GT::Real());

            if (!Config::ReversalZeroesAdjoints) {
              tape.clearAdjoints(end, start);
            }
          }
        } else {
          CODI_EXCEPTION("Evaluation mode not implemented. Mode is: %d.", (int)evalType);
        }
      }

      template<typename Jac>
      static CODI_INLINE void computeJacobian(Position const& start, Position const& end, Identifier const* input,
                                              size_t const inputSize, Identifier const* output, size_t const outputSize,
                                              Jac& jac) {
        computeJacobian(Type::getGlobalTape(), start, end, input, inputSize, output, outputSize, jac);
      }

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

      template<typename Hes, typename Jac = DummyJacobian>
      static CODI_INLINE void computeHessianPrimalValueTapeForward(Tape& tape, Position const& start,
                                                                   Position const& end, Identifier const* input,
                                                                   size_t const inputSize, Identifier const* output,
                                                                   size_t const outputSize, Hes& hes,
                                                                   Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        using GT1st = GT;
        size_t constexpr gradDim1st = GT1st::dim;
        using GT2nd = GradientTraits::TraitsImplementation<typename Real::Gradient>;
        size_t constexpr gradDim2nd = GT2nd::dim;

        // Assume tape that was just recorded
        tape.revertPrimals(start);

        for (size_t j = 0; j < inputSize; j += gradDim2nd) {
          setGradient2ndOnIdentifier(tape, j, input, inputSize, typename GT2nd::Real(1.0));

          // The k = j init is no problem, it will evaluated slightly more elements around the diagonal
          for (size_t k = j; k < inputSize; k += gradDim1st) {
            setGradientOnIdentifier(tape, k, input, inputSize, typename GT1st::Real(1.0));

            tape.evaluateForward(start, end);

            for (size_t i = 0; i < outputSize; i += 1) {
              for (size_t vecPos1st = 0; vecPos1st < gradDim1st && k + vecPos1st < inputSize; vecPos1st += 1) {
                for (size_t vecPos2nd = 0; vecPos2nd < gradDim2nd && j + vecPos2nd < inputSize; vecPos2nd += 1) {
                  hes(i, j + vecPos2nd, k + vecPos1st) =
                      GT2nd::at(GT1st::at(tape.getGradient(output[i]), vecPos1st).gradient(), vecPos2nd);
                  hes(i, k + vecPos1st, j + vecPos2nd) = hes(i, j + vecPos2nd, k + vecPos1st);  // symmetry
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

      template<typename Hes, typename Jac = DummyJacobian>
      static CODI_INLINE void computeHessianPrimalValueTapeReverse(Tape& tape, Position const& start,
                                                                   Position const& end, Identifier const* input,
                                                                   size_t const inputSize, Identifier const* output,
                                                                   size_t const outputSize, Hes& hes,
                                                                   Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        using GT1st = GT;
        size_t constexpr gradDim1st = GT1st::dim;
        using GT2nd = GradientTraits::TraitsImplementation<typename Real::Gradient>;
        size_t constexpr gradDim2nd = GT2nd::dim;

        // Assume tape that was just recorded
        tape.revertPrimals(start);

        for (size_t j = 0; j < inputSize; j += gradDim2nd) {
          setGradient2ndOnIdentifier(tape, j, input, inputSize, typename GT2nd::Real(1.0));

          // propagate the new derivative information
          tape.evaluatePrimal(start, end);

          for (size_t i = 0; i < outputSize; i += gradDim1st) {
            setGradientOnIdentifier(tape, i, output, outputSize, typename GT1st::Real(1.0));

            // propagate the derivatives backward for second order derivatives
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

      template<typename Func, typename VecIn, typename VecOut, typename Hes, typename Jac = DummyJacobian>
      static CODI_INLINE void computeHessianForward(Func func, VecIn& input, VecOut& output, Hes& hes,
                                                    Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        using GT1st = GT;
        size_t constexpr gradDim1st = GT1st::dim;
        using GT2nd = GradientTraits::TraitsImplementation<typename Real::Gradient>;
        size_t constexpr gradDim2nd = GT2nd::dim;

        Tape& tape = Type::getGlobalTape();

        for (size_t j = 0; j < input.size(); j += gradDim2nd) {
          setGradient2ndOnCoDiValue(j, input.data(), input.size(), typename GT2nd::Real(1.0));

          // propagate the new derivative information
          recordTape(func, input, output);

          // The k = j init is no problem, it will evaluated slightly more elements around the diagonal
          for (size_t k = j; k < input.size(); k += gradDim1st) {
            setGradientOnCoDiValue(tape, k, input.data(), input.size(), typename GT1st::Real(1.0));

            // propagate the derivatives forward for second order derivatives
            tape.evaluateForwardKeepState(tape.getZeroPosition(), tape.getPosition());

            for (size_t i = 0; i < output.size(); i += 1) {
              for (size_t vecPos1st = 0; vecPos1st < gradDim1st && k + vecPos1st < input.size(); vecPos1st += 1) {
                for (size_t vecPos2nd = 0; vecPos2nd < gradDim2nd && j + vecPos2nd < input.size(); vecPos2nd += 1) {
                  hes(i, j + vecPos2nd, k + vecPos1st) = GT2nd::at(
                      GT1st::at(tape.getGradient(output[i].getIdentifier()), vecPos1st).gradient(), vecPos2nd);
                  hes(i, k + vecPos1st, j + vecPos2nd) = hes(i, j + vecPos2nd, k + vecPos1st);  // symmetry
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

      template<typename Func, typename VecIn, typename VecOut, typename Hes, typename Jac = DummyJacobian>
      static CODI_INLINE void computeHessianReverse(Func func, VecIn& input, VecOut& output, Hes& hes,
                                                    Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        using GT1st = GT;
        size_t constexpr gradDim1st = GT1st::dim;
        using GT2nd = GradientTraits::TraitsImplementation<typename Real::Gradient>;
        size_t constexpr gradDim2nd = GT2nd::dim;

        Tape& tape = Type::getGlobalTape();

        for (size_t j = 0; j < input.size(); j += gradDim2nd) {
          setGradient2ndOnCoDiValue(j, input.data(), input.size(), typename GT2nd::Real(1.0));

          // propagate the new derivative information
          recordTape(func, input, output);

          for (size_t i = 0; i < output.size(); i += gradDim1st) {
            setGradientOnCoDiValue(tape, i, output.data(), output.size(), typename GT1st::Real(1.0));

            // propagate the derivatives backward for second order derivatives
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

      template<typename T>
      static CODI_INLINE void setGradientOnIdentifier(Tape& tape, size_t const pos, Identifier const* identifiers,
                                                      size_t const size, T value) {
        size_t constexpr gradDim = GT::dim;

        for (size_t curDim = 0; curDim < gradDim && pos + curDim < size; curDim += 1) {
          if (CODI_ENABLE_CHECK(ActiveChecks, 0 != identifiers[pos + curDim])) {
            GT::at(tape.gradient(identifiers[pos + curDim]), curDim) = value;
          }
        }
      }

      template<typename T>
      static CODI_INLINE void setGradient2ndOnIdentifier(Tape& tape, size_t const pos, Identifier const* identifiers,
                                                         size_t const size, T value) {
        using GT2nd = GradientTraits::TraitsImplementation<typename Real::Gradient>;
        size_t constexpr gradDim2nd = GT2nd::dim;

        for (size_t curDim = 0; curDim < gradDim2nd && pos + curDim < size; curDim += 1) {
          GT2nd::at(tape.primal(identifiers[pos + curDim]).gradient(), curDim) = value;
        }
      }

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

      template<typename T>
      static CODI_INLINE void setGradient2ndOnCoDiValue(size_t const pos, Type* identifiers, size_t const size,
                                                        T value) {
        using GT2nd = GradientTraits::TraitsImplementation<typename Real::Gradient>;
        size_t constexpr gradDim2nd = GT2nd::dim;

        for (size_t curDim = 0; curDim < gradDim2nd && pos + curDim < size; curDim += 1) {
          // No check required since this are forward types.
          GT2nd::at(identifiers[pos + curDim].value().gradient(), curDim) = value;
        }
      }

      template<typename Func, typename VecIn, typename VecOut>
      static CODI_INLINE void recordTape(Func func, VecIn& input, VecOut& output) {
        Tape& tape = Type::getGlobalTape();
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
