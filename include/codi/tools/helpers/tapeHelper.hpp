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

#include <vector>

#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../traits/realTraits.hpp"
#include "../../traits/tapeTraits.hpp"
#include "../algorithms.hpp"
#include "../data/hessian.hpp"
#include "../data/jacobian.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief A helper class that allows for a simpler access and management of a CoDiPack tape.
   *
   * This class provides helper functionality to record a tape and to evaluate the forward and reverse mode of AD as
   * well as capabilities to compute the full Jacobian and Hessian. Some functionality is only available with specific
   * CoDiPack types. For Hessian computations, a second order primal value type is required (e.g.
   * HessianComputationType) and for primal re-evaluations a primal value type is required (e.g.
   * RealReversePrimalIndex).
   *
   * The nomenclature and mathematical definitions for the function, the Jacobian, and the Hessian can be found in the
   * section \ref sec_namingConventions. For example, n denotes the number of inputs and m denotes the number of
   * outputs. Function arguments in this class follow the same naming scheme as in the referenced documentation.
   *
   * The general workflow for this class to record the representation of \f$ f \f$ is:
   * \snippet examples/Example_16_TapeHelper.cpp Tape recording
   *
   * The function `func` represents the implementation of \f$ f \f$ and before the function is called all inputs of that
   * function need to be registered on the tape. The same is true for all outputs after the function is called. The
   * methods startRecording() and stopRecording() define the region that will be recorded on the tape and are mandatory.
   * They prepare the helper and the tape structure of CoDiPack for the recording and the evaluation and ensure that
   * everything has the correct state.
   *
   * The order of registerInput() and registerOutput() calls is important. It defines which variable is represented
   * by the first entry, the second, etc. in the gradient vector, primal vector, Jacobian etc.
   *
   * The computation of derivatives can then be done with the functions evalForward(), evalReverse(), evalJacobian(),
   * and evalHessian(). For each of these functions, there is an `eval...At` equivalent that will first perform a primal
   * re-evaluation of the tape on the given inputs (only for primal value tapes) and then performs the reverse
   * evaluation.
   *
   * The function evalPrimal() is used in the `eval...At` methods and can be used to manually re-evaluate for a
   * different choice of input variables.
   *
   * All arguments for the methods can either be created with the corresponding create method or can be created
   * manually. The methods available for that are createGradientVectorInput(), createGradientVectorOutput(),
   * createJacobian(), createHessian(), createPrimalVectorInput() and createPrimalVectorOutput(). Each create function
   * has a corresponding delete function that deletes the objects.
   *
   * The computation of the Hessian could be performed as follows.
   * \snippet examples/Example_16_TapeHelper.cpp Hessian evaluation
   *
   * A simple reverse evaluation works like this.
   * \snippet examples/Example_16_TapeHelper.cpp Reverse evaluation
   *
   * The tape helper can be used to record multiple different tapes. Each time startRecording() is called, the old
   * recording is deleted and a new one is started.
   *
   * For a more detailed example see \ref Example_16_TapeHelper.
   *
   * @tparam T_Type  The CoDiPack type on which the evaluations take place.
   * @tparam T_Impl  The type of the implementing class for the virtual template methods.
   */
  template<typename T_Type, typename T_Impl>
  struct TapeHelperBase {
    public:

      /// See TapeHelperBase.
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);
      using Impl = CODI_DD(T_Impl, TapeHelperBase);  ///< See TapeHelperBase.

      using Real = typename Type::Real;              ///< See LhsExpressionInterface.
      using Identifier = typename Type::Identifier;  ///< See LhsExpressionInterface.
      using Gradient = typename Type::Gradient;      ///< See LhsExpressionInterface.

      using PassiveReal = typename RealTraits::PassiveReal<Real>;  ///< Passive base of the CoDiPack type.

      using JacobianType = Jacobian<PassiveReal>;  ///< Type of the Jacobian.
      using HessianType = Hessian<PassiveReal>;    ///< Type of the Hessian.

    protected:

      using Tape = typename Type::Tape;  ///< Underlying tape type.

      Tape& tape;  ///< Reference to the global tape.

      std::vector<Identifier> inputValues;   ///< Input value identifiers.
      std::vector<Identifier> outputValues;  ///< Input value identifiers.

      bool wasForwardEvaluated;  ///< State of the last evaluation.

    public:

      /// Constructor
      TapeHelperBase() : tape(Type::getTape()), inputValues(), outputValues(), wasForwardEvaluated(false) {}

      /// Destructor
      virtual ~TapeHelperBase() {}

      /**
       * @brief Create a gradient vector that can hold the tangent/adjoint of the input variables.
       *
       * Should only be called after the tape has been recorded.
       * Needs to be deleted with deleteGradientVector.
       *
       * @return Vector of size n.
       */
      Gradient* createGradientVectorInput() {
        return new Gradient[getInputSize()];
      }

      /**
       * @brief Create a gradient vector that can hold the tangent/adjoint of the output variables.
       *
       * Should only be called after the tape has been recorded.
       * Needs to be deleted with deleteGradientVector.
       *
       * @return Vector of size m.
       */
      Gradient* createGradientVectorOutput() {
        return new Gradient[getOutputSize()];
      }

      /**
       * @brief Create a Jacobian that can hold the Jacobian of the recorded tape.
       *
       * Should only be called after the tape has been recorded.
       * Needs to be deleted with deleteJacobian.
       *
       * @return A Jacobian with the size m,n.
       */
      JacobianType& createJacobian() {
        JacobianType* jacPointer = new JacobianType(getOutputSize(), getInputSize());

        return *jacPointer;
      }

      /**
       * @brief Create a Hessian that can hold the Hessian of the recorded tape.
       *
       * Should only be called after the tape has been recorded.
       * Needs to be deleted with deleteHessian.
       *
       * @return A Hessian with the size m,n.
       */
      HessianType& createHessian() {
        HessianType* hesPointer = new HessianType(getOutputSize(), getInputSize());

        return *hesPointer;
      }

      /**
       * @brief Create a primal vector that can hold the primal seeding of the input variables.
       *
       * Should only be called after the tape has been recorded.
       * Needs to be deleted with deletePrimalVector.
       *
       * @return Vector of size n.
       */
      Real* createPrimalVectorInput() {
        return new Real[getInputSize()];
      }

      /**
       * @brief Create a primal vector that can hold the primal result of the output variables.
       *
       * Should only be called after the tape has been recorded.
       * Needs to be deleted with deletePrimalVector.
       *
       * @return Vector of size m.
       */
      Real* createPrimalVectorOutput() {
        return new Real[getOutputSize()];
      }

      /// Delete a gradient vector that was created with createGradientVectorInput or createGradientVectorOutput.
      void deleteGradientVector(Gradient* vec) {
        delete[] vec;

        vec = nullptr;
      }

      ///  Delete the Jacobian that was created with createJacobian function.
      void deleteJacobian(JacobianType& jac) {
        JacobianType* jacPointer = &jac;

        delete jacPointer;
      }

      /// Delete the Hessian that was created with createHessian function.
      void deleteHessian(HessianType& hes) {
        HessianType* hesPointer = &hes;

        delete hesPointer;
      }

      /// Delete a primal vector that was created with createPrimalVectorInput or createPrimalVectorOutput.
      void deletePrimalVector(Real* vec) {
        delete[] vec;

        vec = nullptr;
      }

      /// Get the number of registered inputs. Call after stopRecording().
      /// @return n
      size_t getInputSize() {
        return inputValues.size();
      }

      /// Get the number of registered outputs. Call after stopRecording().
      /// @return m
      size_t getOutputSize() {
        return outputValues.size();
      }

      /**
       * @brief Add an input variable to the tape.
       *
       * Input variables are also known as independent variables.
       *
       * The value is modified such that CoDiPack will recognize it as an active variable. For all variables registered
       * with the method, the derivative is computed.
       *
       * With this tape helper, the sequence of the registerInput calls is important. All primal and derivative vectors
       * will use the first entry for the first value registered, the second for the second one, and so on.
       */
      void registerInput(Type& value) {
        tape.registerInput(value);
        inputValues.push_back(value.getIdentifier());
      }

      /**
       * @brief Add an output variable to the tape.
       *
       * Output variables are also known as dependent variables.
       *
       * The value is modified such that it is safe to use the variable as an output for which the reverse seed can be
       * set.
       *
       * With this tape helper, the sequence of the registerOutput calls is important. All primal and derivative vectors
       * will use the first entry for the first value registered, the second for the second one, and so on.
       */
      void registerOutput(Type& value) {
        tape.registerOutput(value);
        outputValues.push_back(value.getIdentifier());
      }

      /// Start the recording process. Deletes the old tape.
      void startRecording() {
        tape.reset();
        inputValues.clear();
        outputValues.clear();

        tape.setActive();
      }

      /// Stop the recording process.
      void stopRecording() {
        tape.setPassive();
      }

      /**
       * @brief Perform a primal re-evaluation of the tape.
       *
       * The re-evaluation will change the internally stored primal variables of the tape.
       *
       * @param[in] x  The new seeding vector for the primal input variables. The sequence of variables is the same as
       *               for the register input call. The vector should be created with createPrimalVectorInput.
       * @param[out] y  The result of the primal evaluation. The sequence of variables is the same as for the register
       *               output call. If the pointer is a null pointer then the result is not stored. The vector should
       *               be created with createPrimalVectorOutput.
       */
      virtual void evalPrimal(Real const* x, Real* y = nullptr) = 0;

      /**
       * @brief Perform a forward (tangent) evaluation of the recorded tape.
       * @param[in]  x_d  The seeding vector for the input variables (independent variables).
       *                  The vector should be created with createGradientVectorInput.
       * @param[out] y_d  The result vector for the output variables (dependent variables).
       *                  The vector should be created with createGradientVectorOutput.
       */
      CODI_INLINE void evalForward(Gradient const* x_d, Gradient* y_d) {
        changeStateToForwardEvaluation();

        for (size_t j = 0; j < inputValues.size(); j += 1) {
          tape.setGradient(inputValues[j], x_d[j]);
        }

        tape.evaluateForward();

        for (size_t i = 0; i < outputValues.size(); i += 1) {
          y_d[i] = tape.getGradient(outputValues[i]);
          tape.setGradient(outputValues[i], Gradient());
        }
      }

      /**
       * @brief Re-evaluate the tape with new input variables and compute the AD forward mode.
       *
       * This method is a shortcut for calling evalPrimal and evalForward in succession.
       *
       * @param[in]    x  The new seeding vector for the primal input variables. The sequence of variables is the same
       *                  as for the register input call. The vector should be created with createPrimalVectorInput.
       * @param[in]  x_d  The seeding vector for the input variables (independent variables).
       *                  The vector should be created with createGradientVectorInput.
       * @param[out] y_d  The result vector for the output variables (dependent variables).
       *                  The vector should be created with createGradientVectorOutput.
       * @param[out]   y  The result of the primal evaluation. The sequence of variables is the same as for the register
       *                  output call. If the pointer is a null pointer then the result is not stored. The vector should
       *                  be created with createPrimalVectorOutput.
       */
      CODI_INLINE void evalForwardAt(Real const* x, Gradient const* x_d, Gradient* y_d, Real* y = nullptr) {
        evalPrimal(x, y);

        evalForward(x_d, y_d);
      }

      /**
       * @brief Perform a reverse (adjoint) evaluation of the recorded tape.
       * @param[in]  y_b  The seeding vector for the output variables (dependent variables).
       *                  The vector should be created with createGradientVectorOutput.
       * @param[out] x_b  The result vector for the input variables (independent variables).
       *                  The vector should be created with createGradientVectorInput.
       */
      CODI_INLINE void evalReverse(Gradient const* y_b, Gradient* x_b) {
        changeStateToReverseEvaluation();

        for (size_t i = 0; i < outputValues.size(); i += 1) {
          tape.setGradient(outputValues[i], y_b[i]);
        }

        tape.evaluate();

        for (size_t j = 0; j < inputValues.size(); j += 1) {
          x_b[j] = tape.getGradient(inputValues[j]);
          tape.setGradient(inputValues[j], Gradient());
        }

        if (!Config::ReversalZeroesAdjoints) {
          tape.clearAdjoints();
        }
      }

      /**
       * @brief Re-evaluate the tape with new input variables and compute the AD forward mode.
       *
       * This method is a shortcut for calling evalPrimal and evalReverse in succession.
       *
       * @param[in]    x  The new seeding vector for the primal input variables. The sequence of variables is the same
       *                  as for the register input call. The vector should be created with createPrimalVectorInput.
       * @param[out] y_b  The seeding vector for the output variables (dependent variables).
       *                  The vector should be created with createGradientVectorOutput.
       * @param[in]  x_b  The result vector for the input variables (independent variables).
       *                  The vector should be created with createGradientVectorInput.
       * @param[out]   y  The result of the primal evaluation. The sequence of variables is the same as for the register
       *                  output call. If the pointer is a null pointer then the result is not stored. The vector should
       *                  be created with createPrimalVectorOutput.
       */
      CODI_INLINE void evalReverseAt(Real const* x, Gradient const* y_b, Gradient* x_b, Real* y = nullptr) {
        evalPrimal(x, y);

        evalReverse(y_b, x_b);
      }

      /**
       * @brief Evaluates the full Jacobian of the recorded tape.
       *
       * The algorithm will select the best choice for the evaluation, either a forward mode or reverse mode evaluation.
       * It will also use the vector mode if the underlying tape was configured with such a mode.
       *
       * @param[out] jac  The storage for the Jacobian which is evaluated. Must have the correct size and should be
       *                  created with createJacobian.
       */
      CODI_INLINE void evalJacobian(JacobianType& jac) {
        JacobianConvertWrapper<JacobianType> wrapper(jac);

        evalJacobianGen(wrapper);
      }

      /**
       * @brief Re-evaluate the tape with new input variables and compute the full Jacobian at the new inputs.
       *
       * This method is a shortcut for calling evalPrimal and evalJacobian in succession.
       *
       * @param[in]    x  The new seeding vector for the primal input variables. The sequence of variables is the same
       *                  as for the register input call. The vector should be created with createPrimalVectorInput.
       * @param[out] jac  The storage for the Jacobian which is evaluated. Must have the correct size and should be
       *                  created with createJacobian.
       * @param[out]   y  The result of the primal evaluation. The sequence of variables is the same as for the register
       *                  output call. If the pointer is a null pointer then the result is not stored. The vector should
       *                  be created with createPrimalVectorOutput.
       */
      CODI_INLINE void evalJacobianAt(Real const* x, JacobianType& jac, Real* y = nullptr) {
        evalPrimal(x, y);

        evalJacobian(jac);
      }

      /**
       * @brief Evaluates the full Jacobian of the recorded tape with a custom Jacobian type chosen by the user.
       *
       * \copydetails evalJacobian
       *
       * @tparam Jac  Has to implement JacobianInterface.
       */
      template<typename Jac>
      CODI_INLINE void evalJacobianGen(Jac& jac) {
        using Algo = Algorithms<Type>;
        typename Algo::EvaluationType evalType = Algo::getEvaluationChoice(inputValues.size(), outputValues.size());

        if (Algo::EvaluationType::Forward == evalType) {
          changeStateToForwardEvaluation();
        } else if (Algo::EvaluationType::Reverse == evalType) {
          changeStateToReverseEvaluation();
        } else {
          CODI_EXCEPTION("Evaluation type not implemented.");
        }

        Algorithms<Type>::template computeJacobian<Jac, false>(tape, tape.getZeroPosition(), tape.getPosition(),
                                                               inputValues.data(), inputValues.size(),
                                                               outputValues.data(), outputValues.size(), jac);
      }

      /**
       * @brief Evaluates the full Hessian of the recorded tape.
       *
       * The algorithm will select the best choice for the evaluation, either a forward mode or reverse mode evaluation.
       * It will also use the vector mode if the underlying tape was configured with such a mode.
       *
       * @param[out] hes  The storage for the Hessian which is evaluated. Must have the correct size and should be
       *                  created with createHessian.
       * @param[out] jac  If also the Jacobian should be computed alongside the Hessian, a storage for the Jacobian can
       *                  be provided. Must have the correct size and should be created with createJacobian.
       *
       * @tparam Jac  Has to implement JacobianInterface.
       */
      template<typename Jac = DummyJacobian>
      void evalHessian(HessianType& hes, Jac& jac = StaticDummy<DummyJacobian>::dummy);

      /**
       * @brief Re-evaluate the tape with new input variables and compute the full Hessian at the new inputs.
       *
       * This method is a shortcut for calling evalPrimal and evalHessian
       *
       * @param[in]    x  The new seeding vector for the primal input variables. The sequence of variables is the same
       *                  as for the register input call. The vector should be created with createPrimalVectorInput.
       * @param[out] hes  The storage for the Hessian which is evaluated. Needs to have the correct size and should
       *                  be created with createHessian.
       * @param[out]   y  The result of the primal evaluation. The sequence of variables is the same as for the register
       *                  output call. If the pointer is a null pointer then the result is not stored. The vector should
       *                  be created with createPrimalVectorOutput.
       * @param[out] jac  If also the Jacobian should be computed alongside the Hessian, a storage for the Jacobian can
       *                  be provided. Needs to have the correct size and should be created with createJacobian.
       *
       * @tparam Jac  Has to implement JacobianInterface.
       */
      template<typename Jac = DummyJacobian>
      CODI_INLINE void evalHessianAt(Real const* x, HessianType& hes, Real* y = nullptr,
                                     Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        evalPrimal(x, y);

        cast().evalHessian(hes, jac);
      }

    protected:

      /// Cast to the implementing class.
      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

      /// Change state.
      void changeStateToForwardEvaluation() {
        wasForwardEvaluated = true;

        // No cleanup to do.
      }

      /// Change state and clear the adjoints.
      void changeStateToReverseEvaluation() {
        if (wasForwardEvaluated) {
          // Forward evaluation leaves the adjoint vector dirty.

          tape.clearAdjoints();
        }

        wasForwardEvaluated = false;
      }
  };

  // clang-format off
  /// TapeHelper for a CoDiPack type that currently do not support a TapeHelper implementation. Will generate errors.
  ///
  /// See TapeHelperBase for details.
  ///
  /// @tparam T_Type  A CoDiPack type that does not support TapeHelper.
  template<typename T_Type>
  struct TapeHelperNoImpl : public TapeHelperBase<T_Type, TapeHelperNoImpl<T_Type>> {
    public:

      CODI_STATIC_ASSERT(false && std::is_void<T_Type>::value, "Tape helper not implemented for this tape.");

      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See TapeHelperBase.
      using Real = typename Type::Real;                           ///< See TapeHelperBase.

      using Base = TapeHelperBase<Type, TapeHelperNoImpl<Type>>;  ///< Base class abbreviation.

      /// Missing implementation will yield linker errors.
      virtual void evalPrimal(Real const* x, Real* y = nullptr) CODI_DD(= 0;, {})

      /// Missing implementation will yield linker errors.
      template<typename Jac = DummyJacobian>
      void evalHessian(typename Base::HessianType& hes, Jac& jac = StaticDummy<DummyJacobian>::dummy);
  };
  // clang-format on

  /// TapeHelper implementation for the Jacobian taping strategy.
  ///
  /// See TapeHelperBase for details.
  ///
  /// @tparam T_Type  The CoDiPack type on which the evaluations take place.
  template<typename T_Type>
  struct TapeHelperJacobi : public TapeHelperBase<T_Type, TapeHelperJacobi<T_Type>> {
    public:

      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See TapeHelperBase.
      using Real = typename Type::Real;                           ///< See TapeHelperBase.

      using Base = TapeHelperBase<Type, TapeHelperJacobi<Type>>;  ///< Base class abbreviation.

      /// Throws an exception since primal evaluations are not support by Jacobian tapes.
      virtual void evalPrimal(Real const* x, Real* y = nullptr) {
        CODI_UNUSED(x, y);

        CODI_EXCEPTION(
            "No primal evaluation for Jacobian tapes. "
            "Please use codi::RealReversePrimal or codi::RealReversePrimalIndex types for this kind of functionality.");
      }

      /// Throws an exception since primal evaluations are not supported by Jacobian tapes.
      template<typename Jac = DummyJacobian>
      void evalHessian(typename Base::HessianType& hes, Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        CODI_UNUSED(hes, jac);

        CODI_EXCEPTION(
            "No direct Hessian evaluation for Jacobian tapes. "
            "Please use codi::RealReversePrimal or codi::RealReversePrimalIndex types for this kind of functionality "
            "or the EvaluationHelper class.");
      }
  };

  /// TapeHelper implementation for the Jacobian taping strategy.
  ///
  /// See TapeHelperBase for details.
  ///
  /// @tparam T_Type  The CoDiPack type on which the evaluations take place.
  template<typename T_Type>
  struct TapeHelperPrimal : public TapeHelperBase<T_Type, TapeHelperPrimal<T_Type>> {
    public:

      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See TapeHelperBase.
      using Real = typename Type::Real;                           ///< See TapeHelperBase.

      using Base = TapeHelperBase<Type, TapeHelperPrimal<Type>>;  ///< Base class abbreviation.

      /// \copydoc TapeHelperBase::evalPrimal
      virtual void evalPrimal(Real const* x, Real* y = nullptr) {
        for (size_t j = 0; j < this->inputValues.size(); j += 1) {
          this->tape.primal(this->inputValues[j]) = x[j];
        }

        this->tape.evaluatePrimal();

        if (nullptr != y) {
          for (size_t i = 0; i < this->outputValues.size(); i += 1) {
            y[i] = this->tape.primal(this->outputValues[i]);
          }
        }
      }

      /// \copydoc TapeHelperBase::evalHessian
      template<typename Jac = DummyJacobian>
      void evalHessian(typename Base::HessianType& hes, Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        using Algo = Algorithms<Type>;
        typename Algo::EvaluationType evalType =
            Algo::getEvaluationChoice(this->inputValues.size(), this->outputValues.size());

        if (Algo::EvaluationType::Forward == evalType) {
          this->changeStateToForwardEvaluation();
        } else if (Algo::EvaluationType::Reverse == evalType) {
          this->changeStateToReverseEvaluation();
        } else {
          CODI_EXCEPTION("Evaluation type not implemented.");
        }

        Algorithms<Type>::computeHessianPrimalValueTape(
            this->tape, this->tape.getZeroPosition(), this->tape.getPosition(), this->inputValues.data(),
            this->inputValues.size(), this->outputValues.data(), this->outputValues.size(), hes, jac);
      }
  };

  /// See TapeHelperBase.
  template<typename Type, typename = void>
  struct TapeHelper : public TapeHelperNoImpl<Type> {};

#ifndef DOXYGEN_DISABLE
  /// See TapeHelperBase.
  template<typename Type>
  struct TapeHelper<Type, TapeTraits::EnableIfJacobianTape<typename Type::Tape>> : public TapeHelperJacobi<Type> {};

  /// See TapeHelperBase.
  template<typename Type>
  struct TapeHelper<Type, TapeTraits::EnableIfPrimalValueTape<typename Type::Tape>> : public TapeHelperPrimal<Type> {};
#endif

}
