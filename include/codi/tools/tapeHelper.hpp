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

#include "evaluationHelper.hpp"
#include "../configure.h"
#include "../tapes/tapeTraits.hpp"

#include <vector>

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief A helper class that allows for a much simpler access and management of a CoDiPack tape.
   *
   *
   * This class provides helper functionality to record a tape and to evaluate the forward and reverse mode of AD as
   * well as the capability to compute the full Jacobian and Hessian. Some functionality is only available with specific
   * CoDiPack types. For Hessian computations a second order primal value type is required (e.g. HessianComputationType)
   * and for primal reevaluations a primal value type is required (e.g. RealReversePrimalIndex).
   *
   * The nomenclature and mathematical definitions for the function, the Jacobian and the Hessian can be found in the
   * \ref FunctionDefinition "Algorithms" documentation. Function arguments in this class follow the same naming scheme
   * as in the referenced documentation.
   *
   * The general workflow for this class to record the representation of \f$ f \f$ is:
   *
   * \code{.cpp}
   *  TapeHelper<HessianComputationType> th;
   *
   *  th.startRecording();
   *
   *  //for each input
   *  th.registerInput(input);
   *
   *  func();
   *
   *  //for each output
   *  th.registerOutput(output);
   *
   *  th.stopRecording();
   * \endcode
   *
   * The function func represents the implementation of \f$ f \f$ and before the function is called all inputs of that
   * function need to be registered to the tape. The same is true for all outputs after the function is called. The
   * methods startRecording() and stopRecording() define the region that will be recorded on the tape and are mandatory.
   * They prepare the helper and the tape structure of CoDiPack for the recording and the evaluation and ensure that
   * everything has the correct state.
   *
   * The order of registerInput() and registerOutput() calls is important. It defines which variable is represented
   * by the first entry, the second, etc. in the gradient vector, primal vector, Jacobian etc.
   *
   * The computation of derivatives can then be done with the functions evalForward(), evalReverse(), evalJacobian(),
   * and evalHessian(). For each of these functions an eval...At method is also available, that will first reevaluate the
   * tape at the given position (only for primal value tapes) and then perform the evaluation.
   *
   * The function evalPrimal() is used in the eval...At methods and can be used to manually reevaluate the tape at the
   * given position.
   *
   * All arguments for the method can either be created with the corresponding create method or can be created manually.
   * The methods available for that are createGradientVectorInput(), createGradientVectorOutput(), createJacobian(),
   * createHessian(), createPrimalVectorInput() and createPrimalVectorOutput(). Each create function has a corresponding delete
   * function that deletes the objects.
   *
   * The computation of the Hessian could then be done as:
   * \code{.cpp}
   *  auto& hes = th.createHessian();
   *
   *  th.evalHessian(hes);
   *
   *  // get results ... = hes(i,j,k);
   *
   *  th.deleteHessian(hes);
   * \endcode
   *
   * A simple reverse evaluation is done as:
   * \code{.cpp}
   *  auto x_d = th.createGradientVectorInput();
   *  auto y_d = th.createGradientVectorOutput();
   *
   *  // Set the seeding
   *  for(size_t i = 0; i < th.getOutputSize(); i += 1) {
   *    y_d[i] = ...;
   *  }
   *
   *  th.evalReverse(y_d, x_d);
   *
   *  // Use the adjoint
   *  for(size_t i = 0; i < th.getInputSize(); i += 1) {
   *    ... = x_d[i];
   *  }
   *
   *  th.deleteGradientVector(x_d);
   *  th.deleteGradientVector(y_d);
   * \endcode
   *
   * The tape helper can be used to record multiple different tapes. Each time startRecording() is called the old
   * recording is deleted and a new one is started.
   *
   * @tparam CoDiType  A CoDiPack type e.g. codi::HessianComputationType, codi::JacobianComputationType, codi::RealReverse
   * @tparam     Impl  The type of the implementing class for the virtual template methods.
   */
  template<typename CoDiType, typename Impl>
  class TapeHelperBase {
    public:
      using Real = typename CoDiType::Real; /**< The floating point calculation type in the CoDiPack types. */
      using GradientData = typename CoDiType::GradientData; /**< The type for the identification of gradients. */
      using GradientValue = typename CoDiType::GradientValue; /**< The type for the gradient computation */

      using PassiveReal = typename TypeTraits<Real>::PassiveReal; /**< The basic type of the CoDiPack type */

      using JacobianType = Jacobian<std::vector<PassiveReal>>; /**< The Jacobian type for the evaluation. */
      using HessianType = Hessian<std::vector<PassiveReal>>; /**< The Hessian type for the evaluation. */

    protected:
      typedef typename CoDiType::TapeType Tape; /**< The type of the tape implementation. */

      Tape& tape; /**< Reference to the global tape .*/

      std::vector<GradientData> inputValues; /**< Storage for the identifiers of the input values. */
      std::vector<GradientData> outputValues; /**< Storage for the identifiers of the output values. */

      bool wasForwardEvaluated; /**< Tape state of the last evaluation */

    public:

      /**
       * @brief Constructor for the helper.
       */
      TapeHelperBase() :
        tape(CoDiType::getGlobalTape()),
        inputValues(),
        outputValues(),
        wasForwardEvaluated(false)
      {}

      /**
       * @brief Destructor for the helper.
       */
      virtual ~TapeHelperBase() {}

      /**
       * @brief Create a gradient vector that can hold the tangent/adjoint of the input variables.
       *
       * * Should only be called after the tape has been recorded.
       * Needs to be deleted with deleteGradientVector.
       *
       * @return Vector of size n
       */
      GradientValue* createGradientVectorInput() {
        return createGradientVector(getInputSize());
      }

      /**
       * @brief Create a gradient vector that can hold the tangent/adjoint of the output variables.
       *
       * Should only be called after the tape has been recorded.
       * Needs to be deleted with deleteGradientVector.
       *
       * @return Vector of size m
       */
      GradientValue* createGradientVectorOutput() {
        return createGradientVector(getOutputSize());
      }

      /**
       * @brief Create a Jacobian that can hold the Jacobian of the recorded tape.
       *
       * Should only be called after the tape has been recorded.
       * Needs to be deleted with deleteJacobian.
       *
       * @return A Jacobian with the size m,n
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
       * @return A Hessian with the size m,n
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
       * @return Vector of size n
       */
      Real* createPrimalVectorInput() {
        return createPrimalVector(getInputSize());
      }

      /**
       * @brief Create a primal vector that can hold the primal result of the output variables.
       *
       * Should only be called after the tape has been recorded.
       * Needs to be deleted with deletePrimalVector.
       *
       * @return Vector of size m
       */
      Real* createPrimalVectorOutput() {
        return createPrimalVector(getOutputSize());
      }

      /**
       * @brief Delete a gradient vector that was created with createGradientVectorInput or createGradientVectorOutput.
       *
       * @param[in,out] vec  Pointer is deleted and set to null.
       */
      void deleteGradientVector(GradientValue* vec) {
        delete [] vec;

        vec = nullptr;
      }

      /**
       * @brief Delete the Jacobian that was created with createJacobian function.
       *
       * @param[in,out] jac  The reference to the Jacobian.
       */
      void deleteJacobian(JacobianType& jac) {
        JacobianType* jacPointer = &jac;

        delete jacPointer;
      }

      /**
       * @brief Delete the Hessian that was created with createHessian function.
       *
       * @param[in,out] hes  The reference to the Hessian.
       */
      void deleteHessian(HessianType& hes) {
        HessianType* hesPointer = &hes;

        delete hesPointer;
      }

      /**
       * @brief Delete a primal vector that was created with createPrimalVectorInput or createPrimalVectorOutput.
       *
       * @param[in,out] vec  Pointer is deleted and set to null.
       */
      void deletePrimalVector(Real* vec) {
        delete [] vec;

        vec = nullptr;
      }

      /**
       * @brief Get the number of inputs from the tape.
       *
       * @return n
       */
      size_t getInputSize() {
        return inputValues.size();
      }

      /**
       * @brief Get the number of outputs from the tape.
       *
       * @return m
       */
      size_t getOutputSize() {
        return outputValues.size();
      }

      /**
       * @brief Add an input variable to the tape.
       *
       * Also know as and independent variable in the context of AD or ADOL-c
       *
       * The value is modified such that CoDiPack will recognize it as an active variable. For all variables registered
       * with the method the derivative is computed.
       *
       * With this tape helper the sequence of the registerInput calls is important. All primal and derivative vectors
       * will use the first entry for the first value registered. The second for the second one etc.
       *
       * @param[in,out] value  The value which is registered as an input.
       */
      void registerInput(CoDiType& value) {
        tape.registerInput(value);
        inputValues.push_back(value.getGradientData());
      }

      /**
       * @brief Add an output variable to the tape.
       *
       * Also know as and dependent variable in the context of AD or ADOL-c
       *
       * The value is modified such that it is safe to use the variable as an output for which the reverse seed can be
       * set.
       *
       * With this tape helper the sequence of the registerOutput calls is important. All primal and derivative vectors
       * will use the first entry for the first value registered. The second for the second one etc.
       *
       * @param[in,out] value  The value which is registered as an output.
       */
      void registerOutput(CoDiType& value) {
        tape.registerOutput(value);
        outputValues.push_back(value.getGradientData());
      }

      /**
       * @brief Start recording of a new tape.
       *
       * Resets all variables of the TapeHelper object for a new tape recording.
       *
       * CoDiPack records only the statements in the program with active variables. If statements or other branching
       * statements are ignored and only the active branch is recorded.
       */
      void startRecording() {
        tape.reset();
        inputValues.clear();
        outputValues.clear();

        tape.setActive();
      }

      /**
       * @brief Stops the recording of a tape.
       *
       * The tape is now finished and no further recordings can take place.
       * A new recording can be started with startRecording.
       */
      void stopRecording() {
        tape.setPassive();
      }

      /**
       * @brief Perform a primal reevaluation of the tape.
       *
       * The reevaluation will change the internally stored primal variables of the tape.
       *
       * @param[in] x  The new seeding vector for the primal input variables. The sequence of variables is the same as
       *               for the register input call. The vector should be created with createPrimalVectorInput.
       * @param[out] y  The result of the primal evaluation. The sequence of variables is the same as for the register
       *               output call. If the pointer is a null pointer then the result is not stored.  The vector should
       *               be created with createPrimalVectorOutput.
       */
      virtual void evalPrimal(Real const* x, Real* y = nullptr) = 0;

      /**
       * @brief Perform a forward (tangent) evaluation of the recorded tape.
       * @param[in]  x_d  The seeding vector for the input variables (independent variables)
       *                  The vector should be created with createGradientVectorInput.
       * @param[out] y_d  The result vector for the output variables (dependent variables).
       *                  The vector should be created with createGradientVectorOutput.
       */
      CODI_INLINE void evalForward(GradientValue const* x_d, GradientValue* y_d) {
        changeStateToForwardEvaluation();

        for(size_t j = 0; j < inputValues.size(); j += 1) {
          tape.setGradient(inputValues[j], x_d[j]);
        }

        tape.evaluateForward();

        for(size_t i = 0; i < outputValues.size(); i += 1) {
          y_d[i] = tape.getGradient(outputValues[i]);
          tape.setGradient(outputValues[i], GradientValue());
        }
      }

      /**
       * @brief Reevaluate the tape at the given position an compute the AD forward mode.
       *
       * This method is a shortcut for calling evalPrimal and evalForward in succession.
       *
       * @param[in]    x  The new seeding vector for the primal input variables. The sequence of variables is the same as
       *                  for the register input call. The vector should be created with createPrimalVectorInput.
       * @param[in]  x_d  The seeding vector for the input variables (independent variables)
       *                  The vector should be created with createGradientVectorInput.
       * @param[out] y_d  The result vector for the output variables (dependent variables)
       *                  The vector should be created with createGradientVectorOutput.
       * @param[out]   y  The result of the primal evaluation. The sequence of variables is the same as for the register
       *                  output call. If the pointer is a null pointer then the result is not stored. The vector should
       *                  be created with createPrimalVectorOutput.
       */
      CODI_INLINE void evalForwardAt(Real const* x, GradientValue const* x_d, GradientValue* y_d, Real* y = nullptr) {
        evalPrimal(x, y);

        evalForward(x_d, y_d);
      }

      /**
       * @brief Perform a reverse (adjoint) evaluation of the recorded tape.
       * @param[in]  y_b  The seeding vector for the output variables (dependent variables).
       *                  The vector should be created with createGradientVectorOutput.
       * @param[out] x_b  The result vector for the input variables (independent variables)
       *                  The vector should be created with createGradientVectorInput.
       */
      CODI_INLINE void evalReverse(GradientValue const* y_b, GradientValue* x_b) {
        changeStateToReverseEvaluation();

        for(size_t i = 0; i < outputValues.size(); i += 1) {
          tape.setGradient(outputValues[i], y_b[i]);
        }

        tape.evaluate();


        for(size_t j = 0; j < inputValues.size(); j += 1) {
          x_b[j] = tape.getGradient(inputValues[j]);
          tape.setGradient(inputValues[j], GradientValue());
        }

        if(!ZeroAdjointReverse) {
          tape.clearAdjoints();
        }
      }

      /**
       * @brief Reevaluate the tape at the given position an compute the AD forward mode.
       *
       * This method is a shortcut for calling evalPrimal and evalReverse in succession.
       *
       * @param[in]    x  The new seeding vector for the primal input variables. The sequence of variables is the same as
       *                  for the register input call. The vector should be created with createPrimalVectorInput.
       * @param[out] y_b  The seeding vector for the output variables (dependent variables)
       *                  The vector should be created with createGradientVectorOutput.
       * @param[in]  x_b  The result vector for the input variables (independent variables)
       *                  The vector should be created with createGradientVectorInput.
       * @param[out]   y  The result of the primal evaluation. The sequence of variables is the same as for the register
       *                  output call. If the pointer is a null pointer then the result is not stored. The vector should
       *                  be created with createPrimalVectorOutput.
       */
      CODI_INLINE void evalReverseAt(Real const* x, GradientValue const* y_b, GradientValue* x_b, Real* y = nullptr) {
        evalPrimal(x, y);

        evalReverse(y_b, x_b);
      }

      /**
       * @brief Evaluates the full Jacobian of the recorded tape.
       *
       * The algorithm will select the best choice for the evaluation. Either a forward mode or reverse mode evaluation.
       * It will also use the vector mode if the underlying tape was configured with such a mode.
       *
       * @param[out] jac  The storage for the Jacobian which is evaluated. Needs to have the correct size and should
       *                  be created with createJacobian.
       */
      CODI_INLINE void evalJacobian(JacobianType& jac) {
        JacobianConvertWrapper<JacobianType> wrapper(jac);

        evalJacobianGen(wrapper);
      }

      /**
       * @brief Reevaluate the tape at a given position and compute the full Jacobian at the new position.
       *
       * This method is a shortcut for calling evalPrimal and evalJacobian
       *
       * @param[in]    x  The new seeding vector for the primal input variables. The sequence of variables is the same as
       *                  for the register input call. The vector should be created with createPrimalVectorInput.
       * @param[out] jac  The storage for the Jacobian which is evaluated. Needs to have the correct size and should
       *                  be created with createJacobian.
       * @param[out]   y  The result of the primal evaluation. The sequence of variables is the same as for the register
       *                  output call. If the pointer is a null pointer then the result is not stored. The vector should
       *                  be created with createPrimalVectorOutput.
       */
      CODI_INLINE void evalJacobianAt(Real const* x, JacobianType& jac, Real* y = nullptr) {
        evalPrimal(x, y);

        evalJacobian(jac);
      }

      /**
       * @brief Evaluates the full Jacobian of the recorded tape with a custom Jacobian type from the user.
       *
       * \copydetails evalJacobian
       *
       * @tparam Jac  Needs to implement the access operators (e.g. operator()) as in the Jacobian class.
       */
      template<typename Jac>
      CODI_INLINE void evalJacobianGen(Jac& jac) {

        using Algo = Algorithms<CoDiType>;
        typename Algo::EvaluationType evalType = Algo::getEvaluationChoice(inputValues.size(), outputValues.size());

        if(Algo::EvaluationType::Forward == evalType) {
          changeStateToForwardEvaluation();
        } else if(Algo::EvaluationType::Reverse == evalType) {
          changeStateToReverseEvaluation();
        } else {
          CODI_EXCEPTION("Evaluation type not implemented.");
        }

        Algorithms<CoDiType>::template computeJacobian<Jac, false>(
              tape, tape.getZeroPosition(), tape.getPosition(),
              inputValues.data(), inputValues.size(),
              outputValues.data(), outputValues.size(),
              jac);
      }

      /**
       * @brief Evaluates the full Hessian of the recorded tape.
       *
       * The algorithm will select the best choice for the evaluation. Either a forward mode or reverse mode evaluation.
       * It will also use the vector mode if the underlying tape was configured with such a mode.
       *
       * @param[out] hes  The storage for the Hessian which is evaluated. Needs to have the correct size and should
       *                  be created with createHessian.
       * @param[out] jac  If also the Jacobian should be computed alongside the Hessian, a storage for the Jacobian can
       *                  be provided. Needs to have the correct size and should be created with createJacobian.
       *
       * @tparam Jac  Needs to implement the access operators (e.g. operator()) as in the Jacobian class.
       */
      template<typename Jac = DummyJacobian>
      void evalHessian(HessianType& hes, Jac& jac = StaticDummy<DummyJacobian>::dummy);

      /**
       * @brief Reevaluate the tape at a given position and compute the full Hessian at the new position.
       *
       * This method is a shortcut for calling evalPrimal and evalHessian
       *
       * @param[in]    x  The new seeding vector for the primal input variables. The sequence of variables is the same as
       *                  for the register input call. The vector should be created with createPrimalVectorInput.
       * @param[out] hes  The storage for the Hessian which is evaluated. Needs to have the correct size and should
       *                  be created with createHessian.
       * @param[out]   y  The result of the primal evaluation. The sequence of variables is the same as for the register
       *                  output call. If the pointer is a null pointer then the result is not stored. The vector should
       *                  be created with createPrimalVectorOutput.
       * @param[out] jac  If also the Jacobian should be computed alongside the Hessian, a storage for the Jacobian can
       *                  be provided. Needs to have the correct size and should be created with createJacobian.
       *
       * @tparam Jac  Needs to implement the access operators (e.g. operator()) as in the Jacobian class.
       */
      template<typename Jac = DummyJacobian>
      CODI_INLINE void evalHessianAt(Real const* x, HessianType& hes, Real* y = nullptr, Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        evalPrimal(x, y);

        cast().evalHessian(hes, jac);
      }

    protected:

      /**
       * Cast to the implementing class for the virtual template functions
       *
       * @return The instance of the implementing class.
       */
      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

      /**
       * @brief Create a gradient vector with the given size.
       *
       * @param[in] size  The size of the created vector.
       * @return The created vector.
       */
      GradientValue* createGradientVector(size_t size) {
        return new GradientValue[size];
      }

      /**
       * @brief Create a primal value vector with the given size.
       *
       * @param[in] size  The size of the created vector.
       * @return The created vector.
       */
      Real* createPrimalVector(size_t size) {
        return new Real[size];
      }

      /**
       * @brief Forward mode evaluation overwrite dirty adjoint vectors.
       */
      void changeStateToForwardEvaluation() {
        wasForwardEvaluated = true;

        // No cleanup to do
      }

      /**
       * @brief Default of this helper is the reverse state which leaves everything in a clean state, so we need to
       *        zero the adjoint vector.
       */
      void changeStateToReverseEvaluation() {
        if(wasForwardEvaluated) {
          // Forward evaluation leaves the adjoint vector dirty.

          tape.clearAdjoints();
        }

        wasForwardEvaluated = false;
      }
  };


  /**
   * @brief Empty implementation of TapeHelperBase.
   *
   * @tparam CoDiType  The CoDiPack type for which the helper works.
   */
  template<typename CoDiType>
  class TapeHelperNoImpl : public TapeHelperBase <CoDiType, TapeHelperNoImpl<CoDiType> > {
    public:
      typedef typename CoDiType::Real Real; /**< The floating point calculation type in the CoDiPack types. */

    private:
      using Impl = TapeHelperNoImpl<CoDiType>; /**< The definition this class */
    public:

      /**
       * * @brief No implementation yields compile time error.
       * @param[in]  x  Unused.
       * @param[out] y  Unused.
       */
      virtual void evalPrimal(Real const* x, Real* y = nullptr) = 0;

      /**
       * @brief No implementation yields compile time error.
       * @param[out] hes  Unused.
       * @param[out] jac  Unused.
       *
       * @tparam Jac Unused.
       */
      template<typename Jac = DummyJacobian>
      void evalHessian(typename TapeHelperBase<CoDiType, Impl>::HessianType& hes, Jac& jac = StaticDummy<DummyJacobian>::dummy);
  };

  /**
   * @brief Tape helper functionality for Jacobian tapes.
   *
   * Jacobian tapes do not support the change of the evaluation point of the tape, therefore an error is thrown if
   * these methods are called.
   *
   * @tparam CoDiType  The CoDiPack type for which the helper works.
   */
  template<typename CoDiType>
  class TapeHelperJacobi : public TapeHelperBase <CoDiType, TapeHelperJacobi<CoDiType> > {
    public:
      typedef typename CoDiType::Real Real; /**< The floating point calculation type in the CoDiPack types. */

    private:
      using Impl = TapeHelperJacobi<CoDiType>; /**< The definition this class */
    public:

      /**
       * @brief Throws an exception since primal evaluations are not support for Jacobian tapes.
       * @param[in]  x  Unused.
       * @param[out] y  Unused.
       */
      virtual void evalPrimal(Real const* x, Real* y = nullptr) {
        CODI_UNUSED(x);
        CODI_UNUSED(y);

        CODI_EXCEPTION(
          "No primal evaluation for Jacobian tapes. "
          "Please use codi::RealReversePrimal or codi::RealReversePrimalIndex types for this kind of functionality.");
      }

      /**
       * @brief Throws an exception since primal evaluations are not support for Jacobian tapes.
       * @param[out] hes  Unused.
       * @param[out] jac  Unused.
       *
       * @tparam Jac Unused.
       */
      template<typename Jac = DummyJacobian>
      void evalHessian(typename TapeHelperBase<CoDiType, Impl>::HessianType& hes, Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        CODI_UNUSED(hes);

        CODI_EXCEPTION(
          "No direct hessian evaluation for Jacobian tapes. "
          "Please use codi::RealReversePrimal or codi::RealReversePrimalIndex types for this kind of functionality "
          "or the EvaluationHelper class.");
      }
  };

  /**
   * @brief Tape helper functionality for primal values tapes.
   *
   * Primal values support the change of the evaluation point of the tape.
   *
   * @tparam CoDiType  The CoDiPack type for which the helper works.
   */
  template<typename CoDiType>
  class TapeHelperPrimal : public TapeHelperBase <CoDiType, TapeHelperPrimal<CoDiType> > {
    public:
      typedef typename CoDiType::Real Real; /**< The floating point calculation type in the CoDiPack types. */

    private:
      using Impl = TapeHelperJacobi<CoDiType>; /**< The definition this class */
    public:

      /**
       * \copydoc TapeHelperBase::evalPrimal
       */
      virtual void evalPrimal(Real const* x, Real* y = nullptr) {
        for(size_t j = 0; j < this->inputValues.size(); j += 1) {
          this->tape.primalValue(this->inputValues[j]) = x[j];
        }

        this->tape.evaluatePrimal();

        if(nullptr != y) {
          for(size_t i = 0; i < this->outputValues.size(); i += 1) {
            y[i] = this->tape.primalValue(this->outputValues[i]);
          }
        }
      }

      /**
       * \copydoc TapeHelperBase::evalHessian
       */
      template<typename Jac = DummyJacobian>
      void evalHessian(typename TapeHelperBase<CoDiType, Impl>::HessianType& hes, Jac& jac = StaticDummy<DummyJacobian>::dummy) {

        using Algo = Algorithms<CoDiType>;
        typename Algo::EvaluationType evalType = Algo::getEvaluationChoice(this->inputValues.size(), this->outputValues.size());

        if(Algo::EvaluationType::Forward == evalType) {
          this->changeStateToForwardEvaluation();
        } else if(Algo::EvaluationType::Reverse == evalType) {
          this->changeStateToReverseEvaluation();
        } else {
          CODI_EXCEPTION("Evaluation type not implemented.");
        }

        Algorithms<CoDiType>::computeHessianPrimalValueTape(
              this->tape, this->tape.getZeroPosition(), this->tape.getPosition(),
              this->inputValues.data(), this->inputValues.size(),
              this->outputValues.data(), this->outputValues.size(),
              hes, jac);
      }
  };

  /** \copydoc codi::TapeHelperBase */
  template<typename CoDiType, typename = void>
  class TapeHelper : public TapeHelperBase<CoDiType, TapeHelperNoImpl<CoDiType>> {};

  /** \copydoc codi::TapeHelperBase */
  template<typename CoDiType>
  class TapeHelper<CoDiType, enableIfJacobianTape<typename CoDiType::TapeType>> : public TapeHelperJacobi<CoDiType> {};

  /** \copydoc codi::TapeHelperBase */
  template<typename CoDiType>
  class TapeHelper<CoDiType, enableIfPrimalValueTape<typename CoDiType::TapeType>> : public TapeHelperPrimal<CoDiType> {};

}
