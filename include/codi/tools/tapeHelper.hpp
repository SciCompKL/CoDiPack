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

  template<typename CoDiType>
  class TapeHelperBase {
    protected:
      typedef typename CoDiType::Real Real; /**< The floating point calculation type in the CoDiPack types. */
      typedef typename CoDiType::GradientData GradientData; /**< The type for the identification of gradients. */
      typedef typename CoDiType::GradientValue GradientValue; /**< The type for the gradient computation */

      typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

      using JacobianType = Jacobian<std::vector<PassiveReal>>;

      /** The type of the tape implementation. */
      //typedef ReverseTapeInterface<Real, GradientData, GradientValue, typename CoDiType::TapeType, typename CoDiType::TapeType::Position> Tape;
      typedef typename CoDiType::TapeType Tape;

      Tape& tape;

      std::vector<GradientData> inputValues;
      std::vector<GradientData> outputValues;

      bool wasForwardEvaluated;

    public:

      TapeHelperBase() :
        tape(CoDiType::getGlobalTape()),
        inputValues(),
        outputValues(),
        wasForwardEvaluated(false)
      {}

      virtual ~TapeHelperBase() {}

      GradientValue* createGradientVectorInput() {
        return createGradientVector(getInputSize());
      }

      GradientValue* createGradientVectorOutput() {
        return createGradientVector(getOutputSize());
      }

      void deleteGradientVector(GradientValue* vec) {
        delete [] vec;
      }

      JacobianType& createJacobian() {
        JacobianType* jacPointer = new JacobianType(getOutputSize(), getInputSize());

        return *jacPointer;
      }

      void deleteJacobian(JacobianType& jac) {
        JacobianType* jacPointer = *jac;

        delete jacPointer;
      }

      size_t getInputSize() {
        return inputValues.size();
      }

      size_t getOutputSize() {
        return outputValues.size();
      }

      void registerInput(CoDiType& value) {
        tape.registerInput(value);
        inputValues.push_back(value.getGradientData());
      }

      void registerOutput(CoDiType& value) {
        tape.registerOutput(value);
        outputValues.push_back(value.getGradientData());
      }

      void startRecording() {
        tape.reset();
        inputValues.clear();
        outputValues.clear();

        tape.setActive();
      }

      void stopRecording() {
        tape.setPassive();
      }

      virtual void evalPrimal(Real const* x, Real* y = nullptr) = 0;

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

      CODI_INLINE void evalReverse(GradientValue const* y_d, GradientValue* x_d) {
        changeStateToReverseEvaluation();

        for(size_t i = 0; i < outputValues.size(); i += 1) {
          tape.setGradient(outputValues[i], y_d[i]);
        }

        tape.evaluate();


        for(size_t j = 0; j < inputValues.size(); j += 1) {
          x_d[j] = tape.getGradient(inputValues[j]);
          tape.setGradient(inputValues[j], GradientValue());
        }

        if(!ZeroAdjointReverse) {
          tape.clearAdjoints();
        }
      }

      CODI_INLINE void evalJacobian(JacobianType& jac) {
        JacobianConvertWrapper<JacobianType> wrapper(jac);

        evalJacobianGen(wrapper);
      }

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

      CODI_INLINE void evalForwardAt(Real const* x, GradientValue const* x_d, GradientValue* y_d, Real* y = nullptr) {
        evalPrimal(x, y);

        evalForward(x_d, y_d);
      }

      CODI_INLINE void evalReverseAt(Real const* x, GradientValue const* y_d, GradientValue* x_d, Real* y = nullptr) {
        evalPrimal(x, y);

        evalReverse(y_d, x_d);
      }

      CODI_INLINE void evalJacobianAt(Real const* x, JacobianType& jac, Real* y = nullptr) {
        evalPrimal(x, y);

        evalJacobian(jac);
      }

    private:
      GradientValue* createGradientVector(size_t size) {
        return new GradientValue[size];
      }

      void changeStateToForwardEvaluation() {
        wasForwardEvaluated = true;

        // No cleanup to do
      }

      void changeStateToReverseEvaluation() {
        if(wasForwardEvaluated) {
          // Forward evaluation leaves the adjoint vector dirty.

          tape.clearAdjoints();
        }

        wasForwardEvaluated = false;
      }
  };

  template<typename CoDiType>
  class TapeHelperJacobi : public TapeHelperBase <CoDiType> {
    private:
      typedef typename CoDiType::Real Real; /**< The floating point calculation type in the CoDiPack types. */


      virtual void evalPrimal(Real const* x, Real* y = nullptr) {
        CODI_UNUSED(x);
        CODI_UNUSED(y);

        CODI_EXCEPTION(
          "No primal evaluation for Jacbian tapes. "
          "Please use codi::RealReversePrimal or codi::RealReversePrimalIndex types for this kind of functionality.");
      }
  };

  template<typename CoDiType>
  class TapeHelperPrimal : public TapeHelperBase <CoDiType> {
    private:
      typedef typename CoDiType::Real Real; /**< The floating point calculation type in the CoDiPack types. */

    public:

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
  };

  template<typename CoDiType, typename = void>
  class TapeHelper;

  template<typename CoDiType>
  class TapeHelper<CoDiType, enableIfJacobianTape<typename CoDiType::TapeType>> : public TapeHelperJacobi<CoDiType> {};

  template<typename CoDiType>
  class TapeHelper<CoDiType, enableIfPrimalValueTape<typename CoDiType::TapeType>> : public TapeHelperJacobi<CoDiType> {};

}
