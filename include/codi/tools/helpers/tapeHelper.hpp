#pragma once

#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../traits/realTraits.hpp"
#include "../../traits/tapeTraits.hpp"
#include "../algorithms.hpp"
#include "../data/jacobian.hpp"
#include "../data/hessian.hpp"

#include <vector>

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Type, typename _Impl>
  struct TapeHelperBase {
    public:

      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
      using Impl = CODI_DECLARE_DEFAULT(_Impl, TapeHelperBase);

      using Real = typename Type::Real;
      using Identifier = typename Type::Identifier;
      using Gradient = typename Type::Gradient;

      using PassiveReal = typename RealTraits::PassiveReal<Real>;

      using JacobianType = Jacobian<PassiveReal>;
      using HessianType = Hessian<PassiveReal>;

    protected:

      using Tape = typename Type::Tape;

      Tape& tape;

      std::vector<Identifier> inputValues;
      std::vector<Identifier> outputValues;

      bool wasForwardEvaluated;

    public:

      TapeHelperBase() :
        tape(Type::getGlobalTape()),
        inputValues(),
        outputValues(),
        wasForwardEvaluated(false)
      {}

      virtual ~TapeHelperBase() {}

      Gradient* createGradientVectorInput() {
        return createGradientVector(getInputSize());
      }

      Gradient* createGradientVectorOutput() {
        return createGradientVector(getOutputSize());
      }

      JacobianType& createJacobian() {
        JacobianType* jacPointer = new JacobianType(getOutputSize(), getInputSize());

        return *jacPointer;
      }

      HessianType& createHessian() {
        HessianType* hesPointer = new HessianType(getOutputSize(), getInputSize());

        return *hesPointer;
      }

      Real* createPrimalVectorInput() {
        return createPrimalVector(getInputSize());
      }

      Real* createPrimalVectorOutput() {
        return createPrimalVector(getOutputSize());
      }

      void deleteGradientVector(Gradient* vec) {
        delete [] vec;

        vec = nullptr;
      }

      void deleteJacobian(JacobianType& jac) {
        JacobianType* jacPointer = &jac;

        delete jacPointer;
      }

      void deleteHessian(HessianType& hes) {
        HessianType* hesPointer = &hes;

        delete hesPointer;
      }

      void deletePrimalVector(Real* vec) {
        delete [] vec;

        vec = nullptr;
      }

      size_t getInputSize() {
        return inputValues.size();
      }

      size_t getOutputSize() {
        return outputValues.size();
      }

      void registerInput(Type& value) {
        tape.registerInput(value);
        inputValues.push_back(value.getIdentifier());
      }

      void registerOutput(Type& value) {
        tape.registerOutput(value);
        outputValues.push_back(value.getIdentifier());
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

      CODI_INLINE void evalForwardAt(Real const* x, Gradient const* x_d, Gradient* y_d, Real* y = nullptr) {
        evalPrimal(x, y);

        evalForward(x_d, y_d);
      }

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

      CODI_INLINE void evalReverseAt(Real const* x, Gradient const* y_b, Gradient* x_b, Real* y = nullptr) {
        evalPrimal(x, y);

        evalReverse(y_b, x_b);
      }

      CODI_INLINE void evalJacobian(JacobianType& jac) {
        JacobianConvertWrapper<JacobianType> wrapper(jac);

        evalJacobianGen(wrapper);
      }

      CODI_INLINE void evalJacobianAt(Real const* x, JacobianType& jac, Real* y = nullptr) {
        evalPrimal(x, y);

        evalJacobian(jac);
      }

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

        Algorithms<Type>::template computeJacobian<Jac, false>(
              tape, tape.getZeroPosition(), tape.getPosition(),
              inputValues.data(), inputValues.size(),
              outputValues.data(), outputValues.size(),
              jac);
      }

      template<typename Jac = DummyJacobian>
      void evalHessian(HessianType& hes, Jac& jac = StaticDummy<DummyJacobian>::dummy);

      template<typename Jac = DummyJacobian>
      CODI_INLINE void evalHessianAt(Real const* x, HessianType& hes, Real* y = nullptr, Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        evalPrimal(x, y);

        cast().evalHessian(hes, jac);
      }

    protected:

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

      Gradient* createGradientVector(size_t size) {
        return new Gradient[size];
      }

      Real* createPrimalVector(size_t size) {
        return new Real[size];
      }

      void changeStateToForwardEvaluation() {
        wasForwardEvaluated = true;

        // No cleanup to do
      }

      void changeStateToReverseEvaluation() {
        if (wasForwardEvaluated) {
          // Forward evaluation leaves the adjoint vector dirty.

          tape.clearAdjoints();
        }

        wasForwardEvaluated = false;
      }
  };


  template<typename _Type>
  struct TapeHelperNoImpl : public TapeHelperBase<_Type, TapeHelperNoImpl<_Type>> {
    public:

      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
      using Real = typename Type::Real;

    private:

      using Impl = TapeHelperNoImpl<Type>;

    public:

      virtual void evalPrimal(Real const* x, Real* y = nullptr) = 0;

      template<typename Jac = DummyJacobian>
      void evalHessian(typename TapeHelperBase<Type, Impl>::HessianType& hes, Jac& jac = StaticDummy<DummyJacobian>::dummy);
  };


  template<typename _Type>
  struct TapeHelperJacobi : public TapeHelperBase<_Type, TapeHelperJacobi<_Type>> {
    public:

      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
      using Real = typename Type::Real;

    private:

      using Impl = TapeHelperJacobi<Type>;

    public:

      virtual void evalPrimal(Real const* x, Real* y = nullptr) {
        CODI_UNUSED(x, y);

        CODI_EXCEPTION(
          "No primal evaluation for Jacobian tapes. "
          "Please use codi::RealReversePrimal or codi::RealReversePrimalIndex types for this kind of functionality.");
      }

      template<typename Jac = DummyJacobian>
      void evalHessian(typename TapeHelperBase<Type, Impl>::HessianType& hes, Jac& jac = StaticDummy<DummyJacobian>::dummy) {
        CODI_UNUSED(hes, jac);

        CODI_EXCEPTION(
          "No direct hessian evaluation for Jacobian tapes. "
          "Please use codi::RealReversePrimal or codi::RealReversePrimalIndex types for this kind of functionality "
          "or the EvaluationHelper class.");
      }
  };

  template<typename _Type>
  struct TapeHelperPrimal : public TapeHelperBase<_Type, TapeHelperPrimal<_Type>> {
    public:

      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
      using Real = typename Type::Real;

    private:

      using Impl = TapeHelperPrimal<Type>;

    public:

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

      template<typename Jac = DummyJacobian>
      void evalHessian(typename TapeHelperBase<Type, Impl>::HessianType& hes, Jac& jac = StaticDummy<DummyJacobian>::dummy) {

        using Algo = Algorithms<Type>;
        typename Algo::EvaluationType evalType = Algo::getEvaluationChoice(this->inputValues.size(), this->outputValues.size());

        if (Algo::EvaluationType::Forward == evalType) {
          this->changeStateToForwardEvaluation();
        } else if (Algo::EvaluationType::Reverse == evalType) {
          this->changeStateToReverseEvaluation();
        } else {
          CODI_EXCEPTION("Evaluation type not implemented.");
        }

        Algorithms<Type>::computeHessianPrimalValueTape(
              this->tape, this->tape.getZeroPosition(), this->tape.getPosition(),
              this->inputValues.data(), this->inputValues.size(),
              this->outputValues.data(), this->outputValues.size(),
              hes, jac);
      }
  };

  template<typename Type, typename = void>
  struct TapeHelper : public TapeHelperBase<Type, TapeHelperNoImpl<Type>> {};

  template<typename Type>
  struct TapeHelper<Type, TapeTraits::EnableIfJacobianTape<typename Type::Tape>> : public TapeHelperJacobi<Type> {};

  template<typename Type>
  struct TapeHelper<Type, TapeTraits::EnableIfPrimalValueTape<typename Type::Tape>> : public TapeHelperPrimal<Type> {};

}
