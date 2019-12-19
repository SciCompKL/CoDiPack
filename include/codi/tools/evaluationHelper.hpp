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

#include <array>
#include <vector>

#include "data/dummyValue.hpp"
#include "data/dummyVector.hpp"
#include "data/jacobian.hpp"
#include "data/hessian.hpp"
#include "data/vectorStorage.hpp"
#include "tapeHelper.hpp"
#include "../gradientTraits.hpp"
#include "../tapes/tapeTraits.hpp"
#include "../../codi.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Basic interface and data storage for all EvaluationHandle implementations.
   *
   * The class performs no resizing of the vectors. The initial size needs to be adequate for all calls to the
   * function object that the user will perform.
   *
   * @tparam             Func  The type for the function object which defines the evaluation logic.
   * @tparam         CoDiType  The CoDiPack type that is used for the derivative evaluation.
   * @tparam  InputVectorType  Vector used for the storage of input arguments.
   * @tparam OutputVectorType  Vector used for the storage of output arguments.
   */
  template <typename Func, typename CoDiType,
            template<typename> class InputVectorType,
            template<typename> class OutputVectorType>
  struct EvaluationHandleBase {

      using InputVector = InputVectorType<CoDiType>; /**< Input vector with CoDiPack type. */
      using OutputVector = OutputVectorType<CoDiType>;  /**< Input vector with CoDiPack type. */

    protected:
      size_t m; /**< Size of the output variables */
      size_t n; /**< Size of the input variables */

      Func& func; /**< The function object for the evaluations. */

      VectorStorage<InputVector> x;  /**< Storage for the input values as the CoDiPack types. */
      VectorStorage<OutputVector> y;  /**< Storage for the output values as the CoDiPack types. */

      DummyVector dummyVector;  /**< Dummy vector if no output is required. */
      DummyJacobian dummyJacobian;  /**< Dummy vector if no Jacobian is required. */

    public:

      /**
       * @brief Construct the base with the maximum amount of input and output values.
       *
       * @param[in] func  The function object for the evaluation.
       * @param[in]    m  The size of the output vector.
       * @param[in]    n  The size of the input vector.
       */
      EvaluationHandleBase(Func& func, size_t m, size_t n) : m(m), n(n), func(func), x(n), y(m),
          dummyVector(), dummyJacobian() {}

      /**
       * @brief Set the primal values from the user provided vector into the CoDiPack ones.
       *
       * @param[in] locX  The primal vector provided by the user.
       *
       * @tparam  VecX  The type of the user vector. It is accessed via the operator[].
       */
      template<typename VecX>
      void setAllPrimals(const VecX& locX);

      /**
       * @brief Store the primal values from the CoDiPack vector into the user vector.
       *
       * @param[out] locY  The primal vector provided by the user.
       *
       * @tparam  VecY  The type of the user vector. It is accessed via the operator[].
       */
      template<typename VecY>
      void getAllPrimals(VecY& locY);

      /**
       * @brief Perform a primal evaluation at the position provided in locX and store the result in locY.
       *
       * @param[in]  locX  The primal vector at which the evaluation is performed.
       * @param[out] locY  The primal vector for the storage of the outputs.
       *
       * @tparam  VecX  The type of the input user vector. It is accessed via the operator[].
       * @tparam  VecY  The type of the output user vector. It is accessed via the operator[].
       */
      template<typename VecX, typename VecY>
      void computePrimal(const VecX& locX, VecY& locY);

      /**
       * @brief Perform a Jacobian evaluation at the position provided in locX and store the result in jac and locY.
       *
       * VecY can be a dummy vector such that the values are not stored.
       *
       * @param[in]  locX  The primal vector at which the evaluation is performed.
       * @param[out]  jac  The Jacobian where the result is stored.
       * @param[out] locY  The primal vector for the storage of the outputs.
       *
       * @tparam  VecX  The type of the input user vector. It is accessed via the operator[].
       * @tparam   Jac  The type of the Jacobian.
       * @tparam  VecY  The type of the output user vector. It is accessed via the operator[].
       */
      template<typename VecX, typename Jac, typename VecY>
      void computeJacobian(const VecX& locX, Jac& jac, VecY& locY);


      /**
       * @brief Perform a Hessian evaluation at the position provided in locX and store the result in hes, jac and locY.
       *
       * VecY can be a dummy vector such that the values are not stored.
       * Jac can be a dummy Jacobian such that the values are not stored.
       *
       * @param[in]  locX  The primal vector at which the evaluation is performed.
       * @param[out]  hes  The Hessian where the result is stored.
       * @param[out]  jac  The Jacobian where the result is stored.
       * @param[out] locY  The primal vector for the storage of the outputs.
       *
       * @tparam  VecX  The type of the input user vector. It is accessed via the operator[].
       * @tparam   Hes  The type of the Hessian.
       * @tparam   Jac  The type of the Jacobian.
       * @tparam  VecY  The type of the output user vector. It is accessed via the operator[].
       */
      template<typename VecX, typename Hes, typename VecY, typename Jac>
      void computeHessian(const VecX& locX, Hes& hes, VecY& locY, Jac& jac);
    protected:

      /**
       * @brief Helper for the evaluation of the function object with the CoDiPack vectors. *
       */
      void eval() {
        func(x.vec, y.vec);
      }
  };

  /**
   * @brief Implementation for forward mode CoDiPack types of EvaluationHandleBase.
   *
   * \copydetails EvaluationHandleBase
   */
  template <typename Func, typename CoDiType,
            template<typename> class InputVectorType,
            template<typename> class OutputVectorType>
  struct ForwardHandle : public EvaluationHandleBase<Func, CoDiType, InputVectorType, OutputVectorType> {


      // Forward constructors of the base.
      using EvaluationHandleBase<Func, CoDiType, InputVectorType, OutputVectorType>::EvaluationHandleBase;

      /** \copydoc EvaluationHandleBase::setAllPrimals */
      template<typename VecX>
      void setAllPrimals(const VecX& locX) {
        codiAssert(locX.size() <= this->x.size());
        for(size_t j = 0; j < locX.size(); j += 1) {
          this->x[j] = locX[j];
        }
      }

      /** \copydoc EvaluationHandleBase::getAllPrimals */
      template<typename VecY>
      void getAllPrimals(VecY& locY) {
        codiAssert(locY.size() <= this->y.size());
        for(size_t i = 0; i < locY.size(); i += 1) {
          locY[i] = this->y[i].getValue();
        }
      }

      /** \copydoc EvaluationHandleBase::computePrimal */
      template<typename VecX, typename VecY>
      void computePrimal(const VecX& locX, VecY& locY) {
        setAllPrimals(locX);

        this->eval();

        getAllPrimals(locY);
      }

      /**
       * \copybrief EvaluationHandleBase::computeJacobian
       *
       * The vectorization is performed over the input vector. The function object is evaluated n/vecSize times.
       *
       * \copydetails EvaluationHandleBase::computeJacobian
       */
      template<typename VecX, typename Jac, typename VecY>
      void computeJacobian(const VecX& locX, Jac& jac, VecY& locY) {
        setAllPrimals(locX);

        // First order derivatives should always exist
        using GradientTraits1st = GradientValueTraits<typename CoDiType::GradientValue>;
        constexpr size_t VectorSizeFirstOrder = GradientTraits1st::getVectorSize();

        for(size_t j = 0; j < locX.size(); j+= VectorSizeFirstOrder) {
          for(size_t vecPos = 0; vecPos < VectorSizeFirstOrder && j + vecPos < locX.size(); vecPos += 1) {
            GradientTraits1st::at(this->x[j + vecPos].gradient(), vecPos) = 1.0;
          }

          this->eval();

          if(0 == j) {
            getAllPrimals(locY);
          }

          for(size_t i = 0; i < this->y.size(); i += 1) {
            for(size_t vecPos = 0; vecPos < VectorSizeFirstOrder && j + vecPos < locX.size(); vecPos += 1) {
              jac(i, j + vecPos) = GradientTraits1st::at(this->y[i].gradient(), vecPos);
            }
          }

          for(size_t vecPos = 0; vecPos < VectorSizeFirstOrder && j + vecPos < locX.size(); vecPos += 1) {
            GradientTraits1st::at(this->x[j + vecPos].gradient(), vecPos) = 0.0;
          }
        }
      }

      /**
       * \copybrief EvaluationHandleBase::computeHessian
       *
       * The vectorization is performed twice over the input vector. This evaluates the Hessian in a block wise fashion
       * for all output values.. The function object is evaluated n*n/(vecSize1 * vecSize2) times.
       *
       * \copydetails EvaluationHandleBase::computeHessian
       */
      template<typename VecX, typename Hes, typename VecY, typename Jac>
      void computeHessian(const VecX& locX, Hes& hes, VecY& locY, Jac& jac) {
        setAllPrimals(locX);

        // First order derivatives should always exist
        using GradientTraits1st = GradientValueTraits<typename CoDiType::GradientValue>;
        constexpr size_t VectorSizeFirstOrder = GradientTraits1st::getVectorSize();

        // Define these here since not all types have second order derivatives
        using GradientTraits2nd = GradientValueTraits<typename CoDiType::Real::GradientValue>;
        constexpr size_t VectorSizeSecondOrder = GradientTraits2nd::getVectorSize();


        for(size_t k = 0; k < locX.size(); k+= VectorSizeFirstOrder) {

          // Set derivatives from k to k + vecSize_k
          for(size_t vecPos = 0; vecPos < VectorSizeFirstOrder && k + vecPos < locX.size(); vecPos += 1) {
            GradientTraits1st::at(this->x[k + vecPos].gradient(), vecPos).value() = 1.0;
          }

          // The j = k init is no problem, it will evaluated slightly more elements around the diagonal
          for(size_t j = k; j < locX.size(); j+= VectorSizeSecondOrder) {

            // Set derivatives from j to j + vecSize_j
            for(size_t vecPos = 0; vecPos < VectorSizeSecondOrder && j + vecPos < locX.size(); vecPos += 1) {
              GradientTraits2nd::at(this->x[j + vecPos].value().gradient(), vecPos) = 1.0;
            }

            this->eval();

            if(0 == j && 0 == k) {
              getAllPrimals(locY);
            }

            // Extract all hessian values, this populates the hessian from (j,k) to (j + vecSize_j, k + vecSize_k).
            for(size_t i = 0; i < this->y.size(); i += 1) {
              for(size_t vecPos1st = 0; vecPos1st < VectorSizeFirstOrder && k + vecPos1st < locX.size(); vecPos1st += 1) {
                for(size_t vecPos2nd = 0; vecPos2nd < VectorSizeSecondOrder && j + vecPos2nd < locX.size(); vecPos2nd += 1) {
                  auto& firstGrad = GradientTraits1st::at(this->y[i].gradient(), vecPos1st);
                  auto& secondGrad = GradientTraits2nd::at(firstGrad.gradient(), vecPos2nd);

                  hes(i, j + vecPos2nd, k + vecPos1st) = secondGrad;
                  hes(i, k + vecPos1st, j + vecPos2nd) = secondGrad; // symmetry
                }
              }

              if(k == 0) {
                for(size_t vecPos = 0; vecPos < VectorSizeSecondOrder && j + vecPos < locX.size(); vecPos += 1) {
                  jac(i, j + vecPos) = GradientTraits2nd::at(this->y[i].value().gradient(), vecPos);
                }
              }
            }

            // Reset the derivative seeding
            for(size_t vecPos = 0; vecPos < VectorSizeSecondOrder && j + vecPos < locX.size(); vecPos += 1) {
              GradientTraits2nd::at(this->x[j + vecPos].value().gradient(), vecPos) = 0.0;
            }
          }

          // Reset the derivative seeding
          for(size_t vecPos = 0; vecPos < VectorSizeFirstOrder && k + vecPos < locX.size(); vecPos += 1) {
            GradientTraits1st::at(this->x[k + vecPos].gradient(), vecPos).value() = 0.0;
          }
        }
      }
  };

  /**
   * @brief Implementation for reverse mode CoDiPack types of EvaluationHandleBase.
   *
   * The basic evaluation for the reverse mode evaluation types.
   *
   * \copydetails EvaluationHandleBase
   */
  template <typename Func, typename CoDiType,
            template<typename> class InputVectorType,
            template<typename> class OutputVectorType>
  struct ReverseHandleBase : public EvaluationHandleBase<Func, CoDiType, InputVectorType, OutputVectorType> {

    protected:
      TapeHelper<CoDiType> th; /**< Tape helper handle for the algorithms */

    public:

      /** \copydoc EvaluationHandleBase::EvaluationHandleBase */
      ReverseHandleBase(Func& func, size_t m, size_t n) :
        EvaluationHandleBase<Func, CoDiType, InputVectorType, OutputVectorType>(func, m, n),
        th()
      {}

      /** \copydoc EvaluationHandleBase::setAllPrimals
       *
       * @param[in] reg  If the input variables should be registered.
       */
      template<typename VecX>
      void setAllPrimals(const VecX& locX, bool reg) {
        codiAssert(locX.size() <= this->x.size());
        for(size_t j = 0; j < locX.size(); j += 1) {
          this->x[j] = locX[j];

          if(reg) {
            th.registerInput(this->x[j]);
          }
        }
      }

      /**
       * \copydoc EvaluationHandleBase::getAllPrimals
       *
       * @param[in] reg  If the output variables should be registered.
       */
      template<typename VecY>
      void getAllPrimals(VecY& locY, bool reg) {
        codiAssert(locY.size() <= this->y.size());
        for(size_t i = 0; i < this->y.size(); i += 1) {
          if(reg) {
            th.registerOutput(this->y[i]);
          }

          locY[i] = this->y[i].getValue();
        }
      }

      /** \copydoc EvaluationHandleBase::computePrimal */
      template<typename VecX, typename VecY>
      void computePrimal(const VecX& locX, VecY& locY) {
        setAllPrimals(locX, false);

        this->eval();

        getAllPrimals(locY, false);
      }

      /**
       * \copybrief EvaluationHandleBase::computeJacobian
       *
       * The best mode is selected for the evaluation of the Jacobian. If \f$ n < m \f$ the the forward mode is used
       * for the Hessian evaluation and the function object is called \f$ n/vecSize \f$ times.
       * If \f$ m < n \f$ the reverse mode is used and the function object is called \f$ m/vecSize \f$ times.
       *
       * \copydetails EvaluationHandleBase::computeJacobian
       */
      template<typename VecX, typename Jac, typename VecY>
      void computeJacobian(const VecX& locX, Jac& jac, VecY& locY) {
        recordTape(locX, locY);

        th.evalJacobian(jac);
      }

      /**
       * \copybrief EvaluationHandleBase::computeHessian
       *
       * The best mode is selected for the evaluation of the Hessian. If \f$ n < m \f$ the the forward mode is used
       * for both gradient directions and the function object is called \f$ n*n/(vecSize1*vecSize2) \f$ times.
       * If \f$ m < n \f$ the reverse mode is used on the outer type and the forward mode on the inner type.
       * The function object is called \f$ (m/vecSize1) * (n/vecSize2) \f$ times.
       *
       * \copydetails EvaluationHandleBase::computeHessian
       */
      template<typename VecX, typename Hes, typename VecY, typename Jac>
      void computeHessian(const VecX& locX, Hes& hes, VecY& locY, Jac& jac);

    protected:

      /**
       * @brief Helper function for a tape recording.
       *
       * @param[in]  locX  The primal vector at which the evaluation is performed.
       * @param[out] locY  The primal vector for the storage of the outputs.
       *
       * @tparam  VecX  The type of the input user vector. It is accessed via the operator[].
       * @tparam  VecY  The type of the output user vector. It is accessed via the operator[].
       */
      template<typename VecX, typename VecY>
      void recordTape(const VecX& locX, VecY& locY) {
        th.startRecording();
        setAllPrimals(locX, true);

        this->eval();

        getAllPrimals(locY, true);
        th.stopRecording();
      }
  };

  /**
   * @brief Implementation for primal value reverse mode CoDiPack types of EvaluationHandleBase.
   *
   * This tape records the logic behind the function object once and then performs only primal, forward and reverse
   * tape evaluations.
   *
   * \copydetails EvaluationHandleBase
   */
  template <typename Func, typename CoDiType,
            template<typename> class InputVectorType,
            template<typename> class OutputVectorType>
  struct ReverseHandlePrimalValueTapes : public ReverseHandleBase<Func, CoDiType, InputVectorType, OutputVectorType> {

      // Use constructors of the base class.
      using ReverseHandleBase<Func, CoDiType, InputVectorType, OutputVectorType>::ReverseHandleBase;

      /**
       * \copybrief ReverseHandleBase::computeHessian
       *
       * For the primal value tape implementation, the tape is only recorded once and then evaluated multiple times.
       *
       * \copydetails ReverseHandleBase::computeHessian
       */
      template<typename VecX, typename Hes, typename VecY, typename Jac>
      void computeHessian(const VecX& locX, Hes& hes, VecY& locY, Jac& jac) {
        this->recordTape(locX, locY);

        this->th.evalHessian(hes, jac);
      }
  };

  /**
   * @brief Implementation for Jacobian reverse mode CoDiPack types of EvaluationHandleBase.
   *
   * This tape records the logic behind the function object new for every primal, forward and reverse
   * tape evaluations.
   *
   * \copydetails EvaluationHandleBase
   */
  template <typename Func, typename CoDiType,
            template<typename> class InputVectorType,
            template<typename> class OutputVectorType>
  struct ReverseHandleJacobiTapes : public ReverseHandleBase<Func, CoDiType, InputVectorType, OutputVectorType> {

      // Use constructors of the base class.
      using ReverseHandleBase<Func, CoDiType, InputVectorType, OutputVectorType>::ReverseHandleBase;

      /**
       * \copybrief ReverseHandleBase::computeHessian
       *
       * For the Jacobian tape implementation, a new tape is recorded for every evaluation.
       *
       * \copydetails ReverseHandleBase::computeHessian
       */
      template<typename VecX, typename Hes, typename VecY, typename Jac>
      void computeHessian(const VecX& locX, Hes& hes, VecY& locY, Jac& jac) {
        this->setAllPrimals(locX, false);

        Algorithms<CoDiType>::computeHessian(this->func, this->x.vec, this->y.vec, hes, jac);

        this->getAllPrimals(locY, false);
      }
  };

  /** \copydoc EvaluationHandleBase */
  template <typename Func, typename CoDiType,
            template<typename> class InputVectorType,
            template<typename> class OutputVectorType,
            typename = void>
  struct EvaluationHandle : public EvaluationHandleBase<Func, CoDiType, InputVectorType, OutputVectorType> {};

  /** \copydoc ForwardHandle */
  template <typename Func, typename CoDiType,
            template<typename> class InputVectorType,
            template<typename> class OutputVectorType>
  struct EvaluationHandle<Func, CoDiType, InputVectorType, OutputVectorType, enableIfForwardTape<typename CoDiType::TapeType>> :
      public ForwardHandle<Func, CoDiType, InputVectorType, OutputVectorType> {

      using ForwardHandle<Func, CoDiType, InputVectorType, OutputVectorType>::ForwardHandle;
  };

  /** \copydoc ReverseHandleJacobiTapes */
  template <typename Func, typename CoDiType,
            template<typename> class InputVectorType,
            template<typename> class OutputVectorType>
  struct EvaluationHandle<Func, CoDiType, InputVectorType, OutputVectorType, enableIfJacobianTape<typename CoDiType::TapeType>> :
      public ReverseHandleJacobiTapes<Func, CoDiType, InputVectorType, OutputVectorType> {

      using ReverseHandleJacobiTapes<Func, CoDiType, InputVectorType, OutputVectorType>::ReverseHandleJacobiTapes;
  };

  /** \copydoc ReverseHandlePrimalValueTapes */
  template <typename Func, typename CoDiType,
            template<typename> class InputVectorType,
            template<typename> class OutputVectorType>
  struct EvaluationHandle<Func, CoDiType, InputVectorType, OutputVectorType, enableIfPrimalValueTape<typename CoDiType::TapeType>> :
      public ReverseHandlePrimalValueTapes<Func, CoDiType, InputVectorType, OutputVectorType> {

      using ReverseHandlePrimalValueTapes<Func, CoDiType, InputVectorType, OutputVectorType>::ReverseHandlePrimalValueTapes;
  };


  /**
   * @brief Evaluate the primal, Jacobian and Hessian of function objects.
   *
   * This helper provides the means to easily evaluate derivatives for arbitrary function objects. These function objects
   * can be regular functions, lambda functions or structures were operator() is implemented.
   *
   * The nomenclature and mathematical definitions for the function, the Jacobian and the Hessian can be found in the
   * \ref FunctionDefinition "Algorithms" documentation. Function arguments in this class follow the same naming scheme
   * as in the referenced documentation.
   *
   * The algorithms will call the function objects with the vector of inputs and with the vector of outputs as arguments.
   * An example function definition is:
   * \code{.cpp}
   *  void func1(std::vector<ADType> const& x, std::vector<ADType>& y);
   *  // or
   *  void func1(std::array<ADType, n> const& x, std::array<ADType, m>& y);
   * \endcode
   *
   * x is the vector of input values and y is the vector of output values. ADType is the chosen CoDiPack type for the
   * function. For most users this definition will be enough. For more general examples please go to section
   * \ref AdvFuncObjDef.
   *
   * The CoDiPack type can be any ActiveReal type from CoDiPack, for example all types that are defined in codi.hpp. The
   * evaluation helper provides the default CoDiPack type definitions EvaluationHelper::JacobianComputationType and
   * EvaluationHelper::HessianComputationType. These two use the forward mode of algorithmic differentiation and are
   * more appropriate if \f$ m \f$ and \f$ n \f$ are small. They can also be used if \f$ n \f$ is smaller than
   * \f$ m \f$. For cases in in which the dimensions are larger and \f$ m \f$ is smaller than \f$ n \f$ the types
   * codi::JacobianComputationType and codi::HessianComputationType can be used. They use the reverse AD mode for the
   * computation.
   *
   * The most simple example of using the EvaluationHelper is:
   * \code{.cpp}
   *  using ADType = codi::EvaluationHelper::HessianComputationType;
   *
   *  void func(std::vector<ADType> const& x, std::vector<ADType>& y) {
   *    y[0] = x[0] + x[1];
   *    y[1] = x[0] - x[1];
   *    y[2] = x[0] * x[1];
   *    y[3] = x[0] / x[1];
   *  }
   *
   *  void main(int nargs, char** args) {
   *    std::vector<double> x = {3.0, 4.0};
   *    std::vector<double> y(4);
   *
   *    codi::EvaluationHelper eh;
   *    auto jac = eh.createJacobian(4,2);
   *    auto hes = eh.createHessian(4,2);
   *
   *    eh.evalPrimalAndJacobianAndHessian(func, x, y, jac, hes);
   *
   *    // output Jacobian and Hessian
   *    return 0;
   *  }
   * \endcode
   * Since we want to evaluate the Hessian, we use the Hessian type of the EvaluationHelper. The function is defined
   * with this type and the standard vector classes. In the main function we create the vector where we want to call the
   * function and then use the helper to create the storage functions for the Jacobian and Hessian. With
   * \f$ jac(j,i) \f$  the values can be accessed. For the Hessian the values can be accessed with \f$ hes(j,i,k) \f$
   * where \f$ j \f$ is the output dimension and \f$ i \f$ as well as \f$ k \f$ are the input dimensions.
   *
   * The evaluation helper class provides all combinations of evaluation choices that is: evalPrimal(), evalJacobian(),
   * evalHessian(), evalPrimalAndJacobian(), evalPrimalAndHessian(), evalJacobianAndHessian() and
   * evalPrimalAndJacobianAndHessian(). Each of these functions uses the default CoDiPack types in the evaluation helper.
   * In the cases where the primal is not stored, the user has to provide the number of output values manually.
   *
   * If the EvaluationHelper is used to evaluate the same function several times, a higher performance can be achived if
   * a handle for the evaluation is created up front and then used several times. The above example with the handle
   * creation would look like:
   * \code{.cpp}
   *  // func and ADType as before
   *
   *  void main(int nargs, char** args) {
   *    std::vector<double> x = {3.0, 4.0};
   *    std::vector<double> y(4);
   *
   *    codi::EvaluationHelper eh;
   *    auto jac = eh.createJacobian(4,2);
   *    auto hes = eh.createHessian(4,2);
   *
   *    auto handle = eh.createHandleDefault2nd(func, 4, 2);
   *
   *    eh.evalHandlePrimalAndJacobianAndHessian(handle, x, y, jac, hes);
   *
   *    // output Jacobian and Hessian
   *    return 0;
   *  }
   * \endcode
   * The logic nearly stayed the same. Instead of providing the function to the evaluation routine we create a handle
   * upfront and then use this handle in the evalHandle method. Each of the above mentioned eval routines has a
   * corresponding evalHandle method.
   *
   * Each of the create methods has similar create..Fixed method which uses the std::array type instead of the
   * std::vector type for the data management. An example with these methods would be:
   * \code{.cpp}
   *  using ADType = codi::EvaluationHelper::HessianComputationType;
   *
   *  void func(std::array<ADType> const& x, std::array<ADType>& y) {
   *    y[0] = x[0] + x[1];
   *    y[1] = x[0] - x[1];
   *    y[2] = x[0] * x[1];
   *    y[3] = x[0] / x[1];
   *  }
   *
   *  void main(int nargs, char** args) {
   *    std::array<double, 2> x = {3.0, 4.0};
   *    std::array<double, 4> y;
   *
   *    codi::EvaluationHelper eh;
   *    auto jac = eh.createJacobianFixed<4, 2>();
   *    auto hes = eh.createHessianFixed<4, 2>();
   *
   *    auto handle = eh.createHandleDefaultFixed2nd<4, 2>(func);
   *    eh.evalHandleAndJacobianAndHessian(handle, x, y, jac, hes);
   *
   *    // output Jacobian and Hessian
   *    return 0;
   *  }
   * \endcode
   *
   * Until now the default definition for the used CoDiPack types have been used. In order to use an arbitrary CoDiPack
   * type the createHandle(), createHandleFixed(), createHandleFull() methods can be used. The first one uses
   * std::vectors for the storage, the second method std::array and in the third the user can provide the storage class
   * as a template template parameter. The use case for the createHandle() method would look like:
   * \code{.cpp}
   *  auto handle = eh.createHandle<codi::RealReverse>(func, 4, 2);
   * \endcode
   *
   * \section AdvFuncObjDef Advanced function object definitions
   * The function object can also have a template argument for the evaluation type e.g.:
   * \code{.cpp}
   *  struct Func {
   *    template<typename T>
   *    void operator()(std::vector<T> const& x, std::vector<T>& y);
   *  };
   * \endcode
   *
   * There is also no need to specify std::vector as the array class e.g.:
   * \code{.cpp}
   *  struct Func {
   *    template<typename InVec, typename OutVec>
   *    void operator()(InVec const& x, OutVec& y);
   *  };
   * \endcode
   *
   */
  struct EvaluationHelper {

      /** @brief The default type used for first order gradient computation. It is defined as forward vector AD mode of
       * size 4. */
      using JacobianComputationType = RealForwardVec<4>;
      /** @brief The default type used for second order gradient computation. It is defined as forward vector over
       * forward vector AD mode of size 4 and 4. */
      using HessianComputationType = RealForwardGen<RealForwardVec<4>, Direction<RealForwardVec<4>, 4>>;

      /**
       * @brief Type for the default handle for first order gradient computations with a variable vector size.
       * @tparam Func  The function object which describes the evaluation function.
       */
      template <typename Func>
      using DefaultHandle = ForwardHandle<Func, JacobianComputationType, adapters::StdVector, adapters::StdVector>;

      /**
       * @brief Type for the default handle for second order gradient computations with a variable vector size.
       * @tparam Func  The function object which describes the evaluation function.
       */
      template <typename Func>
      using DefaultHandle2nd = ForwardHandle<Func, HessianComputationType, adapters::StdVector, adapters::StdVector>;

      /**
       * @brief Type for the default handle for first order gradient computations with a fixed vector size.
       * @tparam Func  The function object which describes the evaluation function.
       * @tparam    m  The size of the output vector.
       * @tparam    n  The size of the input vector.
       */
      template <typename Func, size_t m, size_t n>
      using DefaultHandleFixed = ForwardHandle<Func, JacobianComputationType, adapters::StdArray<n>::template Type, adapters::StdArray<m>::template Type>;

      /**
       * @brief Type for the default handle for second order gradient computations with a fixed vector size.
       * @tparam Func  The function object which describes the evaluation function.
       * @tparam    m  The size of the output vector.
       * @tparam    n  The size of the input vector.
       */
      template <typename Func, size_t m, size_t n>
      using DefaultHandleFixed2nd = ForwardHandle<Func, HessianComputationType , adapters::StdArray<n>::template Type, adapters::StdArray<m>::template Type>;

      /**
       * @brief Helper function for the creation of a default first order evaluation handle with a variable vector size.
       *
       * @param[in] func  The function object for the evaluation.
       * @param[in]    m  The size of the output vector.
       * @param[in]    n  The size of the input vector.
       *
       * @tparam Func  The function object which describes the evaluation function.
       */
      template<typename Func>
      static CODI_INLINE DefaultHandle<Func> createHandleDefault(Func& func, size_t m, size_t n) {
        return DefaultHandle<Func>(func, m, n);
      }

      /**
       * @brief Helper function for the creation of a default first order evaluation handle with a fixed vector size.
       *
       * @param[in] func  The function object for the evaluation.
       *
       * @tparam    m  The size of the output vector.
       * @tparam    n  The size of the input vector.
       * @tparam Func  The function object which describes the evaluation function.
       */
      template<size_t m, size_t n, typename Func>
      static CODI_INLINE DefaultHandleFixed<Func, m, n> createHandleDefaultFixed(Func& func) {
        return DefaultHandleFixed<Func, m, n>(func, m, n);
      }

      /**
       * @brief Helper function for the creation of a default second order evaluation handle with a variable vector size.
       *
       * @param[in] func  The function object for the evaluation.
       * @param[in]    m  The size of the output vector.
       * @param[in]    n  The size of the input vector.
       *
       * @tparam Func  The function object which describes the evaluation function.
       */
      template<typename Func>
      static CODI_INLINE DefaultHandle2nd<Func> createHandleDefault2nd(Func& func, size_t m, size_t n) {
        return DefaultHandle2nd<Func>(func, m, n);
      }

      /**
       * @brief Helper function for the creation of a default second order evaluation handle with a fixed vector size.
       *
       * @param[in] func  The function object for the evaluation.
       *
       * @tparam    m  The size of the output vector.
       * @tparam    n  The size of the input vector.
       * @tparam Func  The function object which describes the evaluation function.
       */
      template<size_t m, size_t n, typename Func>
      static CODI_INLINE DefaultHandleFixed2nd<Func, m, n> createHandleDefaultFixed2nd(Func& func) {
        return DefaultHandleFixed2nd<Func, m, n>(func, m, n);
      }

      /**
       * @brief Helper function for the creation of an evaluation handle with the specified CoDiPack type and a variable
       * vector size.
       *
       * The CoDiPack type can be an arbitrary one:
       * \code{.cpp}
       *  auto handle = codi::EvaluationHelper::createHandle<codi::RealReverse>(func, m, n);
       * \endcode
       *
       * @param[in] func  The function object for the evaluation.
       * @param[in]    m  The size of the output vector.
       * @param[in]    n  The size of the input vector.
       *
       * @tparam CoDiType  An arbitrary CoDiPack type based on ActiveReal. All definitions in codi.hpp are supported.
       *                   For user developed tapes some sub classes need to be specialized.
       * @tparam     Func  The function object which describes the evaluation function.
       */
      template<typename CoDiType, typename Func>
      static CODI_INLINE EvaluationHandle<Func, CoDiType, adapters::StdVector, adapters::StdVector> createHandle(Func& func, size_t m, size_t n) {
        return EvaluationHandle<Func, CoDiType, adapters::StdVector, adapters::StdVector>(func, m, n);
      }

      /**
       * @brief Helper function for the creation of an evaluation handle with the specified CoDiPack type and a fixed
       * vector size.
       *
       * The CoDiPack type can be an arbitrary one:
       * \code{.cpp}
       *  auto handle = codi::EvaluationHelper::createHandleFixed<codi::RealReverse, m, n>(func);
       * \endcode
       *
       * @param[in] func  The function object for the evaluation.
       *
       * @tparam CoDiType  An arbitrary CoDiPack type based on ActiveReal. All definitions in codi.hpp are supported.
       *                   For user developed tapes some sub classes need to be specialized.
       * @tparam     Func  The function object which describes the evaluation function.
       * @tparam        m  The size of the output vector.
       * @tparam        n  The size of the input vector.
       */
      template<typename CoDiType, size_t m, size_t n, typename Func>
      static CODI_INLINE EvaluationHandle<Func, CoDiType, adapters::StdArray<n>::template Type, adapters::StdArray<m>::template Type> createHandleFixed(Func& func) {
        return EvaluationHandle<Func, CoDiType, adapters::StdArray<n>::template Type, adapters::StdArray<m>::template Type>(func, m, n);
      }

      /**
       * @brief Helper function for the creation of an evaluation handle with the specified CoDiPack type and a vector
       * type specified by the user.
       *
       * The CoDiPack type can be an arbitrary one, the vector types need to be template template parameters:
       * \code{.cpp}
       *  auto handle = codi::EvaluationHelper::createHandleFull<codi::RealReverse, std::adapters::StdVector, codi::adapters::StdArray<m>::template Type>(func, m, n);
       * \endcode
       *
       * @param[in] func  The function object for the evaluation.
       * @param[in]    m  The size of the output vector.
       * @param[in]    n  The size of the input vector.
       *
       * @tparam         CoDiType  An arbitrary CoDiPack type based on ActiveReal. All definitions in codi.hpp are
       *                           supported. For user developed tapes some sub classes need to be specialized.
       * @tparam             Func  The function object which describes the evaluation function.
       * @tparam  InputVectorType  The vector type for vectors of input variables.
       * @tparam OutputVectorType  The vector type for vectors of output variables.
       */
      template<typename CoDiType,
               template<typename> class InputVectorType,
               template<typename> class OutputVectorType,
               typename Func>
      static CODI_INLINE EvaluationHandle<Func, CoDiType, InputVectorType, OutputVectorType> createHandleFull(Func& func, size_t m, size_t n) {
        return EvaluationHandle<Func, CoDiType, InputVectorType, OutputVectorType>(func, m, n);
      }

      /**
       * @brief Create a Jacobian with the given size.
       *
       * @param[in]    m  The size of the output vector.
       * @param[in]    n  The size of the input vector.
       *
       * @tparam T  The storage type of the Jacobian.
       */
      template<typename T = double>
      static CODI_INLINE Jacobian<adapters::StdVector<T>> createJacobian(size_t m, size_t n) {
        return Jacobian<adapters::StdVector<T>>(m, n);
      }

      /**
       * @brief Create a Jacobian with a fixed size.
       *
       * @tparam T  The storage type of the Jacobian.
       * @tparam m  The size of the output vector.
       * @tparam n  The size of the input vector.
       */
      template<size_t m, size_t n, typename T = double>
      static CODI_INLINE Jacobian<typename adapters::StdArray<m*n>::template Type<T>> createJacobianFixed() {
        return Jacobian<typename adapters::StdArray<m*n>::template Type<T>>(m, n);
      }

      /**
       * @brief Create a Hessian with the given size.
       *
       * @param[in]    m  The size of the output vector.
       * @param[in]    n  The size of the input vector.
       *
       * @tparam T  The storage type of the Hessian.
       */
      template<typename T = double>
      static CODI_INLINE Hessian<adapters::StdVector<T>> createHessian(size_t m, size_t n) {
        return Hessian<adapters::StdVector<T>>(m, n);
      }


      /**
       * @brief Create a Hessian with a fixed size.
       *
       * @tparam T  The storage type of the Hessian.
       * @tparam        m  The size of the output vector.
       * @tparam        n  The size of the input vector.
       */
      template<size_t m, size_t n, typename T = double>
      static CODI_INLINE Hessian<typename adapters::StdArray<m*n*n>::template Type<T>> createHessianFixed() {
        return Hessian<typename adapters::StdArray<m*n*n>::template Type<T>>(m, n);
      }

      /**
       * @brief Perform a primal evaluation of the function object with the default first order type.
       *
       * @param[in] func  The function object for the evaluation.
       * @param[in]    x  The vector with the primal values where the function object is evaluated.
       * @param[out]   y  The vector for the result of the primal function evaluation. The vector needs to have the
       *                  correct size allocated.
       *
       * @tparam  Func  The function object for the evaluation.
       * @tparam  VecX  The vector type for the input values.
       * @tparam  VecY  The vector type for the output values.
       */
      template<typename Func, typename VecX, typename VecY>
      static CODI_INLINE void evalPrimal(
          Func& func,
          const VecX& x,
          VecY& y) {
        auto h = createHandleDefault(func, y.size(), x.size());
        evalHandlePrimal(h, x, y);
      }

      /**
       * @brief Compute the Jacobian of the evaluation procedure in the function object.
       *
       * @param[in]  func  The function object for the evaluation.
       * @param[in]     x  The vector with the primal values where the function object is evaluated.
       * @param[in] ySize  The size of the output vector.
       * @param[out]  jac  The Jacobian in which the values are stored.
       *
       * @tparam  Func  The function object for the evaluation.
       * @tparam  VecX  The vector type for the input values.
       * @tparam   Jac  The type of the Jacobian where the results are stored.
       */
      template<typename Func, typename VecX, typename Jac>
      static CODI_INLINE void evalJacobian(
          Func& func,
          const VecX& x,
          const size_t ySize,
          Jac& jac) {
        auto h = createHandleDefault(func, ySize, x.size());
        evalHandleJacobian(h, x, jac);
      }

      /**
       * @brief Compute the Hessian of the evaluation procedure in the function object.
       *
       * @param[in]  func  The function object for the evaluation.
       * @param[in]     x  The vector with the primal values where the function object is evaluated.
       * @param[in] ySize  The size of the output vector.
       * @param[out]  hes  The Hessian in which the values are stored.
       *
       * @tparam  Func  The function object for the evaluation.
       * @tparam  VecX  The vector type for the input values.
       * @tparam   Hes  The type of the Jacobian where the results are stored.
       */
      template<typename Func, typename VecX, typename Hes>
      static CODI_INLINE void evalHessian(
          Func& func,
          const VecX& x,
          const size_t ySize,
          Hes& hes) {
        auto h = createHandleDefault2nd(func, ySize, x.size());
        evalHandleHessian(h, x, hes);
      }

      /**
       * @brief Compute the Jacobian of the evaluation procedure in the function object.
       *
       * In this method the primal values are also stored.
       *
       * @param[in]  func  The function object for the evaluation.
       * @param[in]     x  The vector with the primal values where the function object is evaluated.
       * @param[out]    y  The vector for the result of the primal function evaluation. The vector needs to have the
       *                   correct size allocated.
       * @param[out]  jac  The Jacobian in which the values are stored.
       *
       * @tparam  Func  The function object for the evaluation.
       * @tparam  VecX  The vector type for the input values.
       * @tparam  VecY  The vector type for the output values.
       * @tparam   Jac  The type of the Jacobian where the results are stored.
       */
      template<typename Func, typename VecX, typename VecY, typename Jac>
      static CODI_INLINE void evalPrimalAndJacobian(
          Func& func,
          const VecX& x,
          VecY& y,
          Jac& jac) {
        auto h = createHandleDefault(func, y.size(), x.size());
        evalHandlePrimalAndJacobian(h, x, y, jac);
      }

      /**
       * @brief Compute the Hessian of the evaluation procedure in the function object.
       *
       * In this method the primal values are also stored.
       *
       * @param[in]  func  The function object for the evaluation.
       * @param[in]     x  The vector with the primal values where the function object is evaluated.
       * @param[out]    y  The vector for the result of the primal function evaluation. The vector needs to have the
       *                   correct size allocated.
       * @param[out]  hes  The Hessian in which the values are stored.
       *
       * @tparam  Func  The function object for the evaluation.
       * @tparam  VecX  The vector type for the input values.
       * @tparam  VecY  The vector type for the output values.
       * @tparam   Hes  The type of the Jacobian where the results are stored.
       */
      template<typename Func, typename VecX, typename VecY, typename Hes>
      static CODI_INLINE void evalPrimalAndHessian(
          Func& func,
          const VecX& x,
          VecY& y,
          Hes& hes) {
        auto h = createHandleDefault2nd(func, y.size(), x.size());
        evalHandlePrimalAndHessian(h, x, y, hes);
      }

      /**
       * @brief Compute the Hessian of the evaluation procedure in the function object.
       *
       * In this method the primal values and the Jacobian are also stored.
       *
       * @param[in]  func  The function object for the evaluation.
       * @param[in]     x  The vector with the primal values where the function object is evaluated.
       * @param[out]    y  The vector for the result of the primal function evaluation. The vector needs to have the
       *                   correct size allocated.
       * @param[out]  jac  The Jacobian in which the values are stored.
       * @param[out]  hes  The Hessian in which the values are stored.
       *
       * @tparam  Func  The function object for the evaluation.
       * @tparam  VecX  The vector type for the input values.
       * @tparam  VecY  The vector type for the output values.
       * @tparam   Jac  The type of the Jacobian where the results are stored.
       * @tparam   Hes  The type of the Jacobian where the results are stored.
       */
      template<typename Func, typename VecX, typename VecY, typename Jac, typename Hes>
      static CODI_INLINE void evalPrimalAndJacobianAndHessian(
          Func& func,
          const VecX& x,
          VecY& y,
          Jac& jac,
          Hes& hes) {
        auto h = createHandleDefault2nd(func, y.size(), x.size());
        evalHandlePrimalAndJacobianAndHessian(h, x, y, jac, hes);
      }


      /**
       * @brief Compute the Hessian of the evaluation procedure in the function object.
       *
       * In this method the primal values and the Jacobian are also stored.
       *
       * @param[in]  func  The function object for the evaluation.
       * @param[in]     x  The vector with the primal values where the function object is evaluated.
       * @param[in] ySize  The size of the output vector.
       * @param[out]  jac  The Jacobian in which the values are stored.
       * @param[out]  hes  The Hessian in which the values are stored.
       *
       * @tparam  Func  The function object for the evaluation.
       * @tparam  VecX  The vector type for the input values.
       * @tparam   Jac  The type of the Jacobian where the results are stored.
       * @tparam   Hes  The type of the Jacobian where the results are stored.
       */
      template<typename Func, typename VecX, typename Jac, typename Hes>
      static CODI_INLINE void evalJacobianAndHessian(
          Func& func,
          const VecX& x,
          size_t ySize,
          Jac& jac,
          Hes& hes) {
        auto h = createHandleDefault2nd(func, ySize, x.size());
        evalHandleJacobianAndHessian(h, x, jac, hes);
      }

      /**
       * @brief Perform a primal evaluation of the function object stored in the handle.
       *
       * @param[in] handle  The handle with all data for the evaluation.
       * @param[in]      x  The vector with the primal values where the function object is evaluated.
       * @param[out]     y  The vector for the result of the primal function evaluation. The vector needs to have the
       *                    correct size allocated.
       *
       * @tparam  Func  The function object for the evaluation.
       * @tparam  VecX  The vector type for the input values.
       * @tparam  VecY  The vector type for the output values.
       */
      template<typename Handle, typename VecX, typename VecY>
      static CODI_INLINE void evalHandlePrimal(
          Handle& handle,
          const VecX& x,
          VecY& y) {
        handle.computePrimal(x, y);
      }

      /**
       * @brief Compute the Jacobian of the evaluation procedure of the function object stored in the handle.
       *
       * @param[in] handle  The handle with all data for the evaluation.
       * @param[in]      x  The vector with the primal values where the function object is evaluated.
       * @param[out]   jac  The Jacobian in which the values are stored.
       *
       * @tparam  Handle  The handle type for the data storage and the evaluation.
       * @tparam    VecX  The vector type for the input values.
       * @tparam     Jac  The type of the Jacobian where the results are stored.
       */
      template<typename Handle, typename VecX, typename Jac>
      static CODI_INLINE void evalHandleJacobian(
          Handle& handle,
          const VecX& x,
          Jac& jac) {
        DummyVector dv;
        handle.computeJacobian(x, jac, dv);
      }

      /**
       * @brief Compute the Hessian of the evaluation procedure of the function object stored in the handle.
       *
       * @param[in] handle  The handle with all data for the evaluation.
       * @param[in]      x  The vector with the primal values where the function object is evaluated.
       * @param[out]   hes  The Hessian in which the values are stored.
       *
       * @tparam  Handle  The handle type for the data storage and the evaluation.
       * @tparam    VecX  The vector type for the input values.
       * @tparam     Hes  The type of the Jacobian where the results are stored.
       */
      template<typename Handle, typename VecX, typename Hes>
      static CODI_INLINE void evalHandleHessian(
          Handle& handle,
          const VecX& x,
          Hes& hes) {
        DummyVector dv;
        DummyJacobian dj;
        handle.computeHessian(x, hes, dv, dj);
      }

      /**
       * @brief Compute the Jacobian of the evaluation procedure of the function object stored in the handle.
       *
       * In this method the primal values are also stored.
       *
       * @param[in] handle  The handle with all data for the evaluation.
       * @param[in]      x  The vector with the primal values where the function object is evaluated.
       * @param[out]     y  The vector for the result of the primal function evaluation. The vector needs to have the
       *                    correct size allocated.
       * @param[out]   jac  The Jacobian in which the values are stored.
       *
       * @tparam  Handle  The handle type for the data storage and the evaluation.
       * @tparam    VecX  The vector type for the input values.
       * @tparam    VecY  The vector type for the output values.
       * @tparam     Jac  The type of the Jacobian where the results are stored.
       */
      template<typename Handle, typename VecX, typename VecY, typename Jac>
      static CODI_INLINE void evalHandlePrimalAndJacobian(
          Handle& handle,
          const VecX& x,
          VecY& y,
          Jac& jac) {
        handle.computeJacobian(x, jac, y);
      }

      /**
       * @brief Compute the Hessian of the evaluation procedure of the function object stored in the handle.
       *
       * In this method the primal values are also stored.
       *
       * @param[in] handle  The handle with all data for the evaluation.
       * @param[in]      x  The vector with the primal values where the function object is evaluated.
       * @param[out]     y  The vector for the result of the primal function evaluation. The vector needs to have the
       *                    correct size allocated.
       * @param[out]   hes  The Hessian in which the values are stored.
       *
       * @tparam  Handle  The handle type for the data storage and the evaluation.
       * @tparam    VecX  The vector type for the input values.
       * @tparam    VecY  The vector type for the output values.
       * @tparam     Hes  The type of the Jacobian where the results are stored.
       */
      template<typename Handle, typename VecX, typename VecY, typename Hes>
      static CODI_INLINE void evalHandlePrimalAndHessian(
          Handle& handle,
          const VecX& x,
          VecY& y,
          Hes& hes) {
        DummyJacobian dj;
        handle.computeHessian(x, hes, y, dj);
      }

      /**
       * @brief Compute the Hessian of the evaluation procedure of the function object stored in the handle.
       *
       * In this method the primal values and the Jacobian are also stored.
       *
       * @param[in] handle  The handle with all data for the evaluation.
       * @param[in]      x  The vector with the primal values where the function object is evaluated.
       * @param[out]     y  The vector for the result of the primal function evaluation. The vector needs to have the
       *                    correct size allocated.
       * @param[out]   jac  The Jacobian in which the values are stored.
       * @param[out]   hes  The Hessian in which the values are stored.
       *
       * @tparam  Handle  The handle type for the data storage and the evaluation.
       * @tparam    VecX  The vector type for the input values.
       * @tparam    VecY  The vector type for the output values.
       * @tparam     Jac  The type of the Jacobian where the results are stored.
       * @tparam     Hes  The type of the Jacobian where the results are stored.
       */
      template<typename Handle, typename VecX, typename VecY, typename Jac, typename Hes>
      static CODI_INLINE void evalHandlePrimalAndJacobianAndHessian(
          Handle& handle,
          const VecX& x,
          VecY& y,
          Jac& jac,
          Hes& hes) {
        handle.computeHessian(x, hes, y, jac);
      }

      /**
       * @brief Compute the Hessian of the evaluation procedure of the function object stored in the handle.
       *
       * In this method the primal values and the Jacobian are also stored.
       *
       * @param[in] handle  The handle with all data for the evaluation.
       * @param[in]      x  The vector with the primal values where the function object is evaluated.
       * @param[out]   jac  The Jacobian in which the values are stored.
       * @param[out]   hes  The Hessian in which the values are stored.
       *
       * @tparam  Handle  The handle type for the data storage and the evaluation.
       * @tparam    VecX  The vector type for the input values.
       * @tparam     Jac  The type of the Jacobian where the results are stored.
       * @tparam     Hes  The type of the Jacobian where the results are stored.
       */
      template<typename Handle, typename VecX, typename Jac, typename Hes>
      static CODI_INLINE void evalHandleJacobianAndHessian(
          Handle& handle,
          const VecX& x,
          Jac& jac,
          Hes& hes) {
        DummyVector dv;
        handle.computeHessian(x, hes, dv, jac);
      }
  };

}
