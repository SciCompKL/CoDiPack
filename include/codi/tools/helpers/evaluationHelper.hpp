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

#include <array>
#include <vector>

#include "../../../codi.hpp"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../misc/constructVector.hpp"
#include "../../traits/gradientTraits.hpp"
#include "../../traits/tapeTraits.hpp"
#include "../data/dummy.hpp"
#include "../data/hessian.hpp"
#include "../data/jacobian.hpp"
#include "tapeHelper.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Basic interface and data storage for all EvaluationHandle implementations.
   *
   * The class performs no resizing of the vectors. The initial sizes need to be adequate for all calls to the
   * function object that the user will perform.
   *
   * @tparam T_Func  The type of the function object which defines the evaluation logic.
   * @tparam T_Type  The CoDiPack type that is used for the derivative evaluation.
   * @tparam T_InputStore  Vector used for the storage of input arguments.
   * @tparam T_OutputStore  Vector used for the storage of output arguments.
   */
  template<typename T_Func, typename T_Type, typename T_InputStore = std::vector<T_Type>,
           typename T_OutputStore = std::vector<T_Type>>
  struct EvaluationHandleBase {
    public:

      /// See EvaluationHandleBase.
      using Func = CODI_DD(T_Func, CODI_T(void (*)(T_InputStore const&, T_OutputStore&)));
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);      ///< See EvaluationHandleBase.
      using InputStore = CODI_DD(T_InputStore, std::vector<Type>);    ///< See EvaluationHandleBase.
      using OutputStore = CODI_DD(T_OutputStore, std::vector<Type>);  ///< See EvaluationHandleBase.

    protected:

      size_t m;  ///< Size of the output vector.
      size_t n;  ///< Size of the input vector.

      Func& func;  ///< The function object for the evaluations.

      InputStore x;   ///< Storage for the input arguments.
      OutputStore y;  ///< Storage for the output arguments.

      DummyVector dummyVector;      ///< Used if no output is required.
      DummyJacobian dummyJacobian;  ///< Used if no output is required.

    public:

      /// Constructor
      EvaluationHandleBase(Func& func, size_t m, size_t n)
          : m(m),
            n(n),
            func(func),
            x(constructVector<InputStore>(n)),
            y(constructVector<OutputStore>(m)),
            dummyVector(),
            dummyJacobian() {}

      /// Perform a primal evaluation with the inputs provided in locX and store the result in locY.
      template<typename VecX, typename VecY>
      void computePrimal(VecX const& locX, VecY& locY);

      /// Perform a Jacobian evaluation with the inputs provided in locX and store the result in jac and locY.
      template<typename VecX, typename Jac, typename VecY>
      void computeJacobian(VecX const& locX, Jac& jac, VecY& locY);

      /// Perform a Hessian evaluation with the inputs provided in locX and store the result in hes, jac and locY.
      template<typename VecX, typename Hes, typename VecY, typename Jac>
      void computeHessian(VecX const& locX, Hes& hes, VecY& locY, Jac& jac);

    protected:

      /// Set the primal values from the user provided vector into the CoDiPack ones.
      template<typename VecX>
      void setPrimalInputs(VecX const& locX);

      /// Store the primal values from the CoDiPack vector into the user vector.
      template<typename VecY>
      void getPrimalOutputs(VecY& locY);

      /// Helper for the evaluation of the function object with the internal input and output vector.
      void eval() {
        func(x, y);
      }
  };

  /// Implementation of EvaluationHandleBase for forward mode CoDiPack types.
  ///
  /// \copydetails EvaluationHandleBase
  template<typename T_Func, typename T_Type, typename T_InputStore = std::vector<T_Type>,
           typename T_OutputStore = std::vector<T_Type>>
  struct EvaluationHandleForward : public EvaluationHandleBase<T_Func, T_Type, T_InputStore, T_OutputStore> {
    public:

      /// See EvaluationHandleBase.
      using Func = CODI_DD(T_Func, CODI_T(void (*)(T_InputStore const&, T_OutputStore&)));
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);      ///< See EvaluationHandleBase.
      using InputStore = CODI_DD(T_InputStore, std::vector<Type>);    ///< See EvaluationHandleBase.
      using OutputStore = CODI_DD(T_OutputStore, std::vector<Type>);  ///< See EvaluationHandleBase.

      using Base = EvaluationHandleBase<Func, Type, InputStore, OutputStore>;  ///< Abbreviation for the base class.

      // Forward constructors of the base.
      using Base::EvaluationHandleBase;

      /// \copydoc codi::EvaluationHandleBase::setPrimalInputs
      template<typename VecX>
      void setPrimalInputs(VecX const& locX) {
        codiAssert(locX.size() <= this->x.size());
        for (size_t j = 0; j < locX.size(); j += 1) {
          this->x[j] = locX[j];
        }
      }

      /// \copydoc codi::EvaluationHandleBase::getPrimalOutputs
      template<typename VecY>
      void getPrimalOutputs(VecY& locY) {
        codiAssert(locY.size() <= this->y.size());
        for (size_t i = 0; i < locY.size(); i += 1) {
          locY[i] = RealTraits::getPassiveValue(this->y[i].getValue());
        }
      }

      /// \copydoc codi::EvaluationHandleBase::computePrimal
      template<typename VecX, typename VecY>
      void computePrimal(VecX const& locX, VecY& locY) {
        setPrimalInputs(locX);

        this->eval();

        getPrimalOutputs(locY);
      }

      /// \copydoc codi::EvaluationHandleBase::computeJacobian
      ///
      /// The vectorization is performed over the input vector. The function object is evaluated n/vecSize times.
      template<typename VecX, typename Jac, typename VecY>
      void computeJacobian(VecX const& locX, Jac& jac, VecY& locY) {
        setPrimalInputs(locX);

        JacobianConvertWrapper<Jac> wrapper(jac);

        using GradientTraits1st = GradientTraits::TraitsImplementation<typename Type::Gradient>;
        size_t constexpr VectorSizeFirstOrder = GradientTraits1st::dim;

        for (size_t j = 0; j < locX.size(); j += VectorSizeFirstOrder) {
          for (size_t vecPos = 0; vecPos < VectorSizeFirstOrder && j + vecPos < locX.size(); vecPos += 1) {
            GradientTraits1st::at(this->x[j + vecPos].gradient(), vecPos) = 1.0;
          }

          this->eval();

          if (0 == j) {
            getPrimalOutputs(locY);
          }

          for (size_t i = 0; i < this->y.size(); i += 1) {
            for (size_t vecPos = 0; vecPos < VectorSizeFirstOrder && j + vecPos < locX.size(); vecPos += 1) {
              wrapper(i, j + vecPos) = GradientTraits1st::at(this->y[i].gradient(), vecPos);
            }
          }

          for (size_t vecPos = 0; vecPos < VectorSizeFirstOrder && j + vecPos < locX.size(); vecPos += 1) {
            GradientTraits1st::at(this->x[j + vecPos].gradient(), vecPos) = 0.0;
          }
        }
      }

      /// \copydoc codi::EvaluationHandleBase::computeHessian
      ///
      /// The vectorization is performed twice over the input vector. This evaluates the Hessian in a blockwise fashion
      /// for all output values. The function object is evaluated n*n/(vecSize1 * vecSize2) times.
      template<typename VecX, typename Hes, typename VecY, typename Jac>
      void computeHessian(VecX const& locX, Hes& hes, VecY& locY, Jac& jac) {
        setPrimalInputs(locX);

        using GradientTraits1st = GradientTraits::TraitsImplementation<typename Type::Gradient>;
        size_t constexpr VectorSizeFirstOrder = GradientTraits1st::dim;

        using GradientTraits2nd = GradientTraits::TraitsImplementation<CODI_DD(typename Type::Real::Gradient, double)>;
        size_t constexpr VectorSizeSecondOrder = GradientTraits2nd::dim;

        for (size_t k = 0; k < locX.size(); k += VectorSizeFirstOrder) {
          // Set derivatives from k to k + vecSize_k.
          for (size_t vecPos = 0; vecPos < VectorSizeFirstOrder && k + vecPos < locX.size(); vecPos += 1) {
            GradientTraits1st::at(this->x[k + vecPos].gradient(), vecPos).value() = 1.0;
          }

          // The j = k init is no problem, it will evaluate slightly more elements around the diagonal.
          for (size_t j = k; j < locX.size(); j += VectorSizeSecondOrder) {
            // Set derivatives from j to j + vecSize_j.
            for (size_t vecPos = 0; vecPos < VectorSizeSecondOrder && j + vecPos < locX.size(); vecPos += 1) {
              GradientTraits2nd::at(this->x[j + vecPos].value().gradient(), vecPos) = 1.0;
            }

            this->eval();

            if (0 == j && 0 == k) {
              getPrimalOutputs(locY);
            }

            // Extract all Hessian values, this populates the Hessian from (j,k) to (j + vecSize_j, k + vecSize_k).
            for (size_t i = 0; i < this->y.size(); i += 1) {
              for (size_t vecPos1st = 0; vecPos1st < VectorSizeFirstOrder && k + vecPos1st < locX.size();
                   vecPos1st += 1) {
                for (size_t vecPos2nd = 0; vecPos2nd < VectorSizeSecondOrder && j + vecPos2nd < locX.size();
                     vecPos2nd += 1) {
                  auto& firstGrad = GradientTraits1st::at(this->y[i].gradient(), vecPos1st);
                  auto& secondGrad = GradientTraits2nd::at(firstGrad.gradient(), vecPos2nd);

                  hes(i, j + vecPos2nd, k + vecPos1st) = secondGrad;
                  hes(i, k + vecPos1st, j + vecPos2nd) = secondGrad;  // Symmetry
                }
              }

              if (k == 0) {
                for (size_t vecPos = 0; vecPos < VectorSizeSecondOrder && j + vecPos < locX.size(); vecPos += 1) {
                  jac(i, j + vecPos) = GradientTraits2nd::at(this->y[i].value().gradient(), vecPos);
                }
              }
            }

            // Reset the derivative seeding.
            for (size_t vecPos = 0; vecPos < VectorSizeSecondOrder && j + vecPos < locX.size(); vecPos += 1) {
              GradientTraits2nd::at(this->x[j + vecPos].value().gradient(), vecPos) = 0.0;
            }
          }

          // Reset the derivative seeding.
          for (size_t vecPos = 0; vecPos < VectorSizeFirstOrder && k + vecPos < locX.size(); vecPos += 1) {
            GradientTraits1st::at(this->x[k + vecPos].gradient(), vecPos).value() = 0.0;
          }
        }
      }
  };

  /// @brief Implementation for reverse mode CoDiPack types of EvaluationHandleBase.
  ///
  ///\copydetails EvaluationHandleBase
  template<typename T_Func, typename T_Type, typename T_InputStore = std::vector<T_Type>,
           typename T_OutputStore = std::vector<T_Type>>
  struct EvaluationHandleReverseBase : public EvaluationHandleBase<T_Func, T_Type, T_InputStore, T_OutputStore> {
    public:

      /// See EvaluationHandleBase.
      using Func = CODI_DD(T_Func, CODI_T(void (*)(T_InputStore const&, T_OutputStore&)));
      /// See EvaluationHandleBase.
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);
      using InputStore = CODI_DD(T_InputStore, std::vector<Type>);    ///< See EvaluationHandleBase.
      using OutputStore = CODI_DD(T_OutputStore, std::vector<Type>);  ///< See EvaluationHandleBase.

      using Base = EvaluationHandleBase<Func, Type, InputStore, OutputStore>;  ///< Abbreviation for the base class.

    protected:

      TapeHelper<Type> th;  ///< Manages the evaluations.

    public:

      /// Constructor
      EvaluationHandleReverseBase(Func& func, size_t m, size_t n) : Base(func, m, n), th() {}

      /// \copydoc codi::EvaluationHandleBase::setPrimalInputs
      template<typename VecX>
      void setPrimalInputs(VecX const& locX, bool reg) {
        codiAssert(locX.size() <= this->x.size());
        for (size_t j = 0; j < locX.size(); j += 1) {
          this->x[j] = locX[j];

          if (reg) {
            th.registerInput(this->x[j]);
          }
        }
      }

      /// \copydoc codi::EvaluationHandleBase::getPrimalOutputs
      template<typename VecY>
      void getPrimalOutputs(VecY& locY, bool reg) {
        codiAssert(locY.size() <= this->y.size());
        for (size_t i = 0; i < this->y.size(); i += 1) {
          if (reg) {
            th.registerOutput(this->y[i]);
          }

          locY[i] = RealTraits::getPassiveValue(this->y[i].getValue());
        }
      }

      /// \copydoc codi::EvaluationHandleBase::computePrimal
      template<typename VecX, typename VecY>
      void computePrimal(VecX const& locX, VecY& locY) {
        setPrimalInputs(locX, false);

        this->eval();

        getPrimalOutputs(locY, false);
      }

      /// \copydoc codi::EvaluationHandleBase::computeJacobian
      ///
      /// The best mode is selected for the evaluation of the Jacobian. If \f$ n < m \f$, the forward mode is used
      /// for the Jacobian evaluation and the function object is called \f$ n/vecSize \f$ times.
      /// If \f$ m < n \f$, the reverse mode is used and the function object is called \f$ m/vecSize \f$ times.
      template<typename VecX, typename Jac, typename VecY>
      void computeJacobian(VecX const& locX, Jac& jac, VecY& locY) {
        recordTape(locX, locY);

        th.evalJacobian(jac);
      }

      /// \copydoc codi::EvaluationHandleBase::computeHessian
      template<typename VecX, typename Hes, typename VecY, typename Jac>
      void computeHessian(VecX const& locX, Hes& hes, VecY& locY, Jac& jac);

    protected:

      /// Helper function that records a new tape.
      template<typename VecX, typename VecY>
      void recordTape(VecX const& locX, VecY& locY) {
        th.startRecording();
        setPrimalInputs(locX, true);

        this->eval();

        getPrimalOutputs(locY, true);
        th.stopRecording();
      }
  };

  /**
   * @brief Implementation of EvaluationHandleBase for primal value reverse mode CoDiPack types.
   *
   * This tape records the logic behind the function object once for each primal evaluation point. Afterwards,
   * only primal, forward and reverse tape evaluations are performed until the next primal evaluation point.
   *
   * Primal evaluations without derivative computations are not recorded.
   *
   * \copydetails EvaluationHandleBase
   */
  template<typename T_Func, typename T_Type, typename T_InputStore = std::vector<T_Type>,
           typename T_OutputStore = std::vector<T_Type>>
  struct EvaluationHandleReversePrimalValueTapes
      : public EvaluationHandleReverseBase<T_Func, T_Type, T_InputStore, T_OutputStore> {
    public:

      /// See EvaluationHandleBase.
      using Func = CODI_DD(T_Func, CODI_T(void (*)(T_InputStore const&, T_OutputStore&)));
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);      ///< See EvaluationHandleBase.
      using InputStore = CODI_DD(T_InputStore, std::vector<Type>);    ///< See EvaluationHandleBase.
      using OutputStore = CODI_DD(T_OutputStore, std::vector<Type>);  ///< See EvaluationHandleBase.

      /// Abbreviation for the base class.
      using Base = EvaluationHandleReverseBase<Func, Type, InputStore, OutputStore>;

      // Use constructors of the base class.
      using Base::EvaluationHandleReverseBase;

      /// \copydoc codi::EvaluationHandleBase::computeHessian
      ///
      /// For the primal value tape implementation, the tape is only recorded once and then evaluated multiple times.
      template<typename VecX, typename Hes, typename VecY, typename Jac>
      void computeHessian(VecX const& locX, Hes& hes, VecY& locY, Jac& jac) {
        this->recordTape(locX, locY);

        this->th.evalHessian(hes, jac);
      }
  };

  /**
   * @brief Implementation for Jacobian reverse mode CoDiPack types of EvaluationHandleBase.
   *
   * This tape re-records the logic behind the function object for every forward and reverse
   * tape evaluations. Primal evaluations are not recorded.
   *
   * \copydetails EvaluationHandleBase
   */
  template<typename T_Func, typename T_Type, typename T_InputStore = std::vector<T_Type>,
           typename T_OutputStore = std::vector<T_Type>>
  struct EvaluationHandleReverseJacobianTapes
      : public EvaluationHandleReverseBase<T_Func, T_Type, T_InputStore, T_OutputStore> {
    public:
      /// See EvaluationHandleBase.
      using Func = CODI_DD(T_Func, CODI_T(void (*)(T_InputStore const&, T_OutputStore&)));
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);      ///< See EvaluationHandleBase.
      using InputStore = CODI_DD(T_InputStore, std::vector<Type>);    ///< See EvaluationHandleBase.
      using OutputStore = CODI_DD(T_OutputStore, std::vector<Type>);  ///< See EvaluationHandleBase.

      /// Abbreviation for the base class.
      using Base = EvaluationHandleReverseBase<Func, Type, InputStore, OutputStore>;

      // Use constructors of the base class.
      using Base::EvaluationHandleReverseBase;

      /// \copydoc codi::EvaluationHandleBase::computeHessian
      ///
      /// For the Jacobian tape implementation, a new tape is recorded for every evaluation.
      template<typename VecX, typename Hes, typename VecY, typename Jac>
      void computeHessian(VecX const& locX, Hes& hes, VecY& locY, Jac& jac) {
        this->setPrimalInputs(locX, false);

        Algorithms<Type>::computeHessian(this->func, this->x, this->y, hes, jac);

        this->getPrimalOutputs(locY, false);
      }
  };

  /// See EvaluationHandleBase.
  template<typename T_Func, typename T_Type, typename T_InputStore = std::vector<T_Type>,
           typename T_OutputStore = std::vector<T_Type>, typename = void>
  struct EvaluationHandle : public EvaluationHandleBase<T_Func, T_Type, T_InputStore, T_OutputStore> {};

  /// See EvaluationHandleForward.
  template<typename T_Func, typename T_Type, typename T_InputStore, typename T_OutputStore>
  struct EvaluationHandle<T_Func, T_Type, T_InputStore, T_OutputStore,
                          TapeTraits::EnableIfForwardTape<typename T_Type::Tape>>
      : public EvaluationHandleForward<T_Func, T_Type, T_InputStore, T_OutputStore> {
      using EvaluationHandleForward<T_Func, T_Type, T_InputStore, T_OutputStore>::EvaluationHandleForward;
  };

  /// See EvaluationHandleReverseJacobianTapes.
  template<typename T_Func, typename T_Type, typename T_InputStore, typename T_OutputStore>
  struct EvaluationHandle<T_Func, T_Type, T_InputStore, T_OutputStore,
                          TapeTraits::EnableIfJacobianTape<typename T_Type::Tape>>
      : public EvaluationHandleReverseJacobianTapes<T_Func, T_Type, T_InputStore, T_OutputStore> {
      using EvaluationHandleReverseJacobianTapes<T_Func, T_Type, T_InputStore,
                                                 T_OutputStore>::EvaluationHandleReverseJacobianTapes;
  };

  /// See EvaluationHandleReversePrimalValueTapes.
  template<typename T_Func, typename T_Type, typename T_InputStore, typename T_OutputStore>
  struct EvaluationHandle<T_Func, T_Type, T_InputStore, T_OutputStore,
                          TapeTraits::EnableIfPrimalValueTape<typename T_Type::Tape>>
      : public EvaluationHandleReversePrimalValueTapes<T_Func, T_Type, T_InputStore, T_OutputStore> {
      using EvaluationHandleReversePrimalValueTapes<T_Func, T_Type, T_InputStore,
                                                    T_OutputStore>::EvaluationHandleReversePrimalValueTapes;
  };

  /**
   * @brief Evaluate the primal, Jacobian and Hessian of function objects.
   *
   * This helper provides the means to easily evaluate derivatives of arbitrary function objects. These function
   * objects can be regular functions, lambda functions or structures where operator() is implemented.
   *
   * The nomenclature and mathematical definitions for the function, the Jacobian, and the Hessian can be found in the
   * \ref sec_namingConventions documentation. Function arguments in this class follow the same naming scheme.
   *
   * The algorithms will call the function objects with the vector of inputs and with the vector of outputs as
   * arguments. The function object has to resemble the interface defined by FunctorInterface.
   * An example function definition is:
   * \code{.cpp}
   *  void func1(std::vector<ADType> const& x, std::vector<ADType>& y) { ... }
   *  // or
   *  void func1(std::array<ADType, n> const& x, std::array<ADType, m>& y) { ... }
   * \endcode
   *
   * x is the vector of input values and y is the vector of output values. ADType is the chosen CoDiPack type for the
   * function. For most users, this definition will be enough. For more general examples please go to section
   * \ref AdvFuncObjDef.
   *
   * The CoDiPack type can be any ActiveType type from CoDiPack, for example all types that are defined in codi.hpp. The
   * evaluation helper provides the default CoDiPack type definitions EvaluationHelper::JacobianComputationType and
   * EvaluationHelper::HessianComputationType. These two use the forward mode of algorithmic differentiation and are
   * more appropriate if \f$ m \f$ and \f$ n \f$ are small. They can also be used if \f$ n \f$ is smaller than
   * \f$ m \f$. For cases in in which the dimensions are larger and \f$ m \f$ is smaller than \f$ n \f$, the types
   * codi::JacobianComputationType and codi::HessianComputationType can be used. They use the reverse AD mode for the
   * computation. Please refer to \ref AdvFuncObjDef to see how these types can be used.
   *
   * The most simple example of using the EvaluationHelper is:
   * \snippet examples/evaluationHelper_minimal.cpp EvaluationHelper minimal
   *
   * Since we want to evaluate the Hessian, we use the Hessian type of the EvaluationHelper. The function is defined
   * with this type and the standard vector classes. In the main function, we create the vector on which we want to call
   * the function and then use the helper to create the storage for the Jacobian and Hessian. With
   * \f$ jac(j,i) \f$  the values can be accessed. For the Hessian the values can be accessed with \f$ hes(j,i,k) \f$
   * where \f$ j \f$ is the output dimension and \f$ i \f$ as well as \f$ k \f$ are the input dimensions.
   *
   * The evaluation helper class provides all combinations of evaluation choices, that is: evalPrimal(), evalJacobian(),
   * evalHessian(), evalPrimalAndJacobian(), evalPrimalAndHessian(), evalJacobianAndHessian() and
   * evalPrimalAndJacobianAndHessian(). Each of these functions uses the default CoDiPack types in the evaluation
   * helper. In the cases where the primal is not stored, the user has to provide the number of output values manually.
   *
   * If the EvaluationHelper is used to evaluate the same function several times, a higher performance can be achieved
   * if a handle for the evaluation is created up front and then used several times. The above example with the handle
   * creation would look like this:
   * \snippet examples/evaluationHelper_handle.cpp EvaluationHelper changes
   *
   * The evaluation logic nearly stayed the same, but instead of providing the function to the evaluation routine, we
   * create a handle up front and then use this handle in the evalHandle method. Each of the above mentioned eval
   * routines has a corresponding evalHandle method.
   *
   * Each of the create methods has similar create..Fixed method which uses the std::array type instead of the
   * std::vector type for the data management. These methods can be used if the size is known at compile time. An
   * example with these methods would be:
   * \snippet examples/evaluationHelper_fixed.cpp EvaluationHelper fixed
   *
   * Until now, the default definition for the used CoDiPack types have been used. In order to use an arbitrary CoDiPack
   * type, the createHandle(), createHandleFixed(), and createHandleFull() methods can be used. The first one uses
   * std::vectors for the storage, the second method std::array, and in the third the user can provide the storage class
   * as a template parameter. The use case for the createHandle() method would look like this:
   * \code{.cpp}
   *  auto handle = eh.createHandle<codi::RealReverse>(func, 4, 2);
   * \endcode
   *
   * \section AdvFuncObjDef Advanced function object definitions
   * The function object can also have a template argument for the evaluation type, e.g.:
   * \code{.cpp}
   *  struct Func {
   *    template<typename T>
   *    void operator()(std::vector<T> const& x, std::vector<T>& y);
   *  };
   * \endcode
   *
   * There is also no need to specify std::vector as the array class, e.g.:
   * \code{.cpp}
   *  struct Func {
   *    template<typename InVec, typename OutVec>
   *    void operator()(InVec const& x, OutVec& y);
   *  };
   * \endcode
   *
   */
  struct EvaluationHelper {
    public:

      /// Function object syntax for all Func template arguments.
      /// @tparam VecIn   User defined (default: std::vector).
      /// @tparam VecOut  User defined (default: std::vector).
      template<typename VecIn, typename VecOut>
      using FunctorInterface = void (*)(VecIn const& x, VecOut& y);

      /// The default type used for first order derivative computation. It is defined as forward vector AD mode of
      /// size 4.
      using JacobianComputationType = RealForwardVec<4>;

      /// The default type used for second order derivative computation. It is defined as forward vector over forward
      /// vector AD mode of size 4 and 4. */
      using HessianComputationType = RealForwardGen<RealForwardVec<4>, Direction<RealForwardVec<4>, 4>>;

      /// Type for the default handle for first order derivative computations with a variable vector size.
      /// @tparam Func  See FunctorInterface.
      template<typename Func>
      using DefaultHandle = EvaluationHandleForward<Func, JacobianComputationType>;

      /// Type for the default handle for second order derivative computations with a variable vector size.
      /// @tparam Func  See FunctorInterface.
      template<typename Func>
      using DefaultHandle2nd = EvaluationHandleForward<Func, HessianComputationType>;

      /// Type for the default handle for first order derivative computations with a compile time vector size.
      /// @tparam Func  See FunctorInterface
      /// @tparam    m  The size of the output vector.
      /// @tparam    n  The size of the input vector.
      template<typename Func, size_t m, size_t n>
      using DefaultHandleFixed =
          EvaluationHandleForward<Func, JacobianComputationType, std::array<JacobianComputationType, n>,
                                  std::array<JacobianComputationType, m>>;

      /// Type for the default handle for second order derivative computations with a compile time vector size.
      /// @tparam Func  See FunctorInterface
      /// @tparam    m  The size of the output vector.
      /// @tparam    n  The size of the input vector.
      template<typename Func, size_t m, size_t n>
      using DefaultHandleFixed2nd =
          EvaluationHandleForward<Func, HessianComputationType, std::array<HessianComputationType, n>,
                                  std::array<HessianComputationType, m>>;

      /**
       * @brief Helper function for the creation of a default first order evaluation handle with a variable vector size.
       *
       * @param[in] func  The function object for the evaluation (see FunctorInterface).
       * @param[in]    m  The size of the output vector.
       * @param[in]    n  The size of the input vector.
       *
       * @tparam Func  See FunctorInterface.
       */
      template<typename Func>
      static CODI_INLINE DefaultHandle<Func> createHandleDefault(Func& func, size_t m, size_t n) {
        return DefaultHandle<Func>(func, m, n);
      }

      /**
       * @brief Helper function for the creation of a default first order evaluation handle with a compile time vector
       * size.
       *
       * @param[in] func  The function object for the evaluation (see FunctorInterface).
       *
       * @tparam    m  The size of the output vector.
       * @tparam    n  The size of the input vector.
       * @tparam Func  See FunctorInterface.
       */
      template<size_t m, size_t n, typename Func>
      static CODI_INLINE DefaultHandleFixed<Func, m, n> createHandleDefaultFixed(Func& func) {
        return DefaultHandleFixed<Func, m, n>(func, m, n);
      }

      /**
       * @brief Helper function for the creation of a default second order evaluation handle with a variable vector
       *        size.
       *
       * @param[in] func  The function object for the evaluation (see FunctorInterface).
       * @param[in]    m  The size of the output vector.
       * @param[in]    n  The size of the input vector.
       *
       * @tparam Func  See FunctorInterface.
       */
      template<typename Func>
      static CODI_INLINE DefaultHandle2nd<Func> createHandleDefault2nd(Func& func, size_t m, size_t n) {
        return DefaultHandle2nd<Func>(func, m, n);
      }

      /**
       * @brief Helper function for the creation of a default second order evaluation handle with a compile time vector
       * size.
       *
       * @param[in] func  The function object for the evaluation (see FunctorInterface).
       *
       * @tparam    m  The size of the output vector.
       * @tparam    n  The size of the input vector.
       * @tparam Func  See FunctorInterface.
       */
      template<size_t m, size_t n, typename Func>
      static CODI_INLINE DefaultHandleFixed2nd<Func, m, n> createHandleDefaultFixed2nd(Func& func) {
        return DefaultHandleFixed2nd<Func, m, n>(func, m, n);
      }

      /**
       * @brief Helper function for the creation of an evaluation handle with the specified CoDiPack type and a variable
       *        vector size.
       *
       * The CoDiPack type can be an arbitrary one:
       * \code{.cpp}
       *  auto handle = codi::EvaluationHelper::createHandle<codi::RealReverse>(func, m, n);
       * \endcode
       *
       * @param[in] func  The function object for the evaluation (see FunctorInterface).
       * @param[in]    m  The size of the output vector.
       * @param[in]    n  The size of the input vector.
       *
       * @tparam CoDiType  An arbitrary CoDiPack type based on ActiveType. All definitions in codi.hpp are supported.
       *                   For user developed tapes, EvaluationHandle has to be specialized.
       * @tparam     Func  See FunctorInterface.
       */
      template<typename Type, typename Func>
      static CODI_INLINE EvaluationHandle<Func, Type> createHandle(Func& func, size_t m, size_t n) {
        return EvaluationHandle<Func, Type>(func, m, n);
      }

      /**
       * @brief Helper function for the creation of an evaluation handle with the specified CoDiPack type and a compile
       * time vector size.
       *
       * The CoDiPack type can be an arbitrary one:
       * \code{.cpp}
       *  auto handle = codi::EvaluationHelper::createHandleFixed<codi::RealReverse, m, n>(func);
       * \endcode
       *
       * @param[in] func  The function object for the evaluation (see FunctorInterface).
       *
       * @tparam CoDiType  An arbitrary CoDiPack type based on ActiveType. All definitions in codi.hpp are supported.
       *                   For user developed tapes, EvaluationHandle has to be specialized.
       * @tparam     Func  See FunctorInterface.
       * @tparam        m  The size of the output vector.
       * @tparam        n  The size of the input vector.
       */
      template<typename Type, size_t m, size_t n, typename Func>
      static CODI_INLINE EvaluationHandle<Func, Type, std::array<Type, n>, std::array<Type, m>> createHandleFixed(
          Func& func) {
        return EvaluationHandle<Func, Type, std::array<Type, n>, std::array<Type, m>>(func, m, n);
      }

      /**
       * @brief Helper function for the creation of an evaluation handle with the specified CoDiPack type and storage
       *        types.
       *
       * The CoDiPack type can be an arbitrary one and the storage types must use it for their elements.
       *
       * \code{.cpp}
       *   auto handle = codi::EvaluationHelper::createHandleFull<codi::RealReverse,
       *                                                          std::vector<codi::RealReverse>,
       *                                                          std::array<<codi::RealReverse, m>>(func, m, n);
       * \endcode
       *
       * @param[in] func  The function object for the evaluation (see FunctorInterface).
       * @param[in]    m  The size of the output vector.
       * @param[in]    n  The size of the input vector.
       *
       * @tparam     CoDiType  An arbitrary CoDiPack type based on ActiveType. All definitions in codi.hpp are
       *                       supported. For user developed tapes, EvaluationHandle has to be specialized.
       * @tparam        Func  See FunctorInterface
       * @tparam  InputStore  The storage type for vectors of input variables.
       * @tparam OutputStore  The storage type for vectors of output variables.
       */
      template<typename Type, typename InputStore, typename OutputStore, typename Func>
      static CODI_INLINE EvaluationHandle<Func, Type, InputStore, OutputStore> createHandleFull(Func& func, size_t m,
                                                                                                size_t n) {
        return EvaluationHandle<Func, Type, InputStore, OutputStore>(func, m, n);
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
      static CODI_INLINE Jacobian<T> createJacobian(size_t m, size_t n) {
        return Jacobian<T>(m, n);
      }

      /**
       * @brief Create a Jacobian with a compile time size.
       *
       * @tparam T  The storage type of the Jacobian.
       * @tparam m  The size of the output vector.
       * @tparam n  The size of the input vector.
       */
      template<size_t m, size_t n, typename T = double>
      static CODI_INLINE Jacobian<T, std::array<T, m * n>> createJacobianFixed() {
        return Jacobian<T, std::array<T, m * n>>(m, n);
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
      static CODI_INLINE Hessian<T> createHessian(size_t m, size_t n) {
        return Hessian<T>(m, n);
      }

      /**
       * @brief Create a Hessian with a compile time size.
       *
       * @tparam T  The storage type of the Hessian.
       * @tparam        m  The size of the output vector.
       * @tparam        n  The size of the input vector.
       */
      template<size_t m, size_t n, typename T = double>
      static CODI_INLINE Hessian<T, std::array<T, m * n * n>> createHessianFixed() {
        return Hessian<T, std::array<T, m * n * n>>(m, n);
      }

      /**
       * @brief Perform a primal evaluation of the function object with the default first order type.
       *
       * @param[in] func  The function object for the evaluation (see FunctorInterface).
       * @param[in]    x  The vector with the primal values where the function object is evaluated.
       * @param[out]   y  The vector for the result of the primal function evaluation. The vector must have the
       *                  correct size allocated.
       *
       * @tparam  Func  See FunctorInterface.
       * @tparam  VecX  The vector type for the input values. Element type is e.g. double.
       * @tparam  VecY  The vector type for the output values. Element type is e.g. double.
       */
      template<typename Func, typename VecX, typename VecY>
      static CODI_INLINE void evalPrimal(Func& func, VecX const& x, VecY& y) {
        auto h = createHandleDefault(func, y.size(), x.size());
        evalHandlePrimal(h, x, y);
      }

      /**
       * @brief Compute the Jacobian of the function object.
       *
       * @param[in]  func  The function object for the evaluation (see FunctorInterface).
       * @param[in]     x  The vector with the primal values where the function object is evaluated.
       * @param[in] ySize  The size of the output vector.
       * @param[out]  jac  The Jacobian in which the values are stored.
       *
       * @tparam  Func  See FunctorInterface.
       * @tparam  VecX  The vector type for the input values. Element type is e.g. double.
       * @tparam   Jac  The storage type of the Jacobian. Element type is e.g. double.
       */
      template<typename Func, typename VecX, typename Jac>
      static CODI_INLINE void evalJacobian(Func& func, VecX const& x, size_t const ySize, Jac& jac) {
        auto h = createHandleDefault(func, ySize, x.size());
        evalHandleJacobian(h, x, jac);
      }

      /**
       * @brief Compute the Hessian of the function object.
       *
       * @param[in]  func  The function object for the evaluation (see FunctorInterface).
       * @param[in]     x  The vector with the primal values where the function object is evaluated.
       * @param[in] ySize  The size of the output vector.
       * @param[out]  hes  The Hessian in which the values are stored.
       *
       * @tparam  Func  See FunctorInterface.
       * @tparam  VecX  The vector type for the input values. Element type is e.g. double.
       * @tparam   Hes  The storage type of the Hessian. Element type is e.g. double.
       */
      template<typename Func, typename VecX, typename Hes>
      static CODI_INLINE void evalHessian(Func& func, VecX const& x, size_t const ySize, Hes& hes) {
        auto h = createHandleDefault2nd(func, ySize, x.size());
        evalHandleHessian(h, x, hes);
      }

      /**
       * @brief Compute the primal result and Jacobian of the function object.
       *
       * @param[in]  func  The function object for the evaluation (see FunctorInterface).
       * @param[in]     x  The vector with the primal values where the function object is evaluated.
       * @param[out]    y  The vector for the result of the primal function evaluation. The vector must have the
       *                   correct size allocated.
       * @param[out]  jac  The Jacobian in which the values are stored.
       *
       * @tparam  Func  See FunctorInterface.
       * @tparam  VecX  The vector type for the input values. Element type is e.g. double.
       * @tparam  VecY  The vector type for the output values. Element type is e.g. double.
       * @tparam   Jac  The storage type of the Jacobian. Element type is e.g. double.
       */
      template<typename Func, typename VecX, typename VecY, typename Jac>
      static CODI_INLINE void evalPrimalAndJacobian(Func& func, VecX const& x, VecY& y, Jac& jac) {
        auto h = createHandleDefault(func, y.size(), x.size());
        evalHandlePrimalAndJacobian(h, x, y, jac);
      }

      /**
       * @brief Compute the primal result and Hessian of the function object.
       *
       * @param[in]  func  The function object for the evaluation (see FunctorInterface).
       * @param[in]     x  The vector with the primal values where the function object is evaluated.
       * @param[out]    y  The vector for the result of the primal function evaluation. The vector must have the
       *                   correct size allocated.
       * @param[out]  hes  The Hessian in which the values are stored.
       *
       * @tparam  Func  See FunctorInterface.
       * @tparam  VecX  The vector type for the input values. Element type is e.g. double.
       * @tparam  VecY  The vector type for the output values. Element type is e.g. double.
       * @tparam   Hes  The storage type of the Hessian. Element type is e.g. double.
       */
      template<typename Func, typename VecX, typename VecY, typename Hes>
      static CODI_INLINE void evalPrimalAndHessian(Func& func, VecX const& x, VecY& y, Hes& hes) {
        auto h = createHandleDefault2nd(func, y.size(), x.size());
        evalHandlePrimalAndHessian(h, x, y, hes);
      }

      /**
       * @brief Compute the primal result, Jacobian, and Hessian of the function object.
       *
       * @param[in]  func  The function object for the evaluation (see FunctorInterface).
       * @param[in]     x  The vector with the primal values where the function object is evaluated.
       * @param[out]    y  The vector for the result of the primal function evaluation. The vector must have the
       *                   correct size allocated.
       * @param[out]  jac  The Jacobian in which the values are stored.
       * @param[out]  hes  The Hessian in which the values are stored.
       *
       * @tparam  Func  See FunctorInterface.
       * @tparam  VecX  The vector type for the input values. Element type is e.g. double.
       * @tparam  VecY  The vector type for the output values. Element type is e.g. double.
       * @tparam   Jac  The storage type of the Jacobian. Element type is e.g. double.
       * @tparam   Hes  The storage type of the Hessian. Element type is e.g. double.
       */
      template<typename Func, typename VecX, typename VecY, typename Jac, typename Hes>
      static CODI_INLINE void evalPrimalAndJacobianAndHessian(Func& func, VecX const& x, VecY& y, Jac& jac, Hes& hes) {
        auto h = createHandleDefault2nd(func, y.size(), x.size());
        evalHandlePrimalAndJacobianAndHessian(h, x, y, jac, hes);
      }

      /**
       * @brief Compute the Jacobian and Hessian of the function object.
       *
       * @param[in]  func  The function object for the evaluation (see FunctorInterface).
       * @param[in]     x  The vector with the primal values where the function object is evaluated.
       * @param[in] ySize  The size of the output vector.
       * @param[out]  jac  The Jacobian in which the values are stored.
       * @param[out]  hes  The Hessian in which the values are stored.
       *
       * @tparam  Func  See FunctorInterface
       * @tparam  VecX  The vector type for the input values. Element type is e.g. double.
       * @tparam   Jac  The storage type of the Jacobian. Element type is e.g. double.
       * @tparam   Hes  The storage type of the Hessian. Element type is e.g. double.
       */
      template<typename Func, typename VecX, typename Jac, typename Hes>
      static CODI_INLINE void evalJacobianAndHessian(Func& func, VecX const& x, size_t ySize, Jac& jac, Hes& hes) {
        auto h = createHandleDefault2nd(func, ySize, x.size());
        evalHandleJacobianAndHessian(h, x, jac, hes);
      }

      /**
       * @brief Perform a primal evaluation of the function object stored in the handle.
       *
       * @param[in] handle  The handle with all data for the evaluation.
       * @param[in]     x  The vector with the primal values where the function object is evaluated.
       * @param[out]    y  The vector for the result of the primal function evaluation. The vector must have the
       *                   correct size allocated.
       *
       * @tparam  Handle  The handle type for the data storage and the evaluation.
       * @tparam    VecX  The vector type for the input values. Element type is e.g. double.
       * @tparam    VecY  The vector type for the output values. Element type is e.g. double.
       */
      template<typename Handle, typename VecX, typename VecY>
      static CODI_INLINE void evalHandlePrimal(Handle& handle, VecX const& x, VecY& y) {
        handle.computePrimal(x, y);
      }

      /**
       * @brief Compute the Jacobian of the function object stored in the handle.
       *
       * @param[in] handle  The handle with all data for the evaluation.
       * @param[in]     x  The vector with the primal values where the function object is evaluated.
       * @param[out]  jac  The Jacobian in which the values are stored.
       *
       * @tparam  Handle  The handle type for the data storage and the evaluation.
       * @tparam    VecX  The vector type for the input values. Element type is e.g. double.
       * @tparam     Jac  The storage type of the Jacobian. Element type is e.g. double.
       */
      template<typename Handle, typename VecX, typename Jac>
      static CODI_INLINE void evalHandleJacobian(Handle& handle, VecX const& x, Jac& jac) {
        DummyVector dv;
        handle.computeJacobian(x, jac, dv);
      }

      /**
       * @brief Compute the Hessian of the function object stored in the handle.
       *
       * @param[in] handle  The handle with all data for the evaluation.
       * @param[in]     x  The vector with the primal values where the function object is evaluated.
       * @param[out]  hes  The Hessian in which the values are stored.
       *
       * @tparam  Handle  The handle type for the data storage and the evaluation.
       * @tparam    VecX  The vector type for the input values. Element type is e.g. double.
       * @tparam     Hes  The storage type of the Hessian. Element type is e.g. double.
       */
      template<typename Handle, typename VecX, typename Hes>
      static CODI_INLINE void evalHandleHessian(Handle& handle, VecX const& x, Hes& hes) {
        DummyVector dv;
        DummyJacobian dj;
        handle.computeHessian(x, hes, dv, dj);
      }

      /**
       * @brief Compute the primal result and Jacobian of the function object stored in the handle.
       *
       * @param[in] handle  The handle with all data for the evaluation.
       * @param[in]     x  The vector with the primal values where the function object is evaluated.
       * @param[out]    y  The vector for the result of the primal function evaluation. The vector must have the
       *                   correct size allocated.
       * @param[out]  jac  The Jacobian in which the values are stored.
       *
       * @tparam  Handle  The handle type for the data storage and the evaluation.
       * @tparam    VecX  The vector type for the input values. Element type is e.g. double.
       * @tparam    VecY  The vector type for the output values. Element type is e.g. double.
       * @tparam     Jac  The storage type of the Jacobian. Element type is e.g. double.
       */
      template<typename Handle, typename VecX, typename VecY, typename Jac>
      static CODI_INLINE void evalHandlePrimalAndJacobian(Handle& handle, VecX const& x, VecY& y, Jac& jac) {
        handle.computeJacobian(x, jac, y);
      }

      /**
       * @brief Compute the primal result and Hessian of the function object.
       *
       * @param[in] handle  The handle with all data for the evaluation.
       * @param[in]     x  The vector with the primal values where the function object is evaluated.
       * @param[out]    y  The vector for the result of the primal function evaluation. The vector must have the
       *                   correct size allocated.
       * @param[out]  hes  The Hessian in which the values are stored.
       *
       * @tparam  Handle  The handle type for the data storage and the evaluation.
       * @tparam    VecX  The vector type for the input values. Element type is e.g. double.
       * @tparam    VecY  The vector type for the output values. Element type is e.g. double.
       * @tparam     Hes  The storage type of the Hessian. Element type is e.g. double.
       */
      template<typename Handle, typename VecX, typename VecY, typename Hes>
      static CODI_INLINE void evalHandlePrimalAndHessian(Handle& handle, VecX const& x, VecY& y, Hes& hes) {
        DummyJacobian dj;
        handle.computeHessian(x, hes, y, dj);
      }

      /**
       * @brief Compute the primal result, Jacobian, and Hessian of the function object.
       *
       * @param[in] handle  The handle with all data for the evaluation.
       * @param[in]     x  The vector with the primal values where the function object is evaluated.
       * @param[out]    y  The vector for the result of the primal function evaluation. The vector must have the
       *                   correct size allocated.
       * @param[out]  jac  The Jacobian in which the values are stored.
       * @param[out]  hes  The Hessian in which the values are stored.
       *
       * @tparam  Handle  The handle type for the data storage and the evaluation.
       * @tparam    VecX  The vector type for the input values. Element type is e.g. double.
       * @tparam    VecY  The vector type for the output values. Element type is e.g. double.
       * @tparam     Jac  The storage type of the Jacobian. Element type is e.g. double.
       * @tparam     Hes  The storage type of the Hessian. Element type is e.g. double.
       */
      template<typename Handle, typename VecX, typename VecY, typename Jac, typename Hes>
      static CODI_INLINE void evalHandlePrimalAndJacobianAndHessian(Handle& handle, VecX const& x, VecY& y, Jac& jac,
                                                                    Hes& hes) {
        handle.computeHessian(x, hes, y, jac);
      }

      /**
       * @brief Compute the Hessian of the evaluation procedure in the function object.
       *
       * In this method the Jacobian is also stored.
       *
       * @param[in] handle  The handle with all data for the evaluation.
       * @param[in]     x  The vector with the primal values where the function object is evaluated.
       * @param[out]  jac  The Jacobian in which the values are stored.
       * @param[out]  hes  The Hessian in which the values are stored.
       *
       * @tparam  Handle  The handle type for the data storage and the evaluation.
       * @tparam    VecX  The vector type for the input values. Element type is e.g. double.
       * @tparam     Jac  The storage type of the Jacobian. Element type is e.g. double.
       * @tparam     Hes  The storage type of the Hessian. Element type is e.g. double.
       */
      template<typename Handle, typename VecX, typename Jac, typename Hes>
      static CODI_INLINE void evalHandleJacobianAndHessian(Handle& handle, VecX const& x, Jac& jac, Hes& hes) {
        DummyVector dv;
        handle.computeHessian(x, hes, dv, jac);
      }
  };

}
