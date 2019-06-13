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

  /*
   * Function y = f(x), x \in \R^n, y \in \R^m
   * Jacobian \frac{\d f}{\d x} \in \R^{m \cross n}
   * Hessian \frac{\d^2 f}{\d^2 x} \in \R^{m \cross n \cross n}
   */

  template <typename Func, typename CoDiType,
            template<typename> class InputVectorType,
            template<typename> class OutputVectorType>
  struct EvaluationHandlerBase {
      using InputVector = InputVectorType<CoDiType>;
      using OutputVector = OutputVectorType<CoDiType>;

    protected:
      size_t m;
      size_t n;

      Func func;

      VectorStorage<InputVector> x;
      VectorStorage<OutputVector> y;

      DummyVector dummyVector;
      DummyJacobian dummyJacobian;

    public:

      EvaluationHandlerBase(Func func, size_t m, size_t n) : m(m), n(n), func(func), x(n), y(m),
          dummyVector(), dummyJacobian() {}

    protected:

      void eval() {
        func(x.vec, y.vec);
      }
  };

  template <typename Func, typename CoDiType,
            template<typename> class InputVectorType,
            template<typename> class OutputVectorType>
  struct ForwardHandle : public EvaluationHandlerBase<Func, CoDiType, InputVectorType, OutputVectorType> {

      ForwardHandle(Func func, size_t m, size_t n) : EvaluationHandlerBase<Func, CoDiType, InputVectorType, OutputVectorType>(func, m, n) {}

      template<typename VecX>
      void setAllPrimals(const VecX& locX) {
        codiAssert(locX.size() <= this->x.size());
        for(size_t j = 0; j < locX.size(); j += 1) {
          this->x[j] = locX[j];
        }
      }

      template<typename VecY>
      void getAllPrimals(VecY& locY) {
        codiAssert(locY.size() <= this->y.size());
        for(size_t i = 0; i < locY.size(); i += 1) {
          locY[i] = this->y[i];
        }
      }

      template<typename VecX, typename VecY>
      void computePrimal(const VecX& locX, VecY& locY) {
        setAllPrimals(locX);

        this->eval();

        getAllPrimals(locY);
      }

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

          // The j = k init is no prolbem, it will evaluated slightly more elements around the diagonal
          for(size_t j = k; j < locX.size(); j+= VectorSizeSecondOrder) {

            // Set derivatives from j to j + vecSize_j
            for(size_t vecPos = 0; vecPos < VectorSizeSecondOrder && j + vecPos < locX.size(); vecPos += 1) {
              GradientTraits2nd::at(this->x[j + vecPos].value().gradient(), vecPos) = 1.0;
            }

            this->eval();

            if(0 == j && 0 == k) {
              getAllPrimals(locY);
            }

            // Extract all hessian values, this pobulates the hessian from (j,k) to (j + vecSize_j, k + vecSize_k).
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

            // Reset the deriative seeding
            for(size_t vecPos = 0; vecPos < VectorSizeSecondOrder && j + vecPos < locX.size(); vecPos += 1) {
              GradientTraits2nd::at(this->x[j + vecPos].value().gradient(), vecPos) = 0.0;
            }
          }

          // Reset the deriative seeding
          for(size_t vecPos = 0; vecPos < VectorSizeFirstOrder && k + vecPos < locX.size(); vecPos += 1) {
            GradientTraits1st::at(this->x[k + vecPos].gradient(), vecPos).value() = 0.0;
          }
        }
      }
  };

  template <typename Func, typename CoDiType,
            template<typename> class InputVectorType,
            template<typename> class OutputVectorType>
  struct ReverseHandle : public EvaluationHandlerBase<Func, CoDiType, InputVectorType, OutputVectorType> {

      TapeHelper<CoDiType> th;

      ReverseHandle(Func func, size_t m, size_t n) :
        EvaluationHandlerBase<Func, CoDiType, InputVectorType, OutputVectorType>(func, m, n),
        th()
      {}

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

      template<typename VecX, typename VecY>
      void computePrimal(const VecX& locX, VecY& locY) {
        setAllPrimals(locX, false);

        this->eval();

        getAllPrimals(locY, false);
      }

      template<typename VecX, typename Jac, typename VecY>
      void computeJacobian(const VecX& locX, Jac& jac, VecY& locY) {
        recordTape(locX, locY);

        th.evalJacobian(jac);
      }

    private:

      template<typename VecX, typename VecY>
      void recordTape(const VecX& locX, VecY& locY) {
        th.startRecording();
        setAllPrimals(locX, true);

        this->eval();

        getAllPrimals(locY, true);
        th.stopRecording();
      }

      template<typename VecX, typename Hes, typename VecY, typename Jac>
      void computeHessian(const VecX& locX, Hes& hes, VecY& locY, Jac& jac) {
        // TODO
      }
  };

  template <typename Func, typename CoDiType,
            template<typename> class InputVectorType,
            template<typename> class OutputVectorType,
            typename = void>
  struct HandleSwitch;

  template <typename Func, typename CoDiType,
            template<typename> class InputVectorType,
            template<typename> class OutputVectorType>
  struct HandleSwitch<Func, CoDiType, InputVectorType, OutputVectorType, enableIfForwardTape<typename CoDiType::TapeType>> :
      public ForwardHandle<Func, CoDiType, InputVectorType, OutputVectorType> {

      HandleSwitch(Func func, size_t m, size_t n) : ForwardHandle<Func, CoDiType, InputVectorType, OutputVectorType>(func, m, n) {}
  };

  template <typename Func, typename CoDiType,
            template<typename> class InputVectorType,
            template<typename> class OutputVectorType>
  struct HandleSwitch<Func, CoDiType, InputVectorType, OutputVectorType, enableIfReverseTape<typename CoDiType::TapeType>> :
      public ReverseHandle<Func, CoDiType, InputVectorType, OutputVectorType> {

      HandleSwitch(Func func, size_t m, size_t n) : ReverseHandle<Func, CoDiType, InputVectorType, OutputVectorType>(func, m, n) {}
  };


  struct EvaluationHelper {

      template<size_t s>
      struct StdArrayAdapter {
          template <typename T>
          using Type = std::array<T, s>;
      };

      template <typename T>
      using StdVectorAdapter = std::vector<T>;


      template <typename Func>
      using DefaultHandle = ForwardHandle<Func, RealForwardVec<4>, StdVectorAdapter, StdVectorAdapter>;

      template <typename Func>
      using DefaultHandle2nd = ForwardHandle<Func, RealForwardGen<RealForwardVec<4>, Direction<RealForwardVec<4>, 4>>, StdVectorAdapter, StdVectorAdapter>;

      template <typename Func, size_t m, size_t n>
      using DefaultHandleFixed = ForwardHandle<Func, RealForwardVec<4>, StdArrayAdapter<n>::template Type, StdArrayAdapter<m>::template Type>;

      template <typename Func, size_t m, size_t n>
      using DefaultHandleFixed2nd = ForwardHandle<Func, RealForwardGen<RealForwardVec<4>, Direction<RealForwardVec<4>, 4>> , StdArrayAdapter<n>::template Type, StdArrayAdapter<m>::template Type>;


      template<typename Func>
      static CODI_INLINE DefaultHandle<Func> createDefault(const Func& func, size_t m, size_t n) {
        return DefaultHandle<Func>(func, m, n);
      }

      template<typename Func, size_t m, size_t n>
      static CODI_INLINE DefaultHandleFixed<Func, m, n> createDefaultFixed(const Func& func) {
        return DefaultHandleFixed<Func, m, n>(func, m, n);
      }

      template<typename Func>
      static CODI_INLINE DefaultHandle2nd<Func> createDefault2nd(const Func& func, size_t m, size_t n) {
        return DefaultHandle2nd<Func>(func, m, n);
      }

      template<typename Func, size_t m, size_t n>
      static CODI_INLINE DefaultHandleFixed2nd<Func, m, n> createDefaultFixed2nd(const Func& func) {
        return DefaultHandleFixed2nd<Func, m, n>(func, m, n);
      }

      template<typename CoDiType, typename Func,
               template<typename> class InputVectorType = StdVectorAdapter,
               template<typename> class OutputVectorType = StdVectorAdapter>
      static CODI_INLINE HandleSwitch<Func, CoDiType, InputVectorType, OutputVectorType> create(const Func& func, size_t m, size_t n) {
        return HandleSwitch<Func, CoDiType, InputVectorType, OutputVectorType>(func, m, n);
      }

      template<typename CoDiType, typename Func,
               size_t m,
               size_t n>
      static CODI_INLINE HandleSwitch<Func, CoDiType, StdArrayAdapter<n>::template Type, StdArrayAdapter<m>::template Type> createFixed(const Func& func) {
        return HandleSwitch<Func, CoDiType, StdArrayAdapter<n>::template Type, StdArrayAdapter<m>::template Type>(func, m, n);
      }

      template<typename T = double>
      static CODI_INLINE Jacobian<StdVectorAdapter<T>> createJacobian(size_t m, size_t n) {
        return Jacobian<StdVectorAdapter<T>>(m, n);
      }

      template<size_t m, size_t n, typename T = double>
      static CODI_INLINE Jacobian<typename StdArrayAdapter<m*n>::template Type<T>> createJacobianFixed() {
        return Jacobian<typename StdArrayAdapter<m*n>::template Type<T>>(m, n);
      }

      template<typename T = double>
      static CODI_INLINE Hessian<StdVectorAdapter<T>> createHessian(size_t m, size_t n) {
        return Hessian<StdVectorAdapter<T>>(m, n);
      }

      template<size_t m, size_t n, typename T = double>
      static CODI_INLINE Hessian<typename StdArrayAdapter<m*n*n>::template Type<T>> createHessianFixed() {
        return Hessian<typename StdArrayAdapter<m*n*n>::template Type<T>>(m, n);
      }

      template<typename Func, typename VecX, typename VecY>
      static CODI_INLINE void evalPrimal(
          Func& func,
          const VecX& x,
          VecY& y) {
        auto h = createDefault(func, y.size(), x.size());
        evalHandlePrimal(h, x, y);
      }

      template<typename Func, typename VecX, typename Jac>
      static CODI_INLINE void evalJacobian(
          Func& func,
          const VecX& x,
          const size_t ySize,
          Jac& jac) {
        auto h = createDefault(func, ySize, x.size());
        evalHandleJacobian(h, x, jac);
      }

      template<typename Func, typename VecX, typename Hes>
      static CODI_INLINE void evalHessian(
          Func& func,
          const VecX& x,
          const size_t ySize,
          Hes& hes) {
        auto h = createDefault2nd(func, ySize, x.size());
        evalHandleHessian(h, x, hes);
      }

      template<typename Func, typename VecX, typename VecY, typename Jac>
      static CODI_INLINE void evalPrimalAndJacobian(
          Func& func,
          const VecX& x,
          VecY& y,
          Jac& jac) {
        auto h = createDefault(func, y.size(), x.size());
        evalHandlePrimalAndJacobian(h, x, y, jac);
      }

      template<typename Func, typename VecX, typename VecY, typename Hes>
      static CODI_INLINE void evalPrimalAndHessian(
          Func& func,
          const VecX& x,
          VecY& y,
          Hes& hes) {
        auto h = createDefault2nd(func, y.size(), x.size());
        evalHandlePrimalAndHessian(h, x, y, hes);
      }

      template<typename Func, typename VecX, typename VecY, typename Jac, typename Hes>
      static CODI_INLINE void evalPrimalAndJacobianAndHessian(
          Func& func,
          const VecX& x,
          VecY& y,
          Jac& jac,
          Hes& hes) {
        auto h = createDefault2nd(func, y.size(), x.size());
        evalPrimalAndJacobianAndHessian(h, x, y, jac, hes);
      }

      template<typename Handle, typename VecX, typename VecY>
      static CODI_INLINE void evalHandlePrimal(
          Handle& handle,
          const VecX& x,
          VecY& y) {
        handle.computePrimal(x, y);
      }

      template<typename Handle, typename VecX, typename Jac>
      static CODI_INLINE void evalHandleJacobian(
          Handle& handle,
          const VecX& x,
          Jac& jac) {
        DummyVector dv;
        handle.computeJacobian(x, jac, dv);
      }

      template<typename Handle, typename VecX, typename Hes>
      static CODI_INLINE void evalHandleHessian(
          Handle& handle,
          const VecX& x,
          Hes& hes) {
        DummyVector dv;
        DummyJacobian dj;
        handle.computeHessian(x, hes, dv, dj);
      }

      template<typename Handle, typename VecX, typename VecY, typename Jac>
      static CODI_INLINE void evalHandlePrimalAndJacobian(
          Handle& handle,
          const VecX& x,
          VecY& y,
          Jac& jac) {
        handle.computeJacobian(x, jac, y);
      }

      template<typename Handle, typename VecX, typename VecY, typename Hes>
      static CODI_INLINE void evalHandlePrimalAndHessian(
          Handle& handle,
          const VecX& x,
          VecY& y,
          Hes& hes) {
        DummyJacobian dj;
        handle.computeHessian(x, hes, y, dj);
      }

      template<typename Handle, typename VecX, typename VecY, typename Jac, typename Hes>
      static CODI_INLINE void evalHandlePrimalAndJacobianAndHessian(
          Handle& handle,
          const VecX& x,
          VecY& y,
          Jac& jac,
          Hes& hes) {
        handle.computeHessian(x, hes, y, jac);
      }
  };

}
