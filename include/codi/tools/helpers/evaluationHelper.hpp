#pragma once

#include <array>
#include <vector>

#include "../../../codi.hpp"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../traits/gradientTraits.hpp"
#include "../../traits/tapeTraits.hpp"
#include "../data/dummy.hpp"
#include "../data/jacobian.hpp"
#include "../data/hessian.hpp"
#include "tapeHelper.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template <typename _Func,
            typename _Type,
            typename _InputStore = std::vector<_Type>,
            typename _OutputStore = std::vector<_Type>
            >
  struct EvaluationHandleBase {
    public:

      using Func = CODI_DECLARE_DEFAULT(_Func, CODI_TEMPLATE(void ()(_InputStore const&, _OutputStore&)));
      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
      using InputStore = CODI_DECLARE_DEFAULT(_InputStore, std::vector<Type>);
      using OutputStore = CODI_DECLARE_DEFAULT(_OutputStore, std::vector<Type>);

    protected:

      size_t m;
      size_t n;

      Func& func;

      InputStore x;
      OutputStore y;

      DummyVector dummyVector;
      DummyJacobian dummyJacobian;

    public:

      EvaluationHandleBase(Func& func, size_t m, size_t n) : m(m), n(n), func(func), x(n), y(m),
          dummyVector(), dummyJacobian() {}

      template<typename VecX>
      void setAllPrimals(VecX const& locX);

      template<typename VecY>
      void getAllPrimals(VecY& locY);

      template<typename VecX, typename VecY>
      void computePrimal(VecX const& locX, VecY& locY);

      template<typename VecX, typename Jac, typename VecY>
      void computeJacobian(VecX const& locX, Jac& jac, VecY& locY);

      template<typename VecX, typename Hes, typename VecY, typename Jac>
      void computeHessian(VecX const& locX, Hes& hes, VecY& locY, Jac& jac);

    protected:

      void eval() {
        func(x, y);
      }
  };

  template <typename _Func,
            typename _Type,
            typename _InputStore = std::vector<_Type>,
            typename _OutputStore = std::vector<_Type>
            >
  struct ForwardHandle : public EvaluationHandleBase<_Func, _Type, _InputStore, _OutputStore> {
    public:

      using Func = CODI_DECLARE_DEFAULT(_Func, CODI_TEMPLATE(void ()(_InputStore const&, _OutputStore&)));
      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
      using InputStore = CODI_DECLARE_DEFAULT(_InputStore, std::vector<Type>);
      using OutputStore = CODI_DECLARE_DEFAULT(_OutputStore, std::vector<Type>);

      using Base = EvaluationHandleBase<Func, Type, InputStore, OutputStore>;

      // Forward constructors of the base.
      using Base::EvaluationHandleBase;

      template<typename VecX>
      void setAllPrimals(VecX const& locX) {
        codiAssert(locX.size() <= this->x.size());
        for (size_t j = 0; j < locX.size(); j += 1) {
          this->x[j] = locX[j];
        }
      }

      template<typename VecY>
      void getAllPrimals(VecY& locY) {
        codiAssert(locY.size() <= this->y.size());
        for (size_t i = 0; i < locY.size(); i += 1) {
          locY[i] = this->y[i].getValue();
        }
      }

      template<typename VecX, typename VecY>
      void computePrimal(VecX const& locX, VecY& locY) {
        setAllPrimals(locX);

        this->eval();

        getAllPrimals(locY);
      }

      template<typename VecX, typename Jac, typename VecY>
      void computeJacobian(VecX const& locX, Jac& jac, VecY& locY) {
        setAllPrimals(locX);

        JacobianConvertWrapper<Jac> wrapper(jac);

        // First order derivatives should always exist
        using GradientTraits1st = GradientTraits::TraitsImplementation<typename Type::Gradient>;
        size_t constexpr VectorSizeFirstOrder = GradientTraits1st::dim;

        for (size_t j = 0; j < locX.size(); j+= VectorSizeFirstOrder) {
          for (size_t vecPos = 0; vecPos < VectorSizeFirstOrder && j + vecPos < locX.size(); vecPos += 1) {
            GradientTraits1st::at(this->x[j + vecPos].gradient(), vecPos) = 1.0;
          }

          this->eval();

          if (0 == j) {
            getAllPrimals(locY);
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

      template<typename VecX, typename Hes, typename VecY, typename Jac>
      void computeHessian(VecX const& locX, Hes& hes, VecY& locY, Jac& jac) {
        setAllPrimals(locX);

        // First order derivatives should always exist
        using GradientTraits1st = GradientTraits::TraitsImplementation<typename Type::Gradient>;
        size_t constexpr VectorSizeFirstOrder = GradientTraits1st::dim;

        // Define these here since not all types have second order derivatives
        using GradientTraits2nd = GradientTraits::TraitsImplementation<typename Type::Real::Gradient>;
        size_t constexpr VectorSizeSecondOrder = GradientTraits2nd::dim;


        for (size_t k = 0; k < locX.size(); k+= VectorSizeFirstOrder) {

          // Set derivatives from k to k + vecSize_k
          for (size_t vecPos = 0; vecPos < VectorSizeFirstOrder && k + vecPos < locX.size(); vecPos += 1) {
            GradientTraits1st::at(this->x[k + vecPos].gradient(), vecPos).value() = 1.0;
          }

          // The j = k init is no problem, it will evaluated slightly more elements around the diagonal
          for (size_t j = k; j < locX.size(); j+= VectorSizeSecondOrder) {

            // Set derivatives from j to j + vecSize_j
            for (size_t vecPos = 0; vecPos < VectorSizeSecondOrder && j + vecPos < locX.size(); vecPos += 1) {
              GradientTraits2nd::at(this->x[j + vecPos].value().gradient(), vecPos) = 1.0;
            }

            this->eval();

            if (0 == j && 0 == k) {
              getAllPrimals(locY);
            }

            // Extract all hessian values, this populates the hessian from (j,k) to (j + vecSize_j, k + vecSize_k).
            for (size_t i = 0; i < this->y.size(); i += 1) {
              for (size_t vecPos1st = 0; vecPos1st < VectorSizeFirstOrder && k + vecPos1st < locX.size(); vecPos1st += 1) {
                for (size_t vecPos2nd = 0; vecPos2nd < VectorSizeSecondOrder && j + vecPos2nd < locX.size(); vecPos2nd += 1) {
                  auto& firstGrad = GradientTraits1st::at(this->y[i].gradient(), vecPos1st);
                  auto& secondGrad = GradientTraits2nd::at(firstGrad.gradient(), vecPos2nd);

                  hes(i, j + vecPos2nd, k + vecPos1st) = secondGrad;
                  hes(i, k + vecPos1st, j + vecPos2nd) = secondGrad; // symmetry
                }
              }

              if (k == 0) {
                for (size_t vecPos = 0; vecPos < VectorSizeSecondOrder && j + vecPos < locX.size(); vecPos += 1) {
                  jac(i, j + vecPos) = GradientTraits2nd::at(this->y[i].value().gradient(), vecPos);
                }
              }
            }

            // Reset the derivative seeding
            for (size_t vecPos = 0; vecPos < VectorSizeSecondOrder && j + vecPos < locX.size(); vecPos += 1) {
              GradientTraits2nd::at(this->x[j + vecPos].value().gradient(), vecPos) = 0.0;
            }
          }

          // Reset the derivative seeding
          for (size_t vecPos = 0; vecPos < VectorSizeFirstOrder && k + vecPos < locX.size(); vecPos += 1) {
            GradientTraits1st::at(this->x[k + vecPos].gradient(), vecPos).value() = 0.0;
          }
        }
      }
  };

  template <typename _Func,
            typename _Type,
            typename _InputStore = std::vector<_Type>,
            typename _OutputStore = std::vector<_Type>
            >
  struct ReverseHandleBase : public EvaluationHandleBase<_Func, _Type, _InputStore, _OutputStore> {
    public:

      using Func = CODI_DECLARE_DEFAULT(_Func, CODI_TEMPLATE(void ()(_InputStore const&, _OutputStore&)));
      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
      using InputStore = CODI_DECLARE_DEFAULT(_InputStore, std::vector<Type>);
      using OutputStore = CODI_DECLARE_DEFAULT(_OutputStore, std::vector<Type>);

      using Base = EvaluationHandleBase<Func, Type, InputStore, OutputStore>;

    protected:

      TapeHelper<Type> th;

    public:

      ReverseHandleBase(Func& func, size_t m, size_t n) :
        Base(func, m, n),
        th()
      {}

      template<typename VecX>
      void setAllPrimals(VecX const& locX, bool reg) {
        codiAssert(locX.size() <= this->x.size());
        for (size_t j = 0; j < locX.size(); j += 1) {
          this->x[j] = locX[j];

          if (reg) {
            th.registerInput(this->x[j]);
          }
        }
      }

      template<typename VecY>
      void getAllPrimals(VecY& locY, bool reg) {
        codiAssert(locY.size() <= this->y.size());
        for (size_t i = 0; i < this->y.size(); i += 1) {
          if (reg) {
            th.registerOutput(this->y[i]);
          }

          locY[i] = this->y[i].getValue();
        }
      }

      template<typename VecX, typename VecY>
      void computePrimal(VecX const& locX, VecY& locY) {
        setAllPrimals(locX, false);

        this->eval();

        getAllPrimals(locY, false);
      }

      template<typename VecX, typename Jac, typename VecY>
      void computeJacobian(VecX const& locX, Jac& jac, VecY& locY) {
        recordTape(locX, locY);

        th.evalJacobian(jac);
      }

      template<typename VecX, typename Hes, typename VecY, typename Jac>
      void computeHessian(VecX const& locX, Hes& hes, VecY& locY, Jac& jac);

    protected:

      template<typename VecX, typename VecY>
      void recordTape(VecX const& locX, VecY& locY) {
        th.startRecording();
        setAllPrimals(locX, true);

        this->eval();

        getAllPrimals(locY, true);
        th.stopRecording();
      }
  };

  template <typename _Func,
            typename _Type,
            typename _InputStore = std::vector<_Type>,
            typename _OutputStore = std::vector<_Type>
            >
  struct ReverseHandlePrimalValueTapes : public ReverseHandleBase<_Func, _Type, _InputStore, _OutputStore> {
    public:
      using Func = CODI_DECLARE_DEFAULT(_Func, CODI_TEMPLATE(void ()(_InputStore const&, _OutputStore&)));
      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
      using InputStore = CODI_DECLARE_DEFAULT(_InputStore, std::vector<Type>);
      using OutputStore = CODI_DECLARE_DEFAULT(_OutputStore, std::vector<Type>);

      using Base = ReverseHandleBase<Func, Type, InputStore, OutputStore>;

      // Use constructors of the base class.
      using Base::ReverseHandleBase;

      template<typename VecX, typename Hes, typename VecY, typename Jac>
      void computeHessian(VecX const& locX, Hes& hes, VecY& locY, Jac& jac) {
        this->recordTape(locX, locY);

        this->th.evalHessian(hes, jac);
      }
  };

  template <typename _Func,
            typename _Type,
            typename _InputStore = std::vector<_Type>,
            typename _OutputStore = std::vector<_Type>
            >
  struct ReverseHandleJacobiTapes : public ReverseHandleBase<_Func, _Type, _InputStore, _OutputStore> {
    public:
      using Func = CODI_DECLARE_DEFAULT(_Func, CODI_TEMPLATE(void ()(_InputStore const&, _OutputStore&)));
      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
      using InputStore = CODI_DECLARE_DEFAULT(_InputStore, std::vector<Type>);
      using OutputStore = CODI_DECLARE_DEFAULT(_OutputStore, std::vector<Type>);

      using Base = ReverseHandleBase<Func, Type, InputStore, OutputStore>;

      // Use constructors of the base class.
      using Base::ReverseHandleBase;

      template<typename VecX, typename Hes, typename VecY, typename Jac>
      void computeHessian(VecX const& locX, Hes& hes, VecY& locY, Jac& jac) {
        this->setAllPrimals(locX, false);

        Algorithms<Type>::computeHessian(this->func, this->x, this->y, hes, jac);

        this->getAllPrimals(locY, false);
      }
  };

  template <typename _Func,
            typename _Type,
            typename _InputStore = std::vector<_Type>,
            typename _OutputStore = std::vector<_Type>,
            typename = void
            >
  struct EvaluationHandle : public EvaluationHandleBase<_Func, _Type, _InputStore, _OutputStore> {};

  template <typename _Func,
            typename _Type,
            typename _InputStore,
            typename _OutputStore
            >
  struct EvaluationHandle<_Func, _Type, _InputStore, _OutputStore, TapeTraits::EnableIfForwardTape<typename _Type::Tape>> :
      public ForwardHandle<_Func, _Type, _InputStore, _OutputStore> {

      using ForwardHandle<_Func, _Type, _InputStore, _OutputStore>::ForwardHandle;
  };

  template <typename _Func,
            typename _Type,
            typename _InputStore,
            typename _OutputStore
            >
  struct EvaluationHandle<_Func, _Type, _InputStore, _OutputStore, TapeTraits::EnableIfJacobianTape<typename _Type::Tape>> :
      public ReverseHandleJacobiTapes<_Func, _Type, _InputStore, _OutputStore> {

      using ReverseHandleJacobiTapes<_Func, _Type, _InputStore, _OutputStore>::ReverseHandleJacobiTapes;
  };

  template <typename _Func,
            typename _Type,
            typename _InputStore,
            typename _OutputStore
            >
  struct EvaluationHandle<_Func, _Type, _InputStore, _OutputStore, TapeTraits::EnableIfPrimalValueTape<typename _Type::Tape>> :
      public ReverseHandlePrimalValueTapes<_Func, _Type, _InputStore, _OutputStore> {

      using ReverseHandlePrimalValueTapes<_Func, _Type, _InputStore, _OutputStore>::ReverseHandlePrimalValueTapes;
  };

  struct EvaluationHelper {
    public:

      using JacobianComputationType = RealForwardVec<4>;
      using HessianComputationType = RealForwardGen<RealForwardVec<4>, Direction<RealForwardVec<4>, 4>>;

      template <typename Func>
      using DefaultHandle = ForwardHandle<Func, JacobianComputationType>;

      template <typename Func>
      using DefaultHandle2nd = ForwardHandle<Func, HessianComputationType>;

      template <typename Func, size_t m, size_t n>
      using DefaultHandleFixed = ForwardHandle<Func, JacobianComputationType, std::array<JacobianComputationType, n>, std::array<JacobianComputationType, m>>;

      template <typename Func, size_t m, size_t n>
      using DefaultHandleFixed2nd = ForwardHandle<Func, HessianComputationType, std::array<HessianComputationType, n>, std::array<HessianComputationType, m>>;

      template<typename Func>
      static CODI_INLINE DefaultHandle<Func> createHandleDefault(Func& func, size_t m, size_t n) {
        return DefaultHandle<Func>(func, m, n);
      }

      template<size_t m, size_t n, typename Func>
      static CODI_INLINE DefaultHandleFixed<Func, m, n> createHandleDefaultFixed(Func& func) {
        return DefaultHandleFixed<Func, m, n>(func, m, n);
      }

      template<typename Func>
      static CODI_INLINE DefaultHandle2nd<Func> createHandleDefault2nd(Func& func, size_t m, size_t n) {
        return DefaultHandle2nd<Func>(func, m, n);
      }

      template<size_t m, size_t n, typename Func>
      static CODI_INLINE DefaultHandleFixed2nd<Func, m, n> createHandleDefaultFixed2nd(Func& func) {
        return DefaultHandleFixed2nd<Func, m, n>(func, m, n);
      }

      template<typename Type, typename Func>
      static CODI_INLINE EvaluationHandle<Func, Type> createHandle(Func& func, size_t m, size_t n) {
        return EvaluationHandle<Func, Type>(func, m, n);
      }

      template<typename Type, size_t m, size_t n, typename Func>
      static CODI_INLINE EvaluationHandle<Func, Type, std::array<Type, n>, std::array<Type, m>> createHandleFixed(Func& func) {
        return EvaluationHandle<Func, Type, std::array<Type, n>, std::array<Type, m>>(func, m, n);
      }

      template<typename Type,
               typename InputStore,
               typename OutputStore,
               typename Func
               >
      static CODI_INLINE EvaluationHandle<Func, Type, InputStore, OutputStore> createHandleFull(Func& func, size_t m, size_t n) {
        return EvaluationHandle<Func, Type, InputStore, OutputStore>(func, m, n);
      }

      template<typename T = double>
      static CODI_INLINE Jacobian<T> createJacobian(size_t m, size_t n) {
        return Jacobian<T>(m, n);
      }

      template<size_t m, size_t n, typename T = double>
      static CODI_INLINE Jacobian<T, std::array<T, m * n>> createJacobianFixed() {
        return Jacobian<T, std::array<T, m * n>>(m, n);
      }

      template<typename T = double>
      static CODI_INLINE Hessian<T> createHessian(size_t m, size_t n) {
        return Hessian<T>(m, n);
      }

      template<size_t m, size_t n, typename T = double>
      static CODI_INLINE Hessian<T, std::array<T, m * n * n>> createHessianFixed() {
        return Hessian<T, std::array<T, m * n * n>>(m, n);
      }

      template<typename Func, typename VecX, typename VecY>
      static CODI_INLINE void evalPrimal(
          Func& func,
          VecX const& x,
          VecY& y) {
        auto h = createHandleDefault(func, y.size(), x.size());
        evalHandlePrimal(h, x, y);
      }

      template<typename Func, typename VecX, typename Jac>
      static CODI_INLINE void evalJacobian(
          Func& func,
          VecX const& x,
          size_t const ySize,
          Jac& jac) {
        auto h = createHandleDefault(func, ySize, x.size());
        evalHandleJacobian(h, x, jac);
      }

      template<typename Func, typename VecX, typename Hes>
      static CODI_INLINE void evalHessian(
          Func& func,
          VecX const& x,
          size_t const ySize,
          Hes& hes) {
        auto h = createHandleDefault2nd(func, ySize, x.size());
        evalHandleHessian(h, x, hes);
      }

      template<typename Func, typename VecX, typename VecY, typename Jac>
      static CODI_INLINE void evalPrimalAndJacobian(
          Func& func,
          VecX const& x,
          VecY& y,
          Jac& jac) {
        auto h = createHandleDefault(func, y.size(), x.size());
        evalHandlePrimalAndJacobian(h, x, y, jac);
      }

      template<typename Func, typename VecX, typename VecY, typename Hes>
      static CODI_INLINE void evalPrimalAndHessian(
          Func& func,
          VecX const& x,
          VecY& y,
          Hes& hes) {
        auto h = createHandleDefault2nd(func, y.size(), x.size());
        evalHandlePrimalAndHessian(h, x, y, hes);
      }

      template<typename Func, typename VecX, typename VecY, typename Jac, typename Hes>
      static CODI_INLINE void evalPrimalAndJacobianAndHessian(
          Func& func,
          VecX const& x,
          VecY& y,
          Jac& jac,
          Hes& hes) {
        auto h = createHandleDefault2nd(func, y.size(), x.size());
        evalHandlePrimalAndJacobianAndHessian(h, x, y, jac, hes);
      }

      template<typename Func, typename VecX, typename Jac, typename Hes>
      static CODI_INLINE void evalJacobianAndHessian(
          Func& func,
          VecX const& x,
          size_t ySize,
          Jac& jac,
          Hes& hes) {
        auto h = createHandleDefault2nd(func, ySize, x.size());
        evalHandleJacobianAndHessian(h, x, jac, hes);
      }

      template<typename Handle, typename VecX, typename VecY>
      static CODI_INLINE void evalHandlePrimal(
          Handle& handle,
          VecX const& x,
          VecY& y) {
        handle.computePrimal(x, y);
      }

      template<typename Handle, typename VecX, typename Jac>
      static CODI_INLINE void evalHandleJacobian(
          Handle& handle,
          VecX const& x,
          Jac& jac) {
        DummyVector dv;
        handle.computeJacobian(x, jac, dv);
      }

      template<typename Handle, typename VecX, typename Hes>
      static CODI_INLINE void evalHandleHessian(
          Handle& handle,
          VecX const& x,
          Hes& hes) {
        DummyVector dv;
        DummyJacobian dj;
        handle.computeHessian(x, hes, dv, dj);
      }

      template<typename Handle, typename VecX, typename VecY, typename Jac>
      static CODI_INLINE void evalHandlePrimalAndJacobian(
          Handle& handle,
          VecX const& x,
          VecY& y,
          Jac& jac) {
        handle.computeJacobian(x, jac, y);
      }

      template<typename Handle, typename VecX, typename VecY, typename Hes>
      static CODI_INLINE void evalHandlePrimalAndHessian(
          Handle& handle,
          VecX const& x,
          VecY& y,
          Hes& hes) {
        DummyJacobian dj;
        handle.computeHessian(x, hes, y, dj);
      }

      template<typename Handle, typename VecX, typename VecY, typename Jac, typename Hes>
      static CODI_INLINE void evalHandlePrimalAndJacobianAndHessian(
          Handle& handle,
          VecX const& x,
          VecY& y,
          Jac& jac,
          Hes& hes) {
        handle.computeHessian(x, hes, y, jac);
      }

      template<typename Handle, typename VecX, typename Jac, typename Hes>
      static CODI_INLINE void evalHandleJacobianAndHessian(
          Handle& handle,
          VecX const& x,
          Jac& jac,
          Hes& hes) {
        DummyVector dv;
        handle.computeHessian(x, hes, dv, jac);
      }
  };

}
