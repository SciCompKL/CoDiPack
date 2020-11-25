#pragma once

#include <utility>

#include "../aux/binomial.hpp"
#include "../aux/compileTimeLoop.hpp"
#include "../aux/exceptions.hpp"
#include "../aux/macros.hpp"
#include "../config.h"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../traits/realTraits.hpp"


/** \copydoc codi::Namespace */
namespace codi {

  namespace HigherOrderAccessImpl {
    CODI_INLINE size_t constexpr maximumDerivatives(size_t selectionDepth, size_t order) {
      return binomial(selectionDepth, order);
    }

    CODI_INLINE size_t constexpr maximumDerivativesPrimalBranch(size_t selectionDepth, size_t order) {
      return binomial(selectionDepth - 1, order);
    }

    CODI_INLINE size_t constexpr isPrimalBranch(size_t selectionDepth, size_t order, size_t l) {
      return l < maximumDerivativesPrimalBranch(selectionDepth, order);
    }

    template<typename Type, size_t selectionDepth, size_t order, size_t l>
    struct CheckCompileTimeValues {

        static_assert (selectionDepth <= RealTraits::MaxDerivativeOrder<Type>(),
                       "Selection depth can not be higher than the maximum derivative order." );
        static_assert (order <= selectionDepth,
                       "Derivative order can not be higher than the selection depth.");
        static_assert (l < maximumDerivatives(selectionDepth, order),
                       "Selected derivative can not be greater than the number of available derivatives for that"
                       "order." );

        static bool constexpr isValid = true;

    };

    template<typename Type, bool constant, size_t selectionDepth, size_t order, size_t l, bool primalBranch = isPrimalBranch(selectionDepth, order, l)>
    struct SelectCompileTime;

    template<typename _Type, bool constant, size_t selectionDepth, size_t order, size_t l>
    struct SelectCompileTime<_Type, constant, selectionDepth, order, l, true> {
      public:
        using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

        using Inner = SelectCompileTime<typename Type::Real, constant, selectionDepth - 1, order, l>;
        using ArgType = typename std::conditional<constant, Type const, Type>::type;
        using RType = typename Inner::RType;

        static_assert (CheckCompileTimeValues<Type, selectionDepth, order, l>::isValid, "Checks inside of type.");

        static RType& select(ArgType& value) {
           return Inner::select(value.value());
        }
    };

    template<typename _Type, bool constant, size_t selectionDepth, size_t order, size_t l>
    struct SelectCompileTime<_Type, constant, selectionDepth, order, l, false> {
      public:
        using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

        using Inner = SelectCompileTime<typename Type::Real, constant, selectionDepth - 1, order - 1, l - maximumDerivativesPrimalBranch(selectionDepth, order)>;
        using ArgType = typename std::conditional<constant, Type const, Type>::type;
        using RType = typename Inner::RType;

        static_assert (CheckCompileTimeValues<Type, selectionDepth, order, l>::isValid, "Checks inside of type.");

        static RType& select(ArgType& value) {
           return Inner::select(value.gradient());
        }
    };

    template<typename Type, bool constant>
    struct SelectCompileTime<Type, constant, 0, 0, 0, true> {
      public:
        using ArgType = typename std::conditional<constant, Type const, Type>::type;
        using RType = ArgType;

        static RType& select(ArgType& value) {
           return value;
        }
    };

    template<typename _Type, bool constant, size_t _selectionDepth>
    struct SelectRunTime {

        using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
        static size_t constexpr selectionDepth = CODI_DECLARE_DEFAULT(_selectionDepth, /*TODO*/ 0);

        static_assert (std::is_same<typename Type::Real, typename Type::Gradient>::value, "CoDiPack type needs to have the same real and gradient value for run time derivative selection.");
        static_assert (selectionDepth <= RealTraits::MaxDerivativeOrder<Type>(), "Selection depth can not be higher than the maximum derivative order" );

        using Inner = SelectRunTime<typename Type::Real, constant, selectionDepth - 1>;
        using ArgType = typename std::conditional<constant, Type const, Type>::type;
        using RType = typename Inner::RType;

        static RType& select(ArgType& v, size_t order, size_t l) {
          size_t const maxDerivativesPrimalBranch = binomial(selectionDepth - 1, order);
          if(l < maxDerivativesPrimalBranch) {
            return Inner::select(v.value(), order, l);
          } else {
            return Inner::select(v.gradient(), order - 1, l - maxDerivativesPrimalBranch);
          }
        }
    };

    template<typename _Type, bool constant>
    struct SelectRunTime<_Type, constant, 0> {

        using Type = CODI_DECLARE_DEFAULT(_Type, double);
        using ArgType = typename std::conditional<constant, Type const, Type>::type;
        using RType = ArgType;

        static RType& select(ArgType& v, size_t order, size_t l) {
          return v;
        }
    };
  }


  template<typename _Type>
  struct HigherOrderAccess {
    public:

      using Type = CODI_DECLARE_DEFAULT(_Type, CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

      template<bool constant, size_t selectionDepth>
      using SelectRunTime = HigherOrderAccessImpl::SelectRunTime<Type, constant, selectionDepth>;

      template<bool constant, size_t selectionDepth, size_t order, size_t l>
      using SelectCompileTime = HigherOrderAccessImpl::SelectCompileTime<Type, constant, selectionDepth, order, l>;


      template<size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static typename SelectRunTime<true, selectionDepth>::RType const& derivative(Type const& v, size_t order, size_t l) {

        checkRuntimeSelection<selectionDepth>(order, l);

        return SelectRunTime<true, selectionDepth>::select(v, order, l);
      }

      template<size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static typename SelectRunTime<false, selectionDepth>::RType& derivative(Type& v, size_t order, size_t l) {

        checkRuntimeSelection<selectionDepth>(order, l);

        return SelectRunTime<false, selectionDepth>::select(v, order, l);
      }

      template<typename Derivative, size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static void setAllDerivatives(Type& v, size_t order, Derivative const& d) {
        size_t const maxDerivatives = binomial(selectionDepth, order);
        for(size_t i = 0; i < maxDerivatives; i += 1) {
          derivative<selectionDepth>(v, order, i) = d;
        }
      }

      template<typename Derivative, size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static void setAllDerivativesForward(Type& v, size_t order, Derivative const& d) {
        HigherOrderAccess<typename Type::Real>::template setAllDerivatives<Derivative, selectionDepth - 1>(v.value(), order, d);
      }

      template<typename Derivative, size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static void setAllDerivativesReverse(Type& v, size_t order, Derivative const& d) {
        HigherOrderAccess<typename Type::Gradient>::template setAllDerivatives<Derivative, selectionDepth - 1>(v.gradient(), order - 1, d);
      }

      template<size_t order, size_t l, size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static typename SelectCompileTime<true, selectionDepth, order, l>::RType const& derivative(Type const& v) {
        return SelectCompileTime<true, selectionDepth, order, l>::select(v);
      }

      template<size_t order, size_t l, size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static typename SelectCompileTime<false, selectionDepth, order, l>::RType& derivative(Type& v) {
        return SelectCompileTime<false, selectionDepth, order, l>::select(v);
      }

      template<size_t order, typename Derivative, size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static void setAllDerivatives(Type& v, Derivative const& d) {
        CompileTimeLoop<HigherOrderAccessImpl::maximumDerivatives(selectionDepth, order)>::eval(CallSetDerivative<order, Derivative, selectionDepth>{}, v, d);
      }

      template<size_t order, typename Derivative, size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static void setAllDerivativesForward(Type& v, Derivative const& d) {
        HigherOrderAccess<typename Type::Real>::template setAllDerivatives<order, Derivative, selectionDepth - 1>(v.value(), d);
      }

      template<size_t order, typename Derivative, size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static void setAllDerivativesReverse(Type& v, Derivative const& d) {
        HigherOrderAccess<typename Type::Gradient>::template setAllDerivatives<order - 1, Derivative, selectionDepth - 1>(v.gradient(), d);
      }

    private:

      template<size_t selectionDepth>
      static void checkRuntimeSelection(size_t order, size_t l) {
        if(order > selectionDepth) {
          CODI_EXCEPTION("The derivative order must be smaller or equal than the maximum possible derivative. order: %d, max derivative: %d.", order, selectionDepth);
        }

        size_t numberDerivatives = binomial(selectionDepth, order);
        if(l >= numberDerivatives) {
          CODI_EXCEPTION("The selected derivative must be smaller than the maximum number of derivatives. selected: %d, number derivatives: %d.", l, numberDerivatives);
        }
      }

      template<size_t order, typename Derivative, size_t selectionDepth>
      struct CallSetDerivative {
          template<size_t pos>
          void operator()(std::integral_constant<size_t, pos>, Type& v, Derivative const& d) {
            HigherOrderAccess::derivative<order, pos - 1, selectionDepth>(v) = d;
          }
      };
  };
}
