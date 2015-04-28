#pragma once

#include <cmath>

#include "configure.h"
#include "typeTraits.hpp"

namespace codi {

  /**
   * The Expression type from which all other types of expression
   * derive. Each member function simply calls the specialized version
   * of the function according to the expression's true type, which is
   * given by its template argument.
   *
   * @template Real  The data type of the primal values and the gradient values.
   */

  template<typename Real, class A>
  struct Expression {

    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    /**
     * Cast the expression to its true type, given by the template
     * argument
     */
    inline const A& cast() const {
      return static_cast<const A&>(*this);
    }

    /**
     * Calculate the gradient of the mathematical operation that this
     * expression represents and pass the result to its argument.
     * For functions f(a), pass df/da to the argument in the
     * first case and pass multiplier*df/da in the second case.
     */
    inline void calcGradient(Real& gradient) const {
      cast().calcGradient(gradient);
    }

    /**
     * As the previous but multiplying the gradient by "multiplier"
     */
    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      cast().calcGradient(gradient, multiplier);
    }

    /**
     * Return the numerical value of the expression
     */
    inline const Real& getValue() const {
      return cast().getValue();
    }

  private:
    /**
     * Intentionally inaccessible to prevent an expression appearing
     * on the left-hand-side of a statement
     */
    Expression& operator=(const Expression&) = delete;
  };

  /**
   * Now define particular types of expression, using static
   * polymorphism via the Curiously Recurring Template Pattern
   */

  /**
   * Add: an expression plus another expression
   */
  template<typename Real, class A, class B>
  struct Add : public Expression<Real, Add<Real, A, B> > {
    Add(const Expression<Real, A>& a, const Expression<Real, B>& b)
      : a_(a.cast()), b_(b.cast()) { }

    /**
     * If f(a,b) = a + b, df/da = 1 and likewise for df/db so simply
     * call a and b's versions of calcGradient
     */
    inline void calcGradient(Real& gradient) const {
      a_.calcGradient(gradient);
      b_.calcGradient(gradient);
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      a_.calcGradient(gradient, multiplier);
      b_.calcGradient(gradient, multiplier);
    }

    inline Real getValue() const {
      return a_.getValue() + b_.getValue();
    }

  private:
    /**
     * Store constant references to the arguments of the addition
     */
    const A& a_;
    const B& b_;
  };

  /**
   * Overload the addition operator for Expression arguments to return
   * an Add type
   */
  template<typename Real, class A, class B>
  inline
  Add<Real, A, B> operator+(const Expression<Real, A>& a,
                            const Expression<Real, B>& b) {
    return Add<Real, A, B>(a.cast(), b.cast());
  }

  /**
   * Subtract: an expression minus another expression
   */
  template<typename Real, class A, class B>
  struct Subtract : public Expression<Real, Subtract<Real, A, B> > {
    Subtract(const Expression<Real, A>& a, const Expression<Real, B>& b)
      : a_(a.cast()), b_(b.cast()) { }

    /**
     * If f(a,b) = a - b, df/da = 1 and df/db = -1
     */
    inline void calcGradient(Real& gradient) const {
      a_.calcGradient(gradient);
      b_.calcGradient(gradient, -1.0);
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      a_.calcGradient(gradient, multiplier);
      b_.calcGradient(gradient, -multiplier);
    }

    inline Real getValue() const {
      return a_.getValue() - b_.getValue();
    }

  private:
    const A& a_;
    const B& b_;
  };

  /*
   *  Overload subtraction operator for Expression arguments
   */
  template<typename Real, class A, class B>
  inline
  Subtract<Real, A, B> operator-(const Expression<Real, A>& a,
                           const Expression<Real, B>& b) {
    return Subtract<Real, A, B>(a.cast(), b.cast());
  }

  /**
   * Multiply: an expression multiplied by another expression
   */

  template <typename Real, class A, class B>
  struct Multiply : public Expression<Real, Multiply<Real, A,B> > {
    Multiply(const Expression<Real, A>& a, const Expression<Real, B>& b)
      : a_(a.cast()), b_(b.cast()), result_(a_.getValue()*b_.getValue()) { }
    // If f(a,b) = a*b then df/da = b and df/db = a
    inline void calcGradient(Real& gradient) const {
      a_.calcGradient(gradient, b_.getValue());
      b_.calcGradient(gradient, a_.getValue());
    }
    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      a_.calcGradient(gradient, b_.getValue()*multiplier);
      b_.calcGradient(gradient, a_.getValue()*multiplier);
    }
    inline const Real& getValue() const {
      return result_;
    }
  private:
    const A& a_;
    const B& b_;
    Real result_;
  };


  /**
   * Overload multiplication operator for Expression arguments
   */
  template<typename Real, class A, class B>
  inline
  Multiply<Real, A, B> operator*(const Expression<Real, A>& a,
                           const Expression<Real, B>& b) {
    return Multiply<Real, A, B>(a.cast(), b.cast());
  }

  /**
   *  Divide: an expression divided by another expression
   */
  template<typename Real, class A, class B>
  struct Divide : public Expression<Real, Divide<Real, A, B> > {
    Divide(const Expression<Real, A>& a, const Expression<Real, B>& b)
      : a_(a.cast()), b_(b.cast()), one_over_b_(1.0 / b_.getValue()),
        result_(a_.getValue() * one_over_b_) { }

    // If f(a,b) = a/b then df/da = 1/b and df/db = -a/(b*b)
    inline void calcGradient(Real& gradient) const {
      a_.calcGradient(gradient, one_over_b_);
      b_.calcGradient(gradient, -result_ * one_over_b_);
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      Real tmp = multiplier * one_over_b_;
      a_.calcGradient(gradient, tmp);
      b_.calcGradient(gradient, -tmp * result_);
    }

    inline const Real& getValue() const {
      return result_;
    }

  private:
    const A& a_;
    const B& b_;
    Real one_over_b_;
    Real result_;
  };

  /**
   * Overload division operator for Expression arguments
   */
  template<typename Real, class A, class B>
  inline
  Divide<Real, A, B> operator/(const Expression<Real, A>& a,
                         const Expression<Real, B>& b) {
    return Divide<Real, A, B>(a.cast(), b.cast());
  }

  /**
   * ScalarAdd: an expression plus a scalar
   */
  template<typename Real, class A>
  struct ScalarAdd : public Expression<Real, ScalarAdd<Real, A> > {
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    ScalarAdd(const Expression<Real, A>& a, const PassiveReal& b)
      : a_(a.cast()), result_(a_.getValue() + b) { }

    inline void calcGradient(Real& gradient) const {
      a_.calcGradient(gradient);
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      a_.calcGradient(gradient, multiplier);
    }

    inline const Real& getValue() const {
      return result_;
    }

  private:
    const A& a_;
    Real result_;
  };

  /**
   *  Overload addition operator for expression plus scalar and scalar
   * plus expression
   */
  template<typename Real, class A>
  inline
  ScalarAdd<Real, A> operator+(const Expression<Real, A>& a, const typename TypeTraits<Real>::PassiveReal& b) {
    return ScalarAdd<Real, A>(a.cast(), b);
  }

  template<typename Real, class A>
  inline
  ScalarAdd<Real, A> operator+(const typename TypeTraits<Real>::PassiveReal& b, const Expression<Real, A>& a) {
    return ScalarAdd<Real, A>(a.cast(), b);
  }

  /**
   * Overload subtraction operator for expression minus scalar to
   * return a ScalarAdd object
   */
  template<typename Real, class A>
  inline
  ScalarAdd<Real, A> operator-(const Expression<Real, A>& a, const typename TypeTraits<Real>::PassiveReal& b) {
    return ScalarAdd<Real, A>(a.cast(), -b);
  }

  /**
   *  ScalarSubtract: scalar minus expression
   */
  template<typename Real, class B>
  struct ScalarSubtract : public Expression<Real, ScalarSubtract<Real, B> > {
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    ScalarSubtract(const PassiveReal& a, const Expression<Real, B>& b)
      : b_(b.cast()), result_(a - b_.getValue()) { }

    inline void calcGradient(Real& gradient) const {
      b_.calcGradient(gradient, -1.0);
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      b_.calcGradient(gradient, -multiplier);
    }

    inline const Real& getValue() const {
      return result_;
    }

  private:
    const B& b_;
    Real result_;
  };

  /**
   *  Overload subtraction operator for scalar minus expression
   */
  template<typename Real, class B>
  inline
  ScalarSubtract<Real, B> operator-(const typename TypeTraits<Real>::PassiveReal& a, const Expression<Real, B>& b) {
    return ScalarSubtract<Real, B>(a, b.cast());
  }

  /**
   *  ScalarMultiply: expression multiplied by scalar
   */
  template<typename Real, class A>
  struct ScalarMultiply : public Expression<Real, ScalarMultiply<Real, A> > {
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    ScalarMultiply(const Expression<Real, A>& a, const PassiveReal& b)
      : a_(a.cast()), b_(b) { }

    inline void calcGradient(Real& gradient) const {
      a_.calcGradient(gradient, b_);
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      a_.calcGradient(gradient, multiplier * b_);
    }

    inline Real getValue() const {
      return a_.getValue() * b_;
    }

  private:
    const A& a_;
    PassiveReal b_;
  };

  /**
   * Overload multiplication operator for expression multiplied by
   * scalar and scalar multiplied by expression
   */
  template<typename Real, class A>
  inline
  ScalarMultiply<Real, A> operator*(const Expression<Real, A>& a, const typename TypeTraits<Real>::PassiveReal& b) {
    return ScalarMultiply<Real, A>(a.cast(), b);
  }

  template<typename Real, class A>
  inline
  ScalarMultiply<Real, A> operator*(const typename TypeTraits<Real>::PassiveReal& b, const Expression<Real, A>& a) {
    return ScalarMultiply<Real, A>(a.cast(), b);
  }

  /**
   *  Overload division operator for expression divided by scalar to
   * return ScalarMultiply object
   */
  template<typename Real, class A>
  inline
  ScalarMultiply<Real, A> operator/(const Expression<Real, A>& a, const typename TypeTraits<Real>::PassiveReal& b) {
    return ScalarMultiply<Real, A>(a.cast(), 1.0 / b);
  }

  /**
   * ScalarDivide: scalar divided by expression
   */
  template<typename Real, class B>
  struct ScalarDivide : public Expression<Real, ScalarDivide<Real, B> > {
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    ScalarDivide(const PassiveReal& a, const Expression<Real, B>& b)
      : b_(b.cast()), one_over_b_(1.0 / b_.getValue()),
        result_(a * one_over_b_) { }

    /**
     *  If f(a,b) = a/b then df/db = -a/(b*b)
     */
    inline void calcGradient(Real& gradient) const {
      b_.calcGradient(gradient, -result_ * one_over_b_);
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      b_.calcGradient(gradient, -multiplier * result_ * one_over_b_);
    }

    inline const Real& getValue() const {
      return result_;
    }

  private:
    const B& b_;
    Real one_over_b_;
    Real result_;
  };

  /**
   *  Overload division operator for scalar divided by expression
   */
  template<typename Real, class B>
  inline
  ScalarDivide<Real, B> operator/(const typename TypeTraits<Real>::PassiveReal& a, const Expression<Real, B>& b) {
    return ScalarDivide<Real, B>(a, b.cast());
  }

  /**
   *  Conditional operators should behave exactly the same as with
   * non-active arguments so in each of the cases below the getValue()
   * function is called to extract the value of the expression
   */
#define CODI_DEFINE_CONDITIONAL(OPERATOR, OP)      \
template<typename Real, class A, class B>          \
inline              \
bool OPERATOR(const Expression<Real, A>& a,        \
      const Expression<Real, B>& b) {      \
  return a.getValue() OP b.getValue();        \
}                \
                              \
template<typename Real, class A>            \
inline              \
bool OPERATOR(const Expression<Real, A>& a, const Real& b) {  \
  return a.getValue() OP b;          \
}                \
                              \
template<typename Real, class B>            \
inline              \
bool OPERATOR(const Real& a, const Expression<Real, B>& b) {  \
  return a OP b.getValue();          \
}

  CODI_DEFINE_CONDITIONAL(operator==, ==)

  CODI_DEFINE_CONDITIONAL(operator!=, !=)

  CODI_DEFINE_CONDITIONAL(operator>, >)

  CODI_DEFINE_CONDITIONAL(operator<, <)

  CODI_DEFINE_CONDITIONAL(operator>=, >=)

  CODI_DEFINE_CONDITIONAL(operator<=, <=)

#undef CODI_DEFINE_CONDITIONAL

  /*
   *  UnaryMinus: negation of expression
   */
  template<typename Real, class A>
  struct UnaryMinus : public Expression<Real, UnaryMinus<Real, A> > {
    UnaryMinus(const Expression<Real, A>& a)
      : a_(a.cast()) { }

    inline void calcGradient(Real& gradient) const {
      a_.calcGradient(gradient, -1.0);
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      a_.calcGradient(gradient, -multiplier);
    }

    inline Real getValue() const {
      return -a_.getValue();
    }

  private:
    const A& a_;
  };

  /**
   *  Overload unary minus of expression
   */
  template<typename Real, class A>
  inline
  UnaryMinus<Real, A> operator-(const Expression<Real, A>& a) {
    return UnaryMinus<Real, A>(a.cast());
  }

  /**
   *  Unary plus: returns the argument
   */
  template<typename Real, class A>
  inline
  A operator+(const Expression<Real, A>& a) {
    return a;
  }
}
/**
 * End namespace codi
 */

/* predefine the struct and the function for higher order derivatrives */
namespace codi {
  template<typename Real, class A> struct Exp;
  template<typename Real, class A> struct Atanh;
}
template <typename Real, class A>
inline  codi:: Exp<Real, A> exp(const codi::Expression<Real, A>& a);
template <typename Real, class A>
inline  codi:: Atanh<Real, A> atanh(const codi::Expression<Real, A>& a);

namespace codi {
  /**
   *  Exponential of an expression
   */
  template<typename Real, class A>
  struct Exp : public Expression<Real, Exp<Real, A> > {
    Exp(const Expression<Real, A>& a)
      : a_(a.cast()), result_(exp(a.getValue())) { }

    /**
     * If f(a) = exp(a) then df/da = exp(a)
     */
    inline void calcGradient(Real& gradient) const {
      a_.calcGradient(gradient, result_);
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      a_.calcGradient(gradient, result_ * multiplier);
    }

    inline const Real& getValue() const {
      return result_;
    }

  private:
    const A& a_;
    Real result_;
  };

  /**
   *  Exponential of an expression
   */
  using std::log;

  template<typename Real, class A>
  struct Atanh : public Expression<Real, Atanh<Real, A> > {
    Atanh(const Expression<Real, A>& a)
      : a_(a.cast()), result_(.5 * log((1 + a.getValue()) / (1 - a.getValue()))) { }

    inline void calcGradient(Real& gradient) const {
      a_.calcGradient(gradient, 1.0 / (1 - a_.getValue() * a_.getValue()));
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      a_.calcGradient(gradient, 1.0 / (1 - a_.getValue() * a_.getValue()) * multiplier);
    }

    inline const Real& getValue() const {
      return result_;
    }

  private:
    const A& a_;
    Real result_;
  };
}

/*
 *  Overload exp for Expression objects
 */
template<typename Real, class A>
inline
codi::Exp<Real, A> exp(const codi::Expression<Real, A>& a) {
  return codi::Exp<Real, A>(a.cast());
}

template<typename Real, class A>
inline
codi::Atanh<Real, A> atanh(const codi::Expression<Real, A>& a) {
  return codi::Atanh<Real, A>(a.cast());
}

/**
 * Enable mathematical functions when the derivative is most easily
 * written in terms of the return value of the function
 */
#define CODI_DEFINE_UNARY_FUNCTION2(OP, FUNC, DERIVATIVE)  \
/* predefine the struct and the function for higher order derivatrives */\
namespace codi {                                                         \
  template<typename Real, class A> struct OP;                            \
}                                                                        \
template <typename Real, class A>                                        \
inline  codi:: OP<Real, A> FUNC(const codi::Expression<Real, A>& a);     \
                                                                         \
namespace codi {            \
  template<typename Real, class A>            \
  struct OP : public Expression<Real, OP<Real, A> > {      \
    OP(const Expression<Real, A>& a)        \
  : a_(a.cast()), result_(FUNC(a.getValue())) { }    \
    inline void calcGradient(Real& gradient) const {      \
  a_.calcGradient(gradient, DERIVATIVE);      \
    }                \
    inline void calcGradient(Real& gradient,        \
           const Real& multiplier) const {  \
  a_.calcGradient(gradient, (DERIVATIVE)*multiplier);  \
    }                \
    inline const Real& getValue() const {      \
  return result_;            \
    }                \
  private:              \
  const A& a_;            \
  Real result_;            \
  };                \
}                \
template <typename Real, class A>            \
inline              \
codi:: OP<Real, A> FUNC(const codi::Expression<Real, A>& a) {    \
  return codi:: OP<Real, A>(a.cast());        \
}

CODI_DEFINE_UNARY_FUNCTION2(Sqrt, sqrt, result_ != 0.0 ? 0.5 / result_ : (Real)0.0)

CODI_DEFINE_UNARY_FUNCTION2(Tanh, tanh, 1 - result_ * result_)
#undef CODI_DEFINE_UNARY_FUNCTION2


/**
 * Enable unary mathematical functions when the derivative is most
 * easily written in terms of the argument of the function; note that
 * the class is in the codi name space but the function is not.
 */
#define CODI_DEFINE_UNARY_FUNCTION(OP, FUNC, DERIVATIVE)    \
/* predefine the struct and the function for higher order derivatives*/  \
namespace codi {                                                         \
  template<typename Real, class A> struct OP;                            \
  template <typename Real, class A>                                      \
  inline  codi:: OP<Real, A> FUNC(const codi::Expression<Real, A>& a);   \
                                                                         \
  using std::FUNC;                                          \
                                                            \
  template<typename Real, class A>            \
  struct OP : public Expression<Real, OP<Real, A> > {      \
    OP(const Expression<Real, A>& a)        \
  : a_(a.cast()), result_(FUNC(a_.getValue())) { }    \
    inline void calcGradient(Real& gradient) const {      \
  a_.calcGradient(gradient, DERIVATIVE);      \
    }                \
    inline void calcGradient(Real& gradient,        \
           const Real& multiplier) const {  \
  a_.calcGradient(gradient, (DERIVATIVE)*multiplier);  \
    }                \
    inline const Real& getValue() const {      \
  return result_;            \
    }                \
  private:              \
  const A& a_;            \
  Real result_;            \
  };                \
                                                  \
  template <typename Real, class A>            \
  inline              \
  codi:: OP<Real, A> FUNC(const codi::Expression<Real, A>& a) {    \
    return codi:: OP<Real, A>(a.cast());        \
  }   \
}                \
using codi:: FUNC;


/**
 * The final argument is the value of the derivative
 */
CODI_DEFINE_UNARY_FUNCTION(Log, log, 1.0 / a_.getValue())

CODI_DEFINE_UNARY_FUNCTION(Log10, log10, 0.434294481903252 / a_.getValue())

CODI_DEFINE_UNARY_FUNCTION(Sin, sin, cos(a_.getValue()))

CODI_DEFINE_UNARY_FUNCTION(Cos, cos, -sin(a_.getValue()))

CODI_DEFINE_UNARY_FUNCTION(Asin, asin, 1.0 / sqrt(1.0 - a_.getValue() * a_.getValue()))

CODI_DEFINE_UNARY_FUNCTION(Acos, acos, -1.0 / sqrt(1.0 - a_.getValue() * a_.getValue()))

CODI_DEFINE_UNARY_FUNCTION(Atan, atan, 1.0 / (1 + a_.getValue() * a_.getValue()))

CODI_DEFINE_UNARY_FUNCTION(Sinh, sinh, cosh(a_.getValue()))

CODI_DEFINE_UNARY_FUNCTION(Cosh, cosh, sinh(a_.getValue()))

CODI_DEFINE_UNARY_FUNCTION(Abs, abs, (a_.getValue() > 0.0) - (a_.getValue() < 0.0))
//CODI_DEFINE_UNARY_FUNCTION(Abs, abs, a_.getValue() > 0 ? 1.0 : -1.0)
#undef CODI_DEFINE_UNARY_FUNCTION

/**
 * Need both fabs and abs for C compatibility, so get fabs to return
 * an Abs object
 */
template<typename Real, class A>
inline
codi::Abs<Real, A> fabs(const codi::Expression<Real, A>& a) {
  return codi::Abs<Real, A>(a.cast());
}

/**
 * Need to add ceil, floor...
 * Lots more math function in math.h: erf, bessel functions etc...
 */

template<typename Real, class A>
inline bool isinf(const codi::Expression<Real, A>& a) {
  return isinf(a.getValue());
}

template<typename Real, class A>
inline bool isnan(const codi::Expression<Real, A>& a) {
  return isnan(a.getValue());
}

using std::isfinite;
template<typename Real, class A>
inline bool isfinite(const codi::Expression<Real, A>& a) {
  return isfinite(a.getValue());
}

template<typename Real, class A>
inline bool __isinf(const codi::Expression<Real, A>& a) {
  return __isinf(a.getValue());
}

template<typename Real, class A>
inline bool
__isnan(const codi::Expression<Real, A>& a) {
  return __isnan(a.getValue());
}

template<typename Real, class A>
inline
bool
__isfinite(const codi::Expression<Real, A>& a) {
  return __isfinite(a.getValue());
}

using std::floor;
template<typename Real, class A>
inline typename codi::TypeTraits<Real>::PassiveReal floor(const codi::Expression<Real, A>& a) {
  return floor(a.getValue());
}

using std::ceil;
template<typename Real, class A>
inline typename codi::TypeTraits<Real>::PassiveReal ceil(const codi::Expression<Real, A>& a) {
  return ceil(a.getValue());
}


/**
 * Enable mathematical functions when derivative is a little more
 * complicated
 */
#define CODI_DEFINE_UNARY_FUNCTION3(OP, FUNC, DERIV1, DERIV2)            \
/* predefine the struct and the function for higher order derivatrives */\
namespace codi {                                                         \
  template<typename Real, class A> struct OP;                            \
}                                                                        \
template <typename Real, class A>                                        \
inline  codi:: OP<Real, A> FUNC(const codi::Expression<Real, A>& a);     \
                                                                         \
namespace codi {                                                         \
  template<typename Real, class A>                                       \
  struct OP : public Expression<Real, OP<Real, A> > {                    \
    OP(const Expression<Real, A>& a)                                     \
      : a_(a.cast()), result_(FUNC(a_.getValue())) { }                   \
    inline void calcGradient(Real& gradient) const {                     \
      DERIV1;                                                            \
      a_.calcGradient(gradient, DERIV2);      \
    }                \
    inline void calcGradient(Real& gradient, const Real& multiplier) const {  \
      DERIV1;              \
      a_.calcGradient(gradient, (DERIV2)*multiplier);    \
    }                \
    inline const Real& getValue() const {      \
      return result_;            \
    }                \
  private:              \
    const A& a_;            \
    Real result_;            \
  };                \
}                \
template <typename Real, class A>            \
inline  codi:: OP<Real, A> FUNC(const codi::Expression<Real, A>& a) {    \
  return codi:: OP<Real, A>(a.cast());        \
}

CODI_DEFINE_UNARY_FUNCTION3(Tan, tan, Real tmp = 1 / cos(a_.getValue()), tmp * tmp)

#undef CODI_DEFINE_UNARY_FUNCTION3

/* predefine the struct and the function for higher order derivatrives */
namespace codi {
  template<typename Real, class A, class B> struct Pow;
  template<typename Real, class A> struct PowScalarExponent;
  template<typename Real, class A> struct PowScalarBase;
}
template <typename Real, class A, class B>
inline  codi:: Pow<Real, A, B> pow(const codi::Expression<Real, A>& a, const codi::Expression<Real, B>& b);

template<typename Real, class A>
inline codi::PowScalarExponent<Real, A> pow(const codi::Expression<Real, A>& a,
                                            const typename codi::TypeTraits<Real>::PassiveReal& b);

template<typename Real, class B>
inline codi::PowScalarBase<Real, B> pow(const typename codi::TypeTraits<Real>::PassiveReal& a,
                                        const codi::Expression<Real, B>& b);

namespace codi {
  /**
   *  Pow: an expression to the power of another expression
   */
  template<typename Real, class A, class B>
  struct Pow : public Expression<Real, Pow<Real, A, B> > {
    Pow(const Expression<Real, A>& a, const Expression<Real, B>& b)
      : a_(a.cast()), b_(b.cast()),
        result_(pow(a_.getValue(), b_.getValue())) { };

    /**
     *  If f(a,b)=pow(a,b) then df/da=b*pow(a,b-1) and df/db=log(a)*pow(a,b)
     */
    inline void calcGradient(Real& gradient) const {
      a_.calcGradient(gradient, b_.getValue() * pow(a_.getValue(), b_.getValue() - 1.0));
      if (a_.getValue() > 0.0) {
        b_.calcGradient(gradient, log(a_.getValue()) * result_);
      } else {
        b_.calcGradient(gradient, 0.0);
      }
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      a_.calcGradient(gradient, b_.getValue() * pow(a_.getValue(),
                       b_.getValue() - 1.0) * multiplier);
      if (a_.getValue() > 0.0) {
        b_.calcGradient(gradient, log(a_.getValue()) * result_ * multiplier);
      } else {
        b_.calcGradient(gradient, 0.0);
      }
    }

    inline const Real& getValue() const {
      return result_;
    }

  private:
    const A& a_;
    const B& b_;
    Real result_;
  };

  /**
   *  PowScalarExponent: an expression to the power of a scalar
   */
  template<typename Real, class A>
  struct PowScalarExponent : public Expression<Real, PowScalarExponent<Real, A> > {
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    PowScalarExponent(const Expression<Real, A>& a, const PassiveReal& b)
      : a_(a.cast()), b_(b), result_(pow(a_.getValue(), b)) { }

    inline void calcGradient(Real& gradient) const {
      a_.calcGradient(gradient, b_ * pow(a_.getValue(), b_ - 1.0));
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      a_.calcGradient(gradient, b_ * pow(a_.getValue(), b_ - 1.0) * multiplier);
    }

    inline const Real& getValue() const {
      return result_;
    }

  private:
    const A& a_;
    Real b_;
    Real result_;
  };

  /**
   * PowScalarBase: a scalar to the power of an expression
   */
  template<typename Real, class B>
  struct PowScalarBase : public Expression<Real, PowScalarBase<Real, B> > {
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    PowScalarBase(const PassiveReal & a, const Expression<Real, B>& b)
      : a_(a), b_(b.cast()), result_(pow(a_, b_.getValue())) { }

    inline void calcGradient(Real& gradient) const {
      b_.calcGradient(gradient, log(a_) * result_);
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      b_.calcGradient(gradient, log(a_) * result_ * multiplier);
    }

    inline const Real& getValue() const {
      return result_;
    }

  private:
    Real a_;
    const B& b_;
    Real result_;
  };
}
/**
 * End namespace codi
 */

/**
 *  Overload pow for Expression arguments
 */
template<typename Real, class A, class B>
inline codi::Pow<Real, A, B> pow(const codi::Expression<Real, A>& a,
                                 const codi::Expression<Real, B>& b) {
  return codi::Pow<Real, A, B>(a.cast(), b.cast());
}

/**
 *  Overload pow for expression to the power of scalar
 */
template<typename Real, class A>
inline codi::PowScalarExponent<Real, A> pow(const codi::Expression<Real, A>& a,
                                            const typename codi::TypeTraits<Real>::PassiveReal& b) {
  return codi::PowScalarExponent<Real, A>(a.cast(), b);
}

/**
 * Overload pow for scalar to the power of an expression
 */
template<typename Real, class B>
inline codi::PowScalarBase<Real, B> pow(const typename codi::TypeTraits<Real>::PassiveReal& a,
                                        const codi::Expression<Real, B>& b) {
  return codi::PowScalarBase<Real, B>(a, b.cast());
}

/* predefine the struct and the function for higher order derivatrives */
namespace codi {
  template<typename Real, class A, class B> struct Atan2;
  template<typename Real, class A> struct Atan2Scalar1;
  template<typename Real, class A> struct Atan2Scalar2;
}
template <typename Real, class A, class B>
inline  codi:: Atan2<Real, A, B> atan2(const codi::Expression<Real, A>& a, const codi::Expression<Real, B>& b);

template <typename Real, class A>
inline  codi:: Atan2Scalar1<Real, A> atan2(const codi::Expression<Real, A>& a, const typename codi::TypeTraits<Real>::PassiveReal& b);

template<typename Real, class B>
inline codi::Atan2Scalar2<Real, B> atan2(const typename codi::TypeTraits<Real>::PassiveReal& a, const codi::Expression<Real, B>& b);

namespace codi {
// atan2:
  template<typename Real, class A, class B>
  struct Atan2 : public Expression<Real, Atan2<Real, A, B> > {
    Atan2(const Expression<Real, A>& a, const Expression<Real, B>& b)
      : a_(a.cast()), b_(b.cast()),
        result_(atan2(a_.getValue(), b_.getValue())) { };

    inline void calcGradient(Real& gradient) const {
      Real divisor = a_.getValue() * a_.getValue() + b_.getValue() * b_.getValue();
      divisor = 1.0 / divisor;
      a_.calcGradient(gradient, b_.getValue() * divisor);
      b_.calcGradient(gradient, -a_.getValue() * divisor);
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      Real divisor = a_.getValue() * a_.getValue() + b_.getValue() * b_.getValue();
      divisor = 1.0 / divisor;
      divisor *= multiplier;
      a_.calcGradient(gradient, b_.getValue() * divisor);
      b_.calcGradient(gradient, -a_.getValue() * divisor);
    }

    inline const Real& getValue() const {
      return result_;
    }

  private:
    const A& a_;
    const B& b_;
    Real result_;
  };

// atan2Scalar1
  template<typename Real, class A>
  struct Atan2Scalar1 : public Expression<Real, Atan2Scalar1<Real, A> > {
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    Atan2Scalar1(const Expression<Real, A>& a, const PassiveReal& b)
      : a_(a.cast()), b_(b),
        result_(atan2(a_.getValue(), b)) { }

    inline void calcGradient(Real& gradient) const {
      Real divisor = a_.getValue() * a_.getValue() + b_ * b_;
      divisor = 1.0 / divisor;
      a_.calcGradient(gradient, b_ * divisor);
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      Real divisor = a_.getValue() * a_.getValue() + b_ * b_;
      divisor = 1.0 / divisor;
      divisor *= multiplier;
      a_.calcGradient(gradient, b_ * divisor);
    }

    inline const Real& getValue() const {
      return result_;
    }

  private:
    const A& a_;
    Real b_;
    Real result_;
  };

// atan2Scalar2
  template<typename Real, class B>
  struct Atan2Scalar2 : public Expression<Real, Atan2Scalar2<Real, B> > {
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    Atan2Scalar2(const PassiveReal& a, const Expression<Real, B>& b)
      : a_(a), b_(b.cast()),
        result_(atan2(a_, b_.getValue())) { }

    inline void calcGradient(Real& gradient) const {
      Real divisor = a_ * a_ + b_.getValue() * b_.getValue();
      divisor = 1.0 / divisor;
      b_.calcGradient(gradient, -a_ * divisor);
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      Real divisor = a_ * a_ + b_.getValue() * b_.getValue();
      divisor = 1.0 / divisor;
      divisor *= multiplier;
      b_.calcGradient(gradient, -a_ * divisor);
    }

    inline const Real& getValue() const {
      return result_;
    }

  private:
    Real a_;
    const B& b_;
    Real result_;
  };
} // End namespace codi

// Overload pow for Expression arguments
template<typename Real, class A, class B>
inline
codi::Atan2<Real, A, B> atan2(const codi::Expression<Real, A>& a,
                              const codi::Expression<Real, B>& b) {
  return codi::Atan2<Real, A, B>(a.cast(), b.cast());
}

// Overload pow for expression to the power of scalar
template<typename Real, class A>
inline
codi::Atan2Scalar1<Real, A> atan2(const codi::Expression<Real, A>& a,
                                  const typename codi::TypeTraits<Real>::PassiveReal& b) {
  return codi::Atan2Scalar1<Real, A>(a.cast(), b);
}

// Overload pow for scalar to the power of an expression
template<typename Real, class B>
inline codi::Atan2Scalar2<Real, B> atan2(const typename codi::TypeTraits<Real>::PassiveReal& a,
                                         const codi::Expression<Real, B>& b) {
  return codi::Atan2Scalar2<Real, B>(a, b.cast());
}

