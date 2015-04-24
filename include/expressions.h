#pragma once

#include <cmath>

#include "configure.h"

namespace codi {

  /**
   * The Expression type from which all other types of expression
   * derive. Each member function simply calls the specialized version
   * of the function according to the expression's true type, which is
   * given by its template argument.
   */

  template<class A>
  struct Expression {

    /**
     * Cast the expression to its true type, given by the template
     * argument
     */
    const A& cast() const {
      return static_cast<const A&>(*this);
    }

    /**
     * Calculate the gradient of the mathematical operation that this
     * expression represents and pass the result to its argument.
     * For functions f(a), pass df/da to the argument in the
     * first case and pass multiplier*df/da in the second case.
     */
    void calcGradient(Real& gradient) const {
      cast().calcGradient(gradient);
    }

    /**
     * As the previous but multiplying the gradient by "multiplier"
     */
    void calcGradient(Real& gradient, const Real& multiplier) const {
      cast().calcGradient(gradient, multiplier);
    }

    /**
     * Return the numerical value of the expression
     */
    Real getValue() const {
      return cast().getValue();
    }

    /**
     * Calculate the gradient and return the numerical value
     */
    Real value_and_gradient(Real& gradient) const {
      Real a = 0.0;
      cast().calcGradient(a);
      gradient = a;
      return cast().getValue();
    }

  private:
    /**
     * Intentionally inaccessible to prevent an expression appearing
     * on the left-hand-side of a statement
     */
    Expression& operator=(const Expression&) {
      return *this;
    }
  };

  /**
   * Now define particular types of expression, using static
   * polymorphism via the Curiously Recurring Template Pattern
   */

  /**
   * Add: an expression plus another expression
   */
  template<class A, class B>
  struct Add : public Expression<Add<A, B> > {
    Add(const Expression<A>& a, const Expression<B>& b)
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
  template<class A, class B>
  inline
  Add<A, B> operator+(const Expression<A>& a,
                      const Expression<B>& b) {
    return Add<A, B>(a.cast(), b.cast());
  }

  /**
   * Subtract: an expression minus another expression
   */
  template<class A, class B>
  struct Subtract : public Expression<Subtract<A, B> > {
    Subtract(const Expression<A>& a, const Expression<B>& b)
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
  template<class A, class B>
  inline
  Subtract<A, B> operator-(const Expression<A>& a,
                           const Expression<B>& b) {
    return Subtract<A, B>(a.cast(), b.cast());
  }

  /**
   * Multiply: an expression multiplied by another expression
   */
#ifdef CODI_MULTIPLY_PRECOMPUTES_RESULT
  /**
   * The first version precomputes the result, which should be optimal
   * if getValue() is called multiple times
   */
  template <class A, class B>
  struct Multiply : public Expression<Multiply<A,B> > {
    Multiply(const Expression<A>& a, const Expression<B>& b)
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
    inline Real getValue() const {
      return result_;
    }
  private:
    const A& a_;
    const B& b_;
    Real result_;
  };
#else

  /**
   * The second version does not precompute the result
   */
  template<class A, class B>
  struct Multiply : public Expression<Multiply<A, B> > {
    Multiply(const Expression<A>& a, const Expression<B>& b)
      : a_(a.cast()), b_(b.cast()) { }

    // If f(a,b) = a*b then df/da = b and df/db = a
    inline void calcGradient(Real& gradient) const {
      a_.calcGradient(gradient, b_.getValue());
      b_.calcGradient(gradient, a_.getValue());
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      a_.calcGradient(gradient, b_.getValue() * multiplier);
      b_.calcGradient(gradient, a_.getValue() * multiplier);
    }

    inline Real getValue() const {
      return a_.getValue() * b_.getValue();
    }

  private:
    const A& a_;
    const B& b_;
  };

#endif

  /**
   * Overload multiplication operator for Expression arguments
   */
  template<class A, class B>
  inline
  Multiply<A, B> operator*(const Expression<A>& a,
                           const Expression<B>& b) {
    return Multiply<A, B>(a.cast(), b.cast());
  }

  /**
   *  Divide: an expression divided by another expression
   */
  template<class A, class B>
  struct Divide : public Expression<Divide<A, B> > {
    Divide(const Expression<A>& a, const Expression<B>& b)
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

    inline Real getValue() const {
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
  template<class A, class B>
  inline
  Divide<A, B> operator/(const Expression<A>& a,
                         const Expression<B>& b) {
    return Divide<A, B>(a.cast(), b.cast());
  }

  /**
   * ScalarAdd: an expression plus a scalar
   */
  template<class A>
  struct ScalarAdd : public Expression<ScalarAdd<A> > {
    ScalarAdd(const Expression<A>& a, const Real& b)
      : a_(a.cast()), result_(a_.getValue() + b) { }

    inline void calcGradient(Real& gradient) const {
      a_.calcGradient(gradient);
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      a_.calcGradient(gradient, multiplier);
    }

    inline Real getValue() const {
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
  template<class A>
  inline
  ScalarAdd<A> operator+(const Expression<A>& a,
                         const Real& b) {
    return ScalarAdd<A>(a.cast(), b);
  }

  template<class A>
  inline
  ScalarAdd<A> operator+(const Real& b,
                         const Expression<A>& a) {
    return ScalarAdd<A>(a.cast(), b);
  }

  /**
   * Overload subtraction operator for expression minus scalar to
   * return a ScalarAdd object
   */
  template<class A>
  inline
  ScalarAdd<A> operator-(const Expression<A>& a,
                         const Real& b) {
    return ScalarAdd<A>(a.cast(), -b);
  }

  /**
   *  ScalarSubtract: scalar minus expression
   */
  template<class B>
  struct ScalarSubtract : public Expression<ScalarSubtract<B> > {
    ScalarSubtract(const Real& a, const Expression<B>& b)
      : b_(b.cast()), result_(a - b_.getValue()) { }

    inline void calcGradient(Real& gradient) const {
      b_.calcGradient(gradient, -1.0);
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      b_.calcGradient(gradient, -multiplier);
    }

    inline Real getValue() const {
      return result_;
    }

  private:
    const B& b_;
    Real result_;
  };

  /**
   *  Overload subtraction operator for scalar minus expression
   */
  template<class B>
  inline
  ScalarSubtract<B> operator-(const Real& a,
                              const Expression<B>& b) {
    return ScalarSubtract<B>(a, b.cast());
  }

  /**
   *  ScalarMultiply: expression multiplied by scalar
   */
  template<class A>
  struct ScalarMultiply
      : public Expression<ScalarMultiply<A> > {
    ScalarMultiply(const Expression<A>& a, const Real& b)
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
    Real b_;
  };

  /**
   * Overload multiplication operator for expression multiplied by
   * scalar and scalar multiplied by expression
   */
  template<class A>
  inline
  ScalarMultiply<A> operator*(const Expression<A>& a,
                              const Real& b) {
    return ScalarMultiply<A>(a.cast(), b);
  }

  template<class A>
  inline
  ScalarMultiply<A> operator*(const Real& b,
                              const Expression<A>& a) {
    return ScalarMultiply<A>(a.cast(), b);
  }

  /**
   *  Overload division operator for expression divided by scalar to
   * return ScalarMultiply object
   */
  template<class A>
  inline
  ScalarMultiply<A> operator/(const Expression<A>& a,
                              const Real& b) {
    return ScalarMultiply<A>(a.cast(), 1.0 / b);
  }

  /**
   * ScalarDivide: scalar divided by expression
   */
  template<class B>
  struct ScalarDivide : public Expression<ScalarDivide<B> > {
    ScalarDivide(const Real& a, const Expression<B>& b)
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

    inline Real getValue() const {
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
  template<class B>
  inline
  ScalarDivide<B> operator/(const Real& a,
                            const Expression<B>& b) {
    return ScalarDivide<B>(a, b.cast());
  }

  /**
   *  Conditional operators should behave exactly the same as with
   * non-active arguments so in each of the cases below the getValue()
   * function is called to extract the value of the expression
   */
#define CODI_DEFINE_CONDITIONAL(OPERATOR, OP)      \
template <class A, class B>          \
inline              \
bool OPERATOR(const Expression<A>& a,        \
      const Expression<B>& b) {      \
  return a.getValue() OP b.getValue();        \
}                \
                              \
template <class A>            \
inline              \
bool OPERATOR(const Expression<A>& a, const Real& b) {  \
  return a.getValue() OP b;          \
}                \
                              \
template <class B>            \
inline              \
bool OPERATOR(const Real& a, const Expression<B>& b) {  \
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
  template<class A>
  struct UnaryMinus : public Expression<UnaryMinus<A> > {
    UnaryMinus(const Expression<A>& a)
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
  template<class A>
  inline
  UnaryMinus<A> operator-(const Expression<A>& a) {
    return UnaryMinus<A>(a.cast());
  }

  /**
   *  Unary plus: returns the argument
   */
  template<class A>
  inline
  A operator+(const Expression<A>& a) {
    return a;
  }

  /**
   *  Exponential of an expression
   */
  template<class A>
  struct Exp : public Expression<Exp<A> > {
    Exp(const Expression<A>& a)
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

    inline Real getValue() const {
      return result_;
    }

  private:
    const A& a_;
    Real result_;
  };

  /**
   *  Exponential of an expression
   */
  template<class A>
  struct Atanh : public Expression<Atanh<A> > {
    Atanh(const Expression<A>& a)
      : a_(a.cast()), result_(.5 * std::log((1 + a.getValue()) / (1 - a.getValue()))) { }

    inline void calcGradient(Real& gradient) const {
      a_.calcGradient(gradient, 1.0 / (1 - a_.getValue() * a_.getValue()));
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      a_.calcGradient(gradient, 1.0 / (1 - a_.getValue() * a_.getValue()) * multiplier);
    }

    inline Real getValue() const {
      return result_;
    }

  private:
    const A& a_;
    Real result_;
  };

}
/**
 * End namespace codi
 */

/**
 *  It is important to place overloads of mathematical functions in the
 * global namespace.  If exp(Expression) was placed in the codi
 * namespace then the Exp class (which is in the codi namespace would
 * not be able to find the std::exp(double) function due to C++ name
 * look-up rules.
 */

/*
 *  Overload exp for Expression objects
 */
template<class A>
inline
codi::Exp<A> exp(const codi::Expression<A>& a) {
  return codi::Exp<A>(a.cast());
}

template<class A>
inline
codi::Atanh<A> atanh(const codi::Expression<A>& a) {
  return codi::Atanh<A>(a.cast());
}


/**
 * Enable unary mathematical functions when the derivative is most
 * easily written in terms of the argument of the function; note that
 * the class is in the codi name space but the function is not.
 */
#define CODI_DEFINE_UNARY_FUNCTION(OP, FUNC, DERIVATIVE)    \
namespace codi {            \
  template <class A>            \
  struct OP : public Expression<OP<A> > {      \
    OP(const Expression<A>& a)        \
  : a_(a.cast()), result_(std::FUNC(a_.getValue())) { }    \
    inline void calcGradient(Real& gradient) const {      \
  a_.calcGradient(gradient, DERIVATIVE);      \
    }                \
    inline void calcGradient(Real& gradient,        \
           const Real& multiplier) const {  \
  a_.calcGradient(gradient, (DERIVATIVE)*multiplier);  \
    }                \
    inline Real getValue() const {      \
  return result_;            \
    }                \
  private:              \
  const A& a_;            \
  Real result_;            \
  };                \
}                \
template <class A>            \
inline              \
codi:: OP<A> FUNC(const codi::Expression<A>& a) {    \
  return codi:: OP<A>(a.cast());        \
}

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
template<class A>
inline
codi::Abs<A> fabs(const codi::Expression<A>& a) {
  return codi::Abs<A>(a.cast());
}

/**
 * Need to add ceil, floor...
 * Lots more math function in math.h: erf, bessel functions etc...
 */

template<class A>
inline
bool
isinf(const codi::Expression<A>& a) {
  return std::isinf(a.getValue());
}

template<class A>
inline
bool
isnan(const codi::Expression<A>& a) {
  return std::isnan(a.getValue());
}

template<class A>
inline
bool
isfinite(const codi::Expression<A>& a) {
  return std::isfinite(a.getValue());
}

template<class A>
inline
bool
__isinf(const codi::Expression<A>& a) {
  return __isinf(a.getValue());
}

template<class A>
inline
bool
__isnan(const codi::Expression<A>& a) {
  return __isnan(a.getValue());
}

template<class A>
inline
bool
__isfinite(const codi::Expression<A>& a) {
  return __isfinite(a.getValue());
}

template<class A>
inline Real floor(const codi::Expression<A>& a) {
  return std::floor(a.getValue());
}


/**
 * Enable mathematical functions when the derivative is most easily
 * written in terms of the return value of the function
 */
#define CODI_DEFINE_UNARY_FUNCTION2(OP, FUNC, DERIVATIVE)  \
namespace codi {            \
  template <class A>            \
  struct OP : public Expression<OP<A> > {      \
    OP(const Expression<A>& a)        \
  : a_(a.cast()), result_(FUNC(a.getValue())) { }    \
    inline void calcGradient(Real& gradient) const {      \
  a_.calcGradient(gradient, DERIVATIVE);      \
    }                \
    inline void calcGradient(Real& gradient,        \
           const Real& multiplier) const {  \
  a_.calcGradient(gradient, (DERIVATIVE)*multiplier);  \
    }                \
    inline Real getValue() const {      \
  return result_;            \
    }                \
  private:              \
  const A& a_;            \
  Real result_;            \
  };                \
}                \
template <class A>            \
inline              \
codi:: OP<A> FUNC(const codi::Expression<A>& a) {    \
  return codi:: OP<A>(a.cast());        \
}

CODI_DEFINE_UNARY_FUNCTION2(Sqrt, sqrt, result_ != 0 ? 0.5 / result_ : 0.0)

CODI_DEFINE_UNARY_FUNCTION2(Tanh, tanh, 1 - result_ * result_)
#undef CODI_DEFINE_UNARY_FUNCTION2


/**
 * Enable mathematical functions when derivative is a little more
 * complicated
 */
#define CODI_DEFINE_UNARY_FUNCTION3(OP, FUNC, DERIV1, DERIV2)  \
namespace codi {            \
  template <class A>            \
  struct OP : public Expression<OP<A> > {      \
    OP(const Expression<A>& a)        \
  : a_(a.cast()), result_(FUNC(a_.getValue())) { }    \
    inline void calcGradient(Real& gradient) const {      \
  DERIV1;              \
  a_.calcGradient(gradient, DERIV2);      \
    }                \
    inline void calcGradient(Real& gradient,        \
           const Real& multiplier) const {  \
  DERIV1;              \
  a_.calcGradient(gradient, (DERIV2)*multiplier);    \
    }                \
    inline Real getValue() const {      \
  return result_;            \
    }                \
  private:              \
  const A& a_;            \
  Real result_;            \
  };                \
}                \
template <class A>            \
inline              \
codi:: OP<A> FUNC(const codi::Expression<A>& a) {    \
  return codi:: OP<A>(a.cast());        \
}

CODI_DEFINE_UNARY_FUNCTION3(Tan, tan, Real tmp = 1 / cos(a_.getValue()), tmp * tmp)
//CODI_DEFINE_UNARY_FUNCTION3(Tanh, tanh, Real e=exp(2.0*a_.getValue()), 4.0*e/((2.0+2)*e+1.0))

#undef CODI_DEFINE_UNARY_FUNCTION3

namespace codi {
  /**
   *  Pow: an expression to the power of another expression
   */
  template<class A, class B>
  struct Pow : public Expression<Pow<A, B> > {
    Pow(const Expression<A>& a, const Expression<B>& b)
      : a_(a.cast()), b_(b.cast()),
        result_(pow(a_.getValue(), b_.getValue())) { };

    /**
     *  If f(a,b)=pow(a,b) then df/da=b*pow(a,b-1) and df/db=log(a)*pow(a,b)
     */
    inline void calcGradient(Real& gradient) const {
      a_.calcGradient(gradient, b_.getValue() * pow(a_.getValue(), b_.getValue() - 1.0));
      if (a_.getValue() > 0) {
        b_.calcGradient(gradient, log(a_.getValue()) * result_);
      } else {
        b_.calcGradient(gradient, 0.0);
      }
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      a_.calcGradient(gradient, b_.getValue() * pow(a_.getValue(),
                       b_.getValue() - 1.0) * multiplier);
      if (a_.getValue() > 0) {
        b_.calcGradient(gradient, log(a_.getValue()) * result_ * multiplier);
      } else {
        b_.calcGradient(gradient, 0.0);
      }
    }

    inline Real getValue() const {
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
  template<class A>
  struct PowScalarExponent : public Expression<PowScalarExponent<A> > {
    PowScalarExponent(const Expression<A>& a, const Real& b)
      : a_(a.cast()), b_(b), result_(pow(a_.getValue(), b)) { }

    inline void calcGradient(Real& gradient) const {
      a_.calcGradient(gradient, b_ * pow(a_.getValue(), b_ - 1.0));
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      a_.calcGradient(gradient, b_ * pow(a_.getValue(), b_ - 1.0) * multiplier);
    }

    inline Real getValue() const {
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
  template<class B>
  struct PowScalarBase : public Expression<PowScalarBase<B> > {
    PowScalarBase(const Real& a, const Expression<B>& b)
      : a_(a), b_(b.cast()), result_(pow(a_, b_.getValue())) { }

    inline void calcGradient(Real& gradient) const {
      b_.calcGradient(gradient, log(a_) * result_);
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      b_.calcGradient(gradient, log(a_) * result_ * multiplier);
    }

    inline Real getValue() const {
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
template<class A, class B>
inline
codi::Pow<A, B> pow(const codi::Expression<A>& a,
                    const codi::Expression<B>& b) {
  return codi::Pow<A, B>(a.cast(), b.cast());
}

/**
 *  Overload pow for expression to the power of scalar
 */
template<class A>
inline
codi::PowScalarExponent<A> pow(const codi::Expression<A>& a,
                               const Real& b) {
  return codi::PowScalarExponent<A>(a.cast(), b);
}

/**
 * Overload pow for scalar to the power of an expression
 */
template<class B>
inline
codi::PowScalarBase<B> pow(const Real& a,
                           const codi::Expression<B>& b) {
  return codi::PowScalarBase<B>(a, b.cast());
}

namespace codi {
// atan2: 
  template<class A, class B>
  struct Atan2 : public Expression<Atan2<A, B> > {
    Atan2(const Expression<A>& a, const Expression<B>& b)
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

    inline Real getValue() const {
      return result_;
    }

  private:
    const A& a_;
    const B& b_;
    Real result_;
  };

// atan2Scalar1
  template<class A>
  struct Atan2Scalar1 : public Expression<Atan2Scalar1<A> > {
    Atan2Scalar1(const Expression<A>& a, const Real& b)
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

    inline Real getValue() const {
      return result_;
    }

  private:
    const A& a_;
    Real b_;
    Real result_;
  };

// atan2Scalar2
  template<class B>
  struct Atan2Scalar2 : public Expression<Atan2Scalar2<B> > {
    Atan2Scalar2(const Real& a, const Expression<B>& b)
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

    inline Real getValue() const {
      return result_;
    }

  private:
    Real a_;
    const B& b_;
    Real result_;
  };
} // End namespace codi

// Overload pow for Expression arguments
template<class A, class B>
inline
codi::Atan2<A, B> atan2(const codi::Expression<A>& a,
                        const codi::Expression<B>& b) {
  return codi::Atan2<A, B>(a.cast(), b.cast());
}

// Overload pow for expression to the power of scalar
template<class A>
inline
codi::Atan2Scalar1<A> atan2(const codi::Expression<A>& a,
                            const Real& b) {
  return codi::Atan2Scalar1<A>(a.cast(), b);
}

// Overload pow for scalar to the power of an expression
template<class B>
inline
codi::Atan2Scalar2<B> atan2(const Real& a,
                            const codi::Expression<B>& b) {
  return codi::Atan2Scalar2<B>(a, b.cast());
}

