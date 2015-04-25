#pragma once

#include "tapeInterface.hpp"

namespace codi {
  /**
   * @brief Tape for the tangent or forward AD mode
   *
   * This tape implements the forward or tangent AD mode. For each expression the equation
   *
   * \dot y = \frac{df}{dx}(x)\dot x
   *
   * is evaluated. The gradient data type of this tape is just the same as the active type
   * uses for the storage of the floating point values.
   */
  class ForwardEvaluation : public TapeInterface<Real >{
  public:

    /**
     * @brief Evaluates the primal expression and the tangent
     *
     * The store method evaluates the forward AD equation and the primal equation.
     *
     * @param[out]      value  The value of the rhs.
     * @param[out] lhsTangent  The tangent of the lhs.
     * @param[in]         rhs  The expression of the rhs.
     */
    template<typename Rhs>
    inline void store(Real& value, GradientData& lhsTangent, const Rhs& rhs) {
      lhsTangent = 0.0;
      rhs.calcGradient(lhsTangent);
      value = rhs.getValue();
    }

    /**
     * @brief Adds the jacobi to the tangent value of the expression.
     *
     * This method is called for each value on the rhs. The tangent of the value is added to the
     * tangent of the lhs.
     *
     * @param[in/out] lhsTangent  The tangent of the lhs.
     * @param[]            value  Not used
     * @param[in]     curTangent  The tangent of the current rhs value.
     */
    inline void pushJacobi(Real& lhsTangent, const Real& /*value*/, const GradientData& curTangent) {
      lhsTangent += curTangent;
    }

    /**
     * @brief Adds the jacobi to the tangent value of the expression.
     *
     * This method is called for each value on the rhs. The tangent of the value times the jacobi is added to the
     * tangent of the lhs.
     *
     * @param[in/out] lhsTangent  The tangent of the lhs.
     * @param[in]         jacobi  The jacobi value of the operation.
     * @param[]            value  Not used
     * @param[in]     curTangent  The tangent of the current rhs value.
     */
    inline void pushJacobi(Real& lhsTangent, const Real& jacobi, const Real& /*value*/, const GradientData& curTangent) {
      lhsTangent += jacobi * curTangent;
    }

    /**
     * @brief Tangent is set to zero.
     *
     * The tangent is initialized with zero.
     *
     * @param[]      value  Not used.
     * @param[out] tangent  Set to zero.
     */
    inline void initGradientData(Real& /*value*/, GradientData& tangent) {
      tangent = Real(0.0);
    }

    /**
     * @brief Nothing to do.
     */
    inline void destroyGradientData(Real& /*value*/, GradientData& /*tangent*/) {
      /* do nothing */
    }

    /**
     * Sets the gradient data to the tangent value.
     *
     * @param[out]   tangent  The tangent value of the active type.
     * @param[in] newTangent  The new tangent value.
     */
    void setGradient(GradientData& tangent, const Real& newTangent) {
      tangent = newTangent;
    }

    /**
     * Returns the tangent value of the active type.
     *
     * @param[in]  tangent  The gradient data of the active type is the tangent.
     *
     * @return The tangent value of the active type.
     */
    Real getGradient(const GradientData& tangent) const {
      return tangent;
    }
  };

  /**
   * @brief Specialization for store which has a constant value on the rhs.
   *
   * This implementation of store sets the gradient of th active type to zero as the rhs
   * is inactive.
   */
  template<>
  inline void ForwardEvaluation::store<Real>(Real& value, GradientData& tangent, const Real& rhs) {
    tangent = Real(0.0);
    value = rhs;
  }
}


