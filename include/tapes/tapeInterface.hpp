#pragma once

#include "../configure.h"

namespace codi {
  template <typename GradientDataType>
  class TapeInterface {
  public:

    typedef GradientDataType GradientData;

    /**
     * This functions are called from the expression templates. They tell the
     * tape about the operations which is evaluated.
     */

    /**
     * @brief Set the value of the expression and evaluate the gradient.
     *
     * This function is called for every operation which resembles
     *
     *   lhs = rhs;
     *
     * Store sets the primal value of the operation which is stored and the tape can
     * do the evaluation of the gradient and perform other appropriate steps.
     *
     * @param[out]          lhsValue    The primal value of the lhs. This value is set to the value
     *                                  of the right hand side.
     * @param[out]   lhsGradientData    The gradient data of the lhs. The tape will update the gradient data
     *                                  according to the expression on the right hand side.
     * @param[in]                rhs    The right hand side expression of the assignment.
     */
    template<typename Rhs>
    void store(Real& lhsValue, GradientData& lhsGradientData, const Rhs& rhs);

    virtual void pushJacobi(Real& gradient, const Real& value, const GradientData& gradientData) = 0;

    virtual void pushJacobi(Real& gradient, const Real& jacobi, const Real& value, const GradientData& gradientData) = 0;

    virtual void initGradientData(Real& value, GradientData& gradientData) = 0;

    /**
     * Access functions for the gradient information.
     */
    virtual void setGradient(GradientData& value, const Real& gradient) = 0;
    virtual Real getGradient(const GradientData& value) const = 0;

  };
}