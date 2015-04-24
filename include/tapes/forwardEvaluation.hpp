#pragma once

#include "tapeInterface.hpp"

namespace codi {
  class ForwardEvaluation : public TapeInterface<Real >{
  public:

    template<typename Rhs>
    inline void store(Real& value, GradientData& gradientData, const Rhs& rhs) {
      rhs.calcGradient(gradientData);
      value = rhs.getValue();
    }

    inline void pushJacobi(Real& gradient, const Real& /*value*/, const GradientData& gradientData) {
      gradient += gradientData;
    }

    inline void pushJacobi(Real& gradient, const Real& jacobi, const Real& /*value*/, const GradientData& gradientData) {
      gradient += jacobi * gradientData;
    }

    inline void initGradientData(Real& /*value*/, GradientData& gradientData) {
      gradientData = Real(0.0);
    }

    void setGradient(GradientData& gradientData, const Real& gradient) {
      gradientData = gradient;
    }
    Real getGradient(const GradientData& gradientData) const {
      return gradientData;
    }
  };

  template<>
  inline void ForwardEvaluation::store<Real>(Real& value, GradientData& gradientData, const Real& rhs) {
    value = rhs;
    gradientData = Real(0.0);
  }

}


