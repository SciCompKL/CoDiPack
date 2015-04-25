#pragma once

#include "activeReal.hpp"
#include "tapes/forwardEvaluation.hpp"

namespace codi {
  typedef ActiveReal<double, ForwardEvaluation<double> > RealForward;
  typedef ActiveReal<float, ForwardEvaluation<float> > RealForwardFloat;
}