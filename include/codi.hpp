#pragma once

#include "activeReal.hpp"
#include "tapes/forwardEvaluation.hpp"

namespace codi {
  typedef ActiveReal<ForwardEvaluation> RealForward;
}