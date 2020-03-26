#pragma once

#include "codi/config.h"

#include "codi/expressions/real/allOperators.hpp"
#include "codi/expressions/activeType.hpp"
#include "codi/tapes/forwardEvaluation.hpp"
#include "codi/tapes/indices/linearIndexManager.hpp"
#include "codi/tapes/indices/multiUseIndexManager.hpp"
#include "codi/tapes/jacobianLinearTape.hpp"
#include "codi/tapes/jacobianReuseTape.hpp"
#include "codi/tapes/primalValueLinearTape.hpp"
#include "codi/tapes/primalValueReuseTape.hpp"
#include "codi/tapes/statementEvaluators/reverseStatementEvaluator.hpp"
#include "codi/tapes/statementEvaluators/directStatementEvaluator.hpp"
#include "codi/tapes/statementEvaluators/innerStatementEvaluator.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename Real, typename Gradient = Real>
  using RealForwardGen = ActiveType<ForwardEvaluation<Real, Gradient>>;

  using RealForward = RealForwardGen<double, double>;


  template<typename Real, typename Index, typename Gradient = Real>
  using RealReverseGen = ActiveType<JacobianLinearTape<JacobianTapeTypes<Real, Gradient, LinearIndexManager<Index>>>>;

  using RealReverse = RealReverseGen<double, int, double>;

  template<typename Real, typename IndexManager, typename Gradient = Real>
  using RealReverseIndexGen = ActiveType<JacobianReuseTape<JacobianTapeTypes<Real, Gradient, IndexManager>>>;

  using RealReverseIndex = RealReverseIndexGen<double, MultiUseIndexManager<int>, double>;

  template<typename Real, typename IndexManager, template <typename> class StatementEvaluator, typename Gradient = Real>
  using RealReversePrimalGen = ActiveType<PrimalValueLinearTape<PrimalValueTapeTypes<Real, Gradient, IndexManager, StatementEvaluator>>>;

  using RealReversePrimal = RealReversePrimalGen<double, LinearIndexManager<int>, ReverseStatementEvaluator, double>;

  template<typename Real, typename IndexManager, template <typename> class StatementEvaluator, typename Gradient = Real>
  using RealReversePrimalIndexGen = ActiveType<PrimalValueReuseTape<PrimalValueTapeTypes<Real, Gradient, IndexManager, StatementEvaluator>>>;

  using RealReversePrimalIndex = RealReversePrimalIndexGen<double, MultiUseIndexManager<int>, ReverseStatementEvaluator, double>;

}
