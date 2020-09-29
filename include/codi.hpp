#pragma once

#include "codi/config.h"

#include "codi/expressions/real/allOperators.hpp"
#include "codi/expressions/activeType.hpp"
#include "codi/expressions/referenceActiveType.hpp"
#include "codi/tapes/data/blockData.hpp"
#include "codi/tapes/data/chunkedData.hpp"
#include "codi/tapes/forwardEvaluation.hpp"
#include "codi/tapes/indices/linearIndexManager.hpp"
#include "codi/tapes/indices/multiUseIndexManager.hpp"
#include "codi/tapes/jacobianLinearTape.hpp"
#include "codi/tapes/jacobianReuseTape.hpp"
#include "codi/tapes/primalValueLinearTape.hpp"
#include "codi/tapes/primalValueReuseTape.hpp"
#include "codi/tapes/statementEvaluators/directStatementEvaluator.hpp"
#include "codi/tapes/statementEvaluators/innerStatementEvaluator.hpp"
#include "codi/tapes/statementEvaluators/reverseStatementEvaluator.hpp"
#include "codi/tools/data/externalFunctionUserData.hpp"
#include "codi/tools/helpers/customAdjointVectorHelper.hpp"
#include "codi/tools/helpers/externalFunctionHelper.hpp"
#include "codi/tools/helpers/statementPushHelper.hpp"
#include "codi/tools/higherOrderAccess.hpp"
#include "codi/traits/numericLimits.hpp"
#include "codi/traits/tapeTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * \section sec_forwardAD Forward AD Equation
   *
   * \f$ \dot w = dphi/du * \dot u \f$
   *
   * \section sec_reverseAD Reverse AD Equation
   *
   * \f$ \bar u = dphi/du^T * \bar w \f$
   *
   */
  struct Temp {};

  /* TODO names */
  template<typename Chunk, typename NestedData = EmptyData>
  using DefaultBlockData = BlockData<Chunk, NestedData>;

  template<typename Chunk, typename NestedData = EmptyData>
  using DefaultChunkedData = ChunkedData<Chunk, NestedData>;

  template<typename Real, typename Gradient = Real>
  using RealForwardGen = ActiveType<ForwardEvaluation<Real, Gradient>>;

  using RealForward = RealForwardGen<double, double>;


  template<typename Real, typename Gradient = Real, typename Index = int>
  using RealReverseGen = ActiveType<JacobianLinearTape<JacobianTapeTypes<Real, Gradient, LinearIndexManager<Index>, DefaultChunkedData>>>;

  using RealReverse = RealReverseGen<double>;

  template<typename Real, typename Gradient = Real, typename Index = int>
  using RealReverseUncheckedGen = ActiveType<JacobianLinearTape<JacobianTapeTypes<Real, Gradient, LinearIndexManager<Index>, DefaultBlockData>>>;

  using RealReverseUnchecked = RealReverseUncheckedGen<double>;

  template<typename Real, typename Gradient = Real, typename IndexManager = MultiUseIndexManager<int>>
  using RealReverseIndexGen = ActiveType<JacobianReuseTape<JacobianTapeTypes<Real, Gradient, IndexManager, DefaultChunkedData>>>;

  using RealReverseIndex = RealReverseIndexGen<double>;

  template<typename Real, typename Gradient = Real, typename IndexManager = MultiUseIndexManager<int>>
  using RealReverseIndexUncheckedGen = ActiveType<JacobianReuseTape<JacobianTapeTypes<Real, Gradient, IndexManager, DefaultChunkedData>>>;

  using RealReverseIndexUnchecked = RealReverseIndexUncheckedGen<double>;

  template<typename Real, typename Gradient = Real, typename Index = int, template <typename> class StatementEvaluator = InnerStatementEvaluator>
  using RealReversePrimalGen = ActiveType<PrimalValueLinearTape<PrimalValueTapeTypes<Real, Gradient, LinearIndexManager<Index>, StatementEvaluator, DefaultChunkedData>>>;

  using RealReversePrimal = RealReversePrimalGen<double>;

  template<typename Real, typename Gradient = Real, typename Index = int, template <typename> class StatementEvaluator = InnerStatementEvaluator>
  using RealReversePrimalUncheckedGen = ActiveType<PrimalValueLinearTape<PrimalValueTapeTypes<Real, Gradient, LinearIndexManager<Index>, StatementEvaluator, DefaultChunkedData>>>;

  using RealReversePrimalUnchecked = RealReversePrimalUncheckedGen<double>;

  template<typename Real, typename Gradient = Real, typename IndexManager = MultiUseIndexManager<int>, template <typename> class StatementEvaluator = InnerStatementEvaluator>
  using RealReversePrimalIndexGen = ActiveType<PrimalValueReuseTape<PrimalValueTapeTypes<Real, Gradient, IndexManager, StatementEvaluator, DefaultChunkedData>>>;

  using RealReversePrimalIndex = RealReversePrimalIndexGen<double>;

  template<typename Real, typename Gradient = Real, typename IndexManager = MultiUseIndexManager<int>, template <typename> class StatementEvaluator = InnerStatementEvaluator>
  using RealReversePrimalIndexUncheckedGen = ActiveType<PrimalValueReuseTape<PrimalValueTapeTypes<Real, Gradient, IndexManager, StatementEvaluator, DefaultChunkedData>>>;

  using RealReversePrimalIndexUnchecked = RealReversePrimalIndexUncheckedGen<double>;

}
