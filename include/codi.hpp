#pragma once

#include "codi/config.h"

#include "codi/expressions/real/allOperators.hpp"
#include "codi/expressions/activeType.hpp"
#include "codi/expressions/referenceActiveType.hpp"
#include "codi/tapes/data/blockVector.hpp"
#include "codi/tapes/data/chunkVector.hpp"
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
#include "codi/tools/data/externalFunctionData.hpp"
#include "codi/tools/helpers/customAdjointVectorHelper.hpp"
#include "codi/tools/helpers/externalFunctionHelper.hpp"
#include "codi/tools/helpers/statementPushHelper.hpp"
#include "codi/tools/higherOrderAccess.hpp"
#include "codi/traits/numericLimits.hpp"
#include "codi/traits/tapeTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /* TODO names */
  template<typename Chunk, typename NestedVector = EmptyVector>
  using BlockVector = BlockVectorImpl<Chunk, NestedVector>;

  template<typename Chunk, typename NestedVector = EmptyVector>
  using ChunkVector = ChunkVectorImpl<Chunk, NestedVector>;

  template<typename Real, typename Gradient = Real>
  using RealForwardGen = ActiveType<ForwardEvaluation<Real, Gradient>>;

  using RealForward = RealForwardGen<double, double>;


  template<typename Real, typename Gradient = Real, typename Index = int>
  using RealReverseGen = ActiveType<JacobianLinearTape<JacobianTapeTypes<Real, Gradient, LinearIndexManager<Index>, ChunkVector>>>;

  using RealReverse = RealReverseGen<double>;

  template<typename Real, typename Gradient = Real, typename Index = int>
  using RealReverseUncheckedGen = ActiveType<JacobianLinearTape<JacobianTapeTypes<Real, Gradient, LinearIndexManager<Index>, BlockVector>>>;

  using RealReverseUnchecked = RealReverseUncheckedGen<double>;

  template<typename Real, typename Gradient = Real, typename IndexManager = MultiUseIndexManager<int>>
  using RealReverseIndexGen = ActiveType<JacobianReuseTape<JacobianTapeTypes<Real, Gradient, IndexManager, ChunkVector>>>;

  using RealReverseIndex = RealReverseIndexGen<double>;

  template<typename Real, typename Gradient = Real, typename IndexManager = MultiUseIndexManager<int>>
  using RealReverseIndexUncheckedGen = ActiveType<JacobianReuseTape<JacobianTapeTypes<Real, Gradient, IndexManager, ChunkVector>>>;

  using RealReverseIndexUnchecked = RealReverseIndexUncheckedGen<double>;

  template<typename Real, typename Gradient = Real, typename Index = int, template <typename> class StatementEvaluator = InnerStatementEvaluator>
  using RealReversePrimalGen = ActiveType<PrimalValueLinearTape<PrimalValueTapeTypes<Real, Gradient, LinearIndexManager<Index>, StatementEvaluator, ChunkVector>>>;

  using RealReversePrimal = RealReversePrimalGen<double>;

  template<typename Real, typename Gradient = Real, typename Index = int, template <typename> class StatementEvaluator = InnerStatementEvaluator>
  using RealReversePrimalUncheckedGen = ActiveType<PrimalValueLinearTape<PrimalValueTapeTypes<Real, Gradient, LinearIndexManager<Index>, StatementEvaluator, ChunkVector>>>;

  using RealReversePrimalUnchecked = RealReversePrimalUncheckedGen<double>;

  template<typename Real, typename Gradient = Real, typename IndexManager = MultiUseIndexManager<int>, template <typename> class StatementEvaluator = InnerStatementEvaluator>
  using RealReversePrimalIndexGen = ActiveType<PrimalValueReuseTape<PrimalValueTapeTypes<Real, Gradient, IndexManager, StatementEvaluator, ChunkVector>>>;

  using RealReversePrimalIndex = RealReversePrimalIndexGen<double>;

  template<typename Real, typename Gradient = Real, typename IndexManager = MultiUseIndexManager<int>, template <typename> class StatementEvaluator = InnerStatementEvaluator>
  using RealReversePrimalIndexUncheckedGen = ActiveType<PrimalValueReuseTape<PrimalValueTapeTypes<Real, Gradient, IndexManager, StatementEvaluator, ChunkVector>>>;

  using RealReversePrimalIndexUnchecked = RealReversePrimalIndexUncheckedGen<double>;

}
