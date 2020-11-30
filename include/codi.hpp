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
#include "codi/tools/data/direction.hpp"
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

  /// BlockData DataInterface used in all unchecked tapes.
  template<typename Chunk, typename NestedData = EmptyData>
  using DefaultBlockData = BlockData<Chunk, NestedData>;

  /// ChunkData DataInterface used in all regular tapes.
  template<typename Chunk, typename NestedData = EmptyData>
  using DefaultChunkedData = ChunkedData<Chunk, NestedData>;

  /// General forward AD type. See \ref sec_forwardAD for a forward mode AD explanation or \ref ActiveTypeList for a list of all types.
  template<typename Real, typename Gradient = Real>
  using RealForwardGen = ActiveType<ForwardEvaluation<Real, Gradient>>;

  /// Default forward AD type. See \ref sec_forwardAD for a forward mode AD explanation or \ref ActiveTypeList for a list of all types.
  using RealForward = RealForwardGen<double, double>;

  /// General vector forward AD type. See \ref sec_forwardAD for a forward mode AD explanation or \ref ActiveTypeList for a list of all types.
  template<size_t dim>
  using RealForwardVec = RealForwardGen<double, Direction<double, dim>>;

  /// General reverse AD type. See \ref sec_reverseAD for a reverse mode AD explanation or \ref ActiveTypeList for a list of all types.
  ///
  /// Jacobian taping approach with linear index handling.
  template<typename Real, typename Gradient = Real, typename Index = int>
  using RealReverseGen = ActiveType<JacobianLinearTape<JacobianTapeTypes<Real, Gradient, LinearIndexManager<Index>, DefaultChunkedData>>>;

  /// \copydoc codi::RealReverseGen
  using RealReverse = RealReverseGen<double>;

  /// \copydoc codi::RealReverseGen
  template<size_t dim>
  using RealReverseVec = RealReverseGen<double, Direction<double, dim>>;

  /// General unchecked reverse AD type. See \ref sec_reverseAD for a reverse mode AD explanation or \ref ActiveTypeList for a list of all types.
  ///
  /// Requires preallocation of data. See DataManagementTapeInterface.
  ///
  /// Jacobian taping approach with linear index handling.
  template<typename Real, typename Gradient = Real, typename Index = int>
  using RealReverseUncheckedGen = ActiveType<JacobianLinearTape<JacobianTapeTypes<Real, Gradient, LinearIndexManager<Index>, DefaultBlockData>>>;

  /// \copydoc codi::RealReverseUncheckedGen
  using RealReverseUnchecked = RealReverseUncheckedGen<double>;

  /// General reverse AD type. See \ref sec_reverseAD for a reverse mode AD explanation or \ref ActiveTypeList for a list of all types.
  ///
  /// Jacobian taping approach with reuse index handling.
  template<typename Real, typename Gradient = Real, typename IndexManager = MultiUseIndexManager<int>>
  using RealReverseIndexGen = ActiveType<JacobianReuseTape<JacobianTapeTypes<Real, Gradient, IndexManager, DefaultChunkedData>>>;

  /// \copydoc codi::RealReverseIndexGen
  using RealReverseIndex = RealReverseIndexGen<double>;

  /// \copydoc codi::RealReverseIndexGen
  template<size_t dim>
  using RealReverseIndexVec = RealReverseIndexGen<double, Direction<double, dim>>;

  /// General unchecked reverse AD type. See \ref sec_reverseAD for a reverse mode AD explanation or \ref ActiveTypeList for a list of all types.
  ///
  /// Requires preallocation of data. See DataManagementTapeInterface.
  ///
  /// Jacobian taping approach with reuse index handling.
  template<typename Real, typename Gradient = Real, typename IndexManager = MultiUseIndexManager<int>>
  using RealReverseIndexUncheckedGen = ActiveType<JacobianReuseTape<JacobianTapeTypes<Real, Gradient, IndexManager, DefaultChunkedData>>>;

  /// \copydoc codi::RealReverseIndexUncheckedGen
  using RealReverseIndexUnchecked = RealReverseIndexUncheckedGen<double>;

  /// General reverse AD type. See \ref sec_reverseAD for a reverse mode AD explanation or \ref ActiveTypeList for a list of all types.
  ///
  /// Primal value taping approach with linear index handling.
  template<typename Real, typename Gradient = Real, typename Index = int, template <typename> class StatementEvaluator = InnerStatementEvaluator>
  using RealReversePrimalGen = ActiveType<PrimalValueLinearTape<PrimalValueTapeTypes<Real, Gradient, LinearIndexManager<Index>, StatementEvaluator, DefaultChunkedData>>>;

  /// \copydoc codi::RealReversePrimalGen
  using RealReversePrimal = RealReversePrimalGen<double>;

  /// \copydoc codi::RealReversePrimalGen
  template<size_t dim>
  using RealReversePrimalVec = RealReversePrimalGen<double, Direction<double, dim>>;

  /// General unchecked reverse AD type. See \ref sec_reverseAD for a reverse mode AD explanation or \ref ActiveTypeList for a list of all types.
  ///
  /// Requires preallocation of data. See DataManagementTapeInterface.
  ///
  /// Primal value taping approach with linear index handling.
  template<typename Real, typename Gradient = Real, typename Index = int, template <typename> class StatementEvaluator = InnerStatementEvaluator>
  using RealReversePrimalUncheckedGen = ActiveType<PrimalValueLinearTape<PrimalValueTapeTypes<Real, Gradient, LinearIndexManager<Index>, StatementEvaluator, DefaultChunkedData>>>;

  /// \copydoc codi::RealReversePrimalUncheckedGen
  using RealReversePrimalUnchecked = RealReversePrimalUncheckedGen<double>;

  /// General reverse AD type. See \ref sec_reverseAD for a reverse mode AD explanation or \ref ActiveTypeList for a list of all types.
  ///
  /// Primal value taping approach with reuse index handling.
  template<typename Real, typename Gradient = Real, typename IndexManager = MultiUseIndexManager<int>, template <typename> class StatementEvaluator = InnerStatementEvaluator>
  using RealReversePrimalIndexGen = ActiveType<PrimalValueReuseTape<PrimalValueTapeTypes<Real, Gradient, IndexManager, StatementEvaluator, DefaultChunkedData>>>;

  /// \copydoc codi::RealReversePrimalIndexGen
  using RealReversePrimalIndex = RealReversePrimalIndexGen<double>;

  /// \copydoc codi::RealReversePrimalIndexGen
  template<size_t dim>
  using RealReversePrimalIndexVec = RealReversePrimalIndexGen<double, Direction<double, dim>>;

  /// General unchecked reverse AD type. See \ref sec_reverseAD for a reverse mode AD explanation or \ref ActiveTypeList for a list of all types.
  ///
  /// Requires preallocation of data. See DataManagementTapeInterface.
  ///
  /// Primal value taping approach with reuse index handling.
  template<typename Real, typename Gradient = Real, typename IndexManager = MultiUseIndexManager<int>, template <typename> class StatementEvaluator = InnerStatementEvaluator>
  using RealReversePrimalIndexUncheckedGen = ActiveType<PrimalValueReuseTape<PrimalValueTapeTypes<Real, Gradient, IndexManager, StatementEvaluator, DefaultChunkedData>>>;

  /// \copydoc codi::RealReversePrimalIndexUncheckedGen
  using RealReversePrimalIndexUnchecked = RealReversePrimalIndexUncheckedGen<double>;

}
