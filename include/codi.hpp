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
#include "codi/tools/data/jacobian.hpp"
#include "codi/tools/helpers/customAdjointVectorHelper.hpp"
#include "codi/tools/helpers/externalFunctionHelper.hpp"
#include "codi/tools/helpers/preaccumulationHelper.hpp"
#include "codi/tools/helpers/statementPushHelper.hpp"
#include "codi/tools/helpers/tapeHelper.hpp"
#include "codi/tools/higherOrderAccess.hpp"
#include "codi/traits/numericLimits.hpp"
#include "codi/traits/tapeTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * \section sec_forwardAD Forward AD Equation
   *
   * \f$ \dot w = d\phi/du * \dot u \f$
   *
   * \section sec_reverseAD Reverse AD Equation
   *
   * \f$ \bar u = d\phi/du^T * \bar w \f$
   *
   * \section sec_namingConventions Mathematical naming conventions
   *
   * The function is defined by:
   *  \f[ y = f(x) \f]
   * where \f$ x \in \R^n \f$ is the input vector with the size \f$ n \f$ and \f$y \in \R^m \f$ is the output vector
   * with the size \f$ m \f$
   *
   * The Jacobian of \f$ f \f$ is defined by:
   *  \f[ J = \frac{\d f}{\d x} \in \R^{m \times n} \f]
   * The number of rows (\f$ m \f$) represents the number of output variables and the number of columns (\f$ n \f$)
   * represents the number of input variable. The derivative for the i-th output with respect to the j-th input is
   * represented by \f$ J_{i,j} \f$.
   *
   * The Hessian of \f$ f \f$ is defined by:
   *  \f[ H = \frac{\d^2 f}{\d^2 x} \in \R^{m \times n \times n} \f]
   * The first dimension (\f$ m \f$) represents the number of output variables, the second and third dimension (\f$ n \f$) represents the
   * first and second derivative with respect to the input variables.
   * The second derivative for the i-th output with respect to the j-th and k-th input is
   * represented by \f$ H_{i,j,k} \f$.
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
