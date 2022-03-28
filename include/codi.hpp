/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * CoDiPack is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU
 * General Public License along with CoDiPack.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "codi/config.h"
#include "codi/expressions/activeType.hpp"
#include "codi/expressions/activeTypeWrapper.hpp"
#include "codi/expressions/immutableActiveType.hpp"
#include "codi/expressions/real/allOperators.hpp"
#include "codi/expressions/referenceActiveType.hpp"
#include "codi/misc/enumBitset.hpp"
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
#include "codi/tools/data/aggregatedTypeVectorAccessWrapper.hpp"
#include "codi/tools/data/direction.hpp"
#include "codi/tools/data/externalFunctionUserData.hpp"
#include "codi/tools/data/jacobian.hpp"
#include "codi/tools/derivativeAccess.hpp"
#include "codi/tools/helpers/customAdjointVectorHelper.hpp"
#include "codi/tools/helpers/externalFunctionHelper.hpp"
// #include "codi/tools/helpers/evaluationHelper.hpp" // Included at the end of this file.
#include "codi/tools/helpers/linearSystem/linearSystemHandler.hpp"
#include "codi/tools/helpers/preaccumulationHelper.hpp"
#include "codi/tools/helpers/statementPushHelper.hpp"
#include "codi/tools/helpers/tapeHelper.hpp"
#include "codi/traits/computationTraits.hpp"
#include "codi/traits/numericLimits.hpp"
#include "codi/traits/tapeTraits.hpp"

#if CODI_EnableMPI
  #include "codi/tools/mpi/codiMpiTypes.hpp"
#endif

#if CODI_EnableEigen
  #include "codi/tools/helpers/linearSystem/eigenLinearSystem.hpp"
#endif

#if CODI_EnableEnzyme
  #include "codi/tools/helpers/enzymeExternalFunctionHelper.hpp"
#endif

/** \copydoc codi::Namespace */
namespace codi {

#define CODI_MAJOR_VERSION 2
#define CODI_MINOR_VERSION 1
#define CODI_BUILD_VERSION 0
#define CODI_VERSION "2.1.0"

  /// General forward AD type. See \ref sec_forwardAD for a forward mode AD explanation or \ref ActiveTypeList for a
  /// list of all types.
  template<typename Real, typename Gradient = Real>
  using RealForwardGen = ActiveType<ForwardEvaluation<Real, Gradient>>;

  /// Default forward AD type. See \ref sec_forwardAD for a forward mode AD explanation or \ref ActiveTypeList for a
  /// list of all types.
  using RealForward = RealForwardGen<double, double>;

  /// General vector forward AD type. See \ref sec_forwardAD for a forward mode AD explanation or \ref ActiveTypeList
  /// for a list of all types.
  template<size_t dim>
  using RealForwardVec = RealForwardGen<double, Direction<double, dim>>;

  /// General reverse AD type. See \ref sec_reverseAD for a reverse mode AD explanation or \ref ActiveTypeList for a
  /// list of all types.
  ///
  /// Jacobian taping approach with linear index handling.
  template<typename Real, typename Gradient = Real, typename Index = int>
  using RealReverseGen =
      ActiveType<JacobianLinearTape<JacobianTapeTypes<Real, Gradient, LinearIndexManager<Index>, DefaultChunkedData>>>;

  /// \copydoc codi::RealReverseGen
  using RealReverse = RealReverseGen<double>;

  /// \copydoc codi::RealReverseGen
  template<size_t dim>
  using RealReverseVec = RealReverseGen<double, Direction<double, dim>>;

  /// General unchecked reverse AD type. See \ref sec_reverseAD for a reverse mode AD explanation or \ref ActiveTypeList
  /// for a list of all types.
  ///
  /// Requires preallocation of data. See DataManagementTapeInterface.
  ///
  /// Jacobian taping approach with linear index handling.
  template<typename Real, typename Gradient = Real, typename Index = int>
  using RealReverseUncheckedGen =
      ActiveType<JacobianLinearTape<JacobianTapeTypes<Real, Gradient, LinearIndexManager<Index>, DefaultBlockData>>>;

  /// \copydoc codi::RealReverseUncheckedGen
  using RealReverseUnchecked = RealReverseUncheckedGen<double>;

  /// General reverse AD type. See \ref sec_reverseAD for a reverse mode AD explanation or \ref ActiveTypeList for a
  /// list of all types.
  ///
  /// Jacobian taping approach with reuse index handling.
  template<typename Real, typename Gradient = Real, typename IndexManager = MultiUseIndexManager<int>>
  using RealReverseIndexGen =
      ActiveType<JacobianReuseTape<JacobianTapeTypes<Real, Gradient, IndexManager, DefaultChunkedData>>>;

  /// \copydoc codi::RealReverseIndexGen
  using RealReverseIndex = RealReverseIndexGen<double>;

  /// \copydoc codi::RealReverseIndexGen
  template<size_t dim>
  using RealReverseIndexVec = RealReverseIndexGen<double, Direction<double, dim>>;

  /// General unchecked reverse AD type. See \ref sec_reverseAD for a reverse mode AD explanation or \ref ActiveTypeList
  /// for a list of all types.
  ///
  /// Requires preallocation of data. See DataManagementTapeInterface.
  ///
  /// Jacobian taping approach with reuse index handling.
  template<typename Real, typename Gradient = Real, typename IndexManager = MultiUseIndexManager<int>>
  using RealReverseIndexUncheckedGen =
      ActiveType<JacobianReuseTape<JacobianTapeTypes<Real, Gradient, IndexManager, DefaultChunkedData>>>;

  /// \copydoc codi::RealReverseIndexUncheckedGen
  using RealReverseIndexUnchecked = RealReverseIndexUncheckedGen<double>;

  /// General reverse AD type. See \ref sec_reverseAD for a reverse mode AD explanation or \ref ActiveTypeList for a
  /// list of all types.
  ///
  /// Primal value taping approach with linear index handling.
  template<typename Real, typename Gradient = Real, typename Index = int,
           template<typename> class StatementEvaluator = InnerStatementEvaluator>
  using RealReversePrimalGen = ActiveType<PrimalValueLinearTape<
      PrimalValueTapeTypes<Real, Gradient, LinearIndexManager<Index>, StatementEvaluator, DefaultChunkedData>>>;

  /// \copydoc codi::RealReversePrimalGen
  using RealReversePrimal = RealReversePrimalGen<double>;

  /// \copydoc codi::RealReversePrimalGen
  template<size_t dim>
  using RealReversePrimalVec = RealReversePrimalGen<double, Direction<double, dim>>;

  /// General unchecked reverse AD type. See \ref sec_reverseAD for a reverse mode AD explanation or \ref ActiveTypeList
  /// for a list of all types.
  ///
  /// Requires preallocation of data. See DataManagementTapeInterface.
  ///
  /// Primal value taping approach with linear index handling.
  template<typename Real, typename Gradient = Real, typename Index = int,
           template<typename> class StatementEvaluator = InnerStatementEvaluator>
  using RealReversePrimalUncheckedGen = ActiveType<PrimalValueLinearTape<
      PrimalValueTapeTypes<Real, Gradient, LinearIndexManager<Index>, StatementEvaluator, DefaultChunkedData>>>;

  /// \copydoc codi::RealReversePrimalUncheckedGen
  using RealReversePrimalUnchecked = RealReversePrimalUncheckedGen<double>;

  /// General reverse AD type. See \ref sec_reverseAD for a reverse mode AD explanation or \ref ActiveTypeList for a
  /// list of all types.
  ///
  /// Primal value taping approach with reuse index handling.
  template<typename Real, typename Gradient = Real, typename IndexManager = MultiUseIndexManager<int>,
           template<typename> class StatementEvaluator = InnerStatementEvaluator>
  using RealReversePrimalIndexGen = ActiveType<
      PrimalValueReuseTape<PrimalValueTapeTypes<Real, Gradient, IndexManager, StatementEvaluator, DefaultChunkedData>>>;

  /// \copydoc codi::RealReversePrimalIndexGen
  using RealReversePrimalIndex = RealReversePrimalIndexGen<double>;

  /// \copydoc codi::RealReversePrimalIndexGen
  template<size_t dim>
  using RealReversePrimalIndexVec = RealReversePrimalIndexGen<double, Direction<double, dim>>;

  /// General unchecked reverse AD type. See \ref sec_reverseAD for a reverse mode AD explanation or \ref ActiveTypeList
  /// for a list of all types.
  ///
  /// Requires preallocation of data. See DataManagementTapeInterface.
  ///
  /// Primal value taping approach with reuse index handling.
  template<typename Real, typename Gradient = Real, typename IndexManager = MultiUseIndexManager<int>,
           template<typename> class StatementEvaluator = InnerStatementEvaluator>
  using RealReversePrimalIndexUncheckedGen = ActiveType<
      PrimalValueReuseTape<PrimalValueTapeTypes<Real, Gradient, IndexManager, StatementEvaluator, DefaultChunkedData>>>;

  /// \copydoc codi::RealReversePrimalIndexUncheckedGen
  using RealReversePrimalIndexUnchecked = RealReversePrimalIndexUncheckedGen<double>;

  /**
   * @brief A regular CoDiPack type that can be used for Hessian computations in the TapeHelper.
   */
  using HessianComputationType = RealReversePrimalIndexGen<RealForwardVec<4>, Direction<RealForwardVec<4>, 4>>;

  /**
   * @brief A regular CoDiPack type that can be used for Hessian computations in the TapeHelper.
   *
   * This is the scalar version which does not use a vector mode.
   */
  using HessianComputationScalarType = RealReversePrimalIndexGen<RealForward>;

  /**
   * @brief A regular CoDiPack type that can be used for Jacobian computations in the TapeHelper.
   */
  using JacobianComputationType = RealReverseIndexVec<4>;

  /**
   * @brief A regular CoDiPack type that can be used for Jacobian computations in the TapeHelper.
   *
   * This is the scalar version which does not use a vector mode.
   */
  using JacobianComputationScalarType = RealReverseIndex;

}

#include "codi/tools/helpers/evaluationHelper.hpp"

#if CODI_EnableOpenMP
  #include "codi/tools/parallel/openmp/codiOpenMP.hpp"
#endif

#if CODI_EnableOpDiLib
  #include "codi/tools/parallel/openmp/codiOpDiLibTool.hpp"
#endif
