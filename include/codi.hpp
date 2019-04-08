/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#pragma once

#include "codi/activeReal.hpp"
#include "codi/numericLimits.hpp"
#include "codi/referenceActiveReal.hpp"
#include "codi/tapeTypes.hpp"
#include "codi/tapes/forwardEvaluation.hpp"
#include "codi/tapes/jacobiTape.hpp"
#include "codi/tapes/jacobiIndexTape.hpp"
#include "codi/tapes/primalValueTape.hpp"
#include "codi/tapes/primalValueIndexTape.hpp"
#include "codi/tapes/indices/linearIndexHandler.hpp"
#include "codi/tapes/indices/reuseIndexHandler.hpp"
#include "codi/tapes/indices/reuseIndexHandlerUseCount.hpp"
#include "codi/tapes/handles/staticFunctionHandleFactory.hpp"
#include "codi/tapes/handles/staticObjectHandleFactory.hpp"
#include "codi/tools/dataStore.hpp"
#include "codi/tools/derivativeHelper.hpp"
#include "codi/tools/direction.hpp"
#include "codi/tools/externalFunctionHelper.hpp"
#include "codi/tools/preaccumulationHelper.hpp"
#include "codi/tools/statementPushHelper.hpp"
#include "codi/tools/tapeVectorHelper.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief The default forward type in CoDiPack with a generalized calculation type.
   *
   * See the documentation of #RealForward.
   *
   * @tparam     Real  The underlying calculation type for the AD evaluation. Needs to implement all mathematical functions.
   * @tparam Gradient  The type of the derivative values for the AD evaluation. Needs to implement an addition and multiplication operation.
   */
  template<typename Real, typename Gradient = Real>
  using RealForwardGen = ActiveReal<ForwardEvaluation<Real, Gradient> >;

  /**
   * @brief The default forward type in CoDiPack.
   *
   * This type is used to evaluate forward or tangent derivatives. For the function
   * \f[ y = f(x) \f]
   * the forward type calculates in addition to the primal computation
   * \f[ \dot{y} = \frac{df}{dx}(x)\cdot \dot {x}. \f]
   *
   * If you have the following program.
   * \code{.cpp}
   *  double a = 3.0;
   *  double b = a * a;
   * \endcode
   * the differentiated version looks like
   * \code{.cpp}
   *  RealForward a = 3.0;
   *  a.setGradient(1.0);
   *
   *  RealForward b = a * a;
   *  assert(b.getGradient() == 6.0);
   * \endcode
   */
  typedef RealForwardGen<double, double> RealForward;

  /**
   * @brief Vector mode of the #RealForward type.
   *
   * The type of the direction can be accessed with RealForwardVec<dim>::GradientValue<br>
   *
   * See @ref Tutorial6 for details.
   *
   * @tparam dim  The fixed dimension of the vector.
   */
  template<size_t dim>
  using RealForwardVec = RealForwardGen<double, Direction<double, dim> >;

  /**
   * @brief The default reverse type in CoDiPack with a generalized calculation type.
   *
   * See the documentation of #RealReverse.
   *
   * @tparam     Real  The underlying calculation type for the AD evaluation. Needs to implement all mathematical functions.
   * @tparam Gradient  The type of the derivative values for the AD evaluation. Needs to implement an addition and multiplication operation.
   */
  template<typename Real, typename Gradient = Real>
  using RealReverseGen = ActiveReal<JacobiTape<JacobiTapeTypes<ReverseTapeTypes<Real, Gradient, LinearIndexHandler<int> >, ChunkVector > > >;

  /**
   * @brief The default reverse type in CoDiPack.
   *
   * This type is used to evaluate reverse or adjoint derivatives. For the function
   * \f[ y = f(x) \f]
   * the reverse type stores information to calculate
   * \f[ \bar{x} = \frac{df}{dx}^T(x)\cdot \bar{y}. \f]
   *
   * If you have the following program.
   * \code{.cpp}
   *  double a = 3.0;
   *  double b = a * a;
   * \endcode
   * the differentiated version looks like
   * \code{.cpp}
   *  // For convenience get a reference to the global tape.
   *  RealReverse::TapeType& tape = RealReverse::globalTape;
   *  tape.setActive();
   *  RealReverse a = 3.0;
   *  tape.registerInput(a);
   *
   *  RealReverse b = a * a;
   *  tape.registerOutput(b);
   *  tape.setPassive();
   *
   *  b.setGradient(1.0);
   *  tape.evaluate();
   *  assert(a.getGradient() == 6.0);
   * \endcode
   */
  typedef RealReverseGen<double, double > RealReverse;

  /**
   * @brief Vector mode of the #RealReverse type.
   *
   * The type of the direction can be accessed with RealReverseVec<dim>::GradientValue<br>
   *
   * See @ref Tutorial6 for details.
   *
   * @tparam dim  The fixed dimension of the vector.
   */
  template<size_t dim>
  using RealReverseVec = RealReverseGen<double, Direction<double, dim> >;

  /**
   * @brief The reverse type in CoDiPack with a generalized calculation type and an unchecked tape.
   *
   * See the documentation of #RealReverseUnchecked.
   *
   * @tparam     Real  The underlying calculation type for the AD evaluation. Needs to implement all mathematical functions.
   * @tparam Gradient  The type of the derivative values for the AD evaluation. Needs to implement an addition and multiplication operation.
   */
  template<typename Real, typename Gradient = Real>
  using RealReverseUncheckedGen = ActiveReal<JacobiTape<JacobiTapeTypes<ReverseTapeTypes<Real, Gradient, LinearIndexHandler<int> >, SingleChunkVector > > >;

  /**
   * @brief The reverse type in CoDiPack with an unchecked tape.
   *
   * For details on the AD reverse mode see #RealReverse.
   *
   * This reverse type uses a tape which has no checks on the available space. The user has to specify at the beginning of
   * the calculation how many statements will be evaluated during the calculation and how many active arguments occur in
   * all statements.
   *
   * This tape should only be used if you are an experienced AD user or developer. If you are new to AD please use the
   * #RealReverse type. The tape of the #RealReverse type has bounds checking and increases the storage if more is needed.
   *
   * If you have the following program.
   * \code{.cpp}
   *  double a = 3.0;
   *  double b = a * a;
   * \endcode
   * the differentiated version looks like
   * \code{.cpp}
   *  // For convenience get a reference to the global tape.
   *  RealReverseUnchecked::TapeType& tape = RealReverseUnchecked::globalTape;
   *  // The only but very important change with respect to RealReverse. Without this you should get a segmentation fault.
   *  tape.resize(2, 2);
   *  tape.setActive();
   *  RealReverseUnchecked a = 3.0;
   *  tape.registerInput(a);
   *
   *  RealReverseUnchecked b = a * a;
   *  tape.registerOutput(b);
   *  tape.setPassive();
   *
   *  b.setGradient(1.0);
   *  tape.evaluate();
   *  assert(a.getGradient() == 6.0);
   * \endcode
   */
  typedef RealReverseUncheckedGen<double, double> RealReverseUnchecked;

  /**
   * @brief The reverse type in CoDiPack with a generalized calculation type and an index reuse tape.
   *
   * See the documentation of #RealReverseIndex.
   *
   * @tparam     Real  The underlying calculation type for the AD evaluation. Needs to implement all mathematical functions.
   * @tparam Gradient  The type of the derivative values for the AD evaluation. Needs to implement an addition and multiplication operation.
   */
  template<typename Real, typename Gradient = Real>
  using RealReverseIndexGen = ActiveReal<JacobiIndexTape<JacobiIndexTapeTypes<ReverseTapeTypes<Real, Gradient, ReuseIndexHandlerUseCount<int> >, ChunkVector> > >;

  /**
   * @brief A reverse type like the default reverse type in CoDiPack but with index reuse.
   *
   * The difference between this type and the #RealReverse is the handling of the indices.
   * This type stores deleted indices and uses them again on other variables. Usually this
   * tape behaves as the #RealReverse but it is not compatible with c-like memory operations
   * like memset and memcpy.
   *
   */
  typedef RealReverseIndexGen<double, double> RealReverseIndex;

  /**
   * @brief Vector mode of the #RealReverseIndex type.
   *
   * The type of the direction can be accessed with RealReverseIndexVec<dim>::GradientValue<br>
   *
   * See @ref Tutorial6 for details.
   *
   * @tparam dim  The fixed dimension of the vector.
   */
  template<size_t dim>
  using RealReverseIndexVec = RealReverseIndexGen<double, Direction<double, dim> >;

  /**
   * @brief The reverse type in CoDiPack  with a generalized calculation type and an unchecked index reuse tape.
   *
   * See the documentation of #RealReverseIndexUnchecked.
   *
   * @tparam     Real  The underlying calculation type for the AD evaluation. Needs to implement all mathematical functions.
   * @tparam Gradient  The type of the derivative values for the AD evaluation. Needs to implement an addition and multiplication operation.
   */
  template<typename Real, typename Gradient = Real>
  using RealReverseIndexUncheckedGen = ActiveReal<JacobiIndexTape<JacobiIndexTapeTypes<ReverseTapeTypes<Real, Gradient, ReuseIndexHandlerUseCount<int> >, SingleChunkVector> > >;

  /**
   * @brief A reverse type like the unchecked reverse type in CoDiPack but with index reuse.
   *
   * See the documentation of #RealReverseIndex and #RealReverseUnchecked.
   */
  typedef RealReverseIndexUncheckedGen<double, double> RealReverseIndexUnchecked;

  /**
   * @brief The primal value reverse type in CoDiPack with a generalized calculation type.
   *
   * See the documentation of #RealReversePrimal.
   *
   * @tparam     Real  The underlying calculation type for the AD evaluation. Needs to implement all mathematical functions.
   * @tparam Gradient  The type of the derivative values for the AD evaluation. Needs to implement an addition and multiplication operation.
   */
  template<typename Real, typename Gradient = Real>
  using RealReversePrimalGen = ActiveReal<PrimalValueTape<PrimalValueTapeTypes<ReverseTapeTypes<Real, Gradient, LinearIndexHandler<int> >, StaticFunctionHandleFactory, ChunkVector> > >;

  /**
   * @brief A reverse type like the default reverse type in CoDiPack but with primal value taping instead of Jacobian taping.
   *
   * The difference between this type and the #RealReverse is how the derivatives are computed and stored
   * the Jacobian approach computes the derivatives directly and stores them. The primal value taping
   * stores the primal values of the expression. With these primal values it computes the derivatives in the
   * reverse mode.
   */
  typedef RealReversePrimalGen<double, double> RealReversePrimal;

  /**
   * @brief Vector mode of the #RealReversePrimal type.
   *
   * The type of the direction can be accessed with RealReversePrimalVec<dim>::GradientValue<br>
   *
   * See @ref Tutorial6 for details.
   *
   * @tparam dim  The fixed dimension of the vector.
   */
  template<size_t dim>
  using RealReversePrimalVec = RealReversePrimalGen<double, Direction<double, dim> >;

  /**
   * @brief The primal value reverse type in CoDiPack with a generalized calculation type and an unchecked version.
   *
   * See the documentation of #RealReversePrimal.
   *
   * @tparam     Real  The underlying calculation type for the AD evaluation. Needs to implement all mathematical functions.
   * @tparam Gradient  The type of the derivative values for the AD evaluation. Needs to implement an addition and multiplication operation.
   */
  template<typename Real, typename Gradient = Real>
  using RealReversePrimalUncheckedGen = ActiveReal<PrimalValueTape<PrimalValueTapeTypes<ReverseTapeTypes<Real, Gradient, LinearIndexHandler<int> >, StaticFunctionHandleFactory, SingleChunkVector> > >;

  /**
   * @brief The primal value reverse type in CoDiPack with an unchecked tape.
   *
   * See the documentation of #RealReversePrimal and #RealReverseUnchecked.
   */
  typedef RealReversePrimalUncheckedGen<double, double> RealReversePrimalUnchecked;

  /**
   * @brief The primal value reverse type in CoDiPack with an index management and with a generalized calculation type.
   *
   * See the documentation of #RealReversePrimal and #RealReverseIndex.
   *
   * @tparam     Real  The underlying calculation type for the AD evaluation. Needs to implement all mathematical functions.
   * @tparam Gradient  The type of the derivative values for the AD evaluation. Needs to implement an addition and multiplication operation.
   */
  template<typename Real, typename Gradient = Real>
  using RealReversePrimalIndexGen = ActiveReal<PrimalValueIndexTape<IndexPrimalValueTapeTypes<ReverseTapeTypes<Real, Gradient, ReuseIndexHandlerUseCount<int> >, StaticObjectHandleFactory, ChunkVector > > >;


  /**
   * @brief The primal value reverse type in CoDiPack with an index management.
   *
   * See the documentation of #RealReversePrimal and #RealReverseIndex.
   */
  typedef RealReversePrimalIndexGen<double, double> RealReversePrimalIndex;

  /**
   * @brief Vector mode of the #RealReversePrimalIndex type.
   *
   * The type of the direction can be accessed with RealReversePrimalIndexVec<dim>::GradientValue<br>
   *
   * See @ref Tutorial6 for details.
   *
   * @tparam dim  The fixed dimension of the vector.
   */
  template<size_t dim>
  using RealReversePrimalIndexVec = RealReversePrimalIndexGen<double, Direction<double, dim> >;

  /**
   * @brief The primal value reverse type in CoDiPack with an index management in an unchecked version and with a generalized calculation type.
   *
   * See the documentation of #RealReversePrimal, #RealReverseIndex and #RealReverseUnchecked
   *
   * @tparam     Real  The underlying calculation type for the AD evaluation. Needs to implement all mathematical functions.
   * @tparam Gradient  The type of the derivative values for the AD evaluation. Needs to implement an addition and multiplication operation.
   */
  template<typename Real, typename Gradient = Real>
  using RealReversePrimalIndexUncheckedGen = ActiveReal<PrimalValueIndexTape<IndexPrimalValueTapeTypes<ReverseTapeTypes<Real, Gradient, ReuseIndexHandlerUseCount<int> >, StaticObjectHandleFactory, SingleChunkVector > > >;

  /**
   * @brief The primal value reverse type in CoDiPack with an index management in an unchecked version.
   *
   * See the documentation of #RealReversePrimal, #RealReverseIndex and #RealReverseUnchecked
   */
  typedef RealReversePrimalIndexUncheckedGen<double, double> RealReversePrimalIndexUnchecked;
}
