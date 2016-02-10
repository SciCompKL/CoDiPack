/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
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

#include "activeReal.hpp"
#include "tapes/forwardEvaluation.hpp"
#include "tapes/jacobiTape.hpp"
#include "tapes/jacobiIndexTape.hpp"
#include "tapes/indices/linearIndexHandler.hpp"
#include "tapes/indices/reuseIndexHandler.hpp"
#include "tools/dataStore.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {
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
  typedef ActiveReal<double, ForwardEvaluation<double> > RealForward;

  /**
   * @brief The default forward type in CoDiPack with float as the real value type.
   *
   * See the documentation of #RealForward.
   */
  typedef ActiveReal<float, ForwardEvaluation<float> > RealForwardFloat;

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
  typedef ActiveReal<double, JacobiTape<ChunkTapeTypes<double, LinearIndexHandler<int> > > > RealReverse;

  /**
   * @brief The default reverse type in CoDiPack with float as the real value type.
   *
   * See the documentation of #RealReverse.
   */
  typedef ActiveReal<float, JacobiTape<ChunkTapeTypes<float, LinearIndexHandler<int> > > > RealReverseFloat;

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
  typedef ActiveReal<double, JacobiTape<SimpleTapeTypes<double, LinearIndexHandler<int> > > > RealReverseUnchecked;

  /**
   * @brief The reverse type in CoDiPack with float as the real value type and an unchecked tape.
   *
   * See the documentation of #RealReverseUnchecked.
   */
  typedef ActiveReal<float, JacobiTape<SimpleTapeTypes<float, LinearIndexHandler<int> > > > RealReverseUncheckedFloat;


  /**
   * @brief A reverse type like the default reverse type in CoDiPack but with index reuse.
   *
   * The difference between this type and the #RealReverse is the handling of the indices.
   * This type stores deleted indices and uses them again on other variables. Usually this
   * tape behaves as the #RealReverse but it is not compatible with c-like memory operations
   * like memset and memcpy.
   *
   */
  typedef ActiveReal<double, JacobiIndexTape<ChunkIndexTapeTypes<double, ReuseIndexHandler<int> > > > RealReverseIndex;

  /**
   * @brief The reverse type in CoDiPack with float as the real value type and an index reuse tape.
   *
   * See the documentation of #RealReverseIndex.
   */
  typedef ActiveReal<float, JacobiIndexTape<ChunkIndexTapeTypes<float, ReuseIndexHandler<int> > > > RealReverseIndexFloat;

  /**
   * @brief A reverse type like the unchecked reverse type in CoDiPack but with index reuse.
   *
   * See the documentation of #RealReverseIndex and #RealReverseUnchecked.
   */
  typedef ActiveReal<double, JacobiIndexTape<SimpleIndexTapeTypes<double, ReuseIndexHandler<int> > > > RealReverseIndexUnchecked;

  /**
   * @brief The reverse type in CoDiPack with float as the real value type and an unchecked index reuse tape.
   *
   * See the documentation of #RealReverseIndexUnchecked.
   */
  typedef ActiveReal<float, JacobiIndexTape<SimpleIndexTapeTypes<float, ReuseIndexHandler<int> > > > RealReverseIndexUncheckedFloat;


}
