/**
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing, TU Kaiserslautern
 *
 * This file is part of CoDiPack.
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
 * Authors: TODO
 */
#pragma once

#include "activeReal.hpp"
#include "tapes/forwardEvaluation.hpp"
#include "tapes/simpleTape.hpp"
#include "tapes/chunkTape.hpp"
#include "tools/DataStore.hpp"

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
   * the differentitated version looks like
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
   * the differentitated version looks like
   * \code{.cpp}
   *  // For convinience get a reference to the global tape.
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
  typedef ActiveReal<double, ChunkTape<double, int> > RealReverse;

  /**
   * @brief The default reverse type in CoDiPack with float as the real value type.
   *
   * See the documentation of #RealReverse.
   */
  typedef ActiveReal<float, ChunkTape<float, int> > RealReverseFloat;

  /**
   * @brief The reverse type in CoDiPack with an unchecked tape.
   *
   * For details on the AD reverse mode see #RealReverse.
   *
   * This reverse type uses a tape which has no checks on the available space. The user has to specify at the beginning of
   * the calculation how many statements will be evaluted during the calculation and how many active arguments occour in
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
   * the differentitated version looks like
   * \code{.cpp}
   *  // For convinience get a reference to the global tape.
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
  typedef ActiveReal<double, SimpleTape<double, int> > RealReverseUnchecked;

  /**
   * @brief The reverse type in CoDiPack with float as the real value type and an unchecked tape.
   *
   * See the documentation of #RealReverseUnchecked.
   */
  typedef ActiveReal<float, SimpleTape<float, int> > RealReverseUncheckedFloat;


}
