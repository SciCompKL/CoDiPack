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

#include <type_traits>

#include "forwardEvaluation.hpp"
#include "jacobiIndexTape.hpp"
#include "jacobiTape.hpp"
#include "primalValueIndexTape.hpp"
#include "primalValueTape.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  template<typename Tape>
  using isForwardTape = std::is_base_of<ForwardEvaluation<typename Tape::TapeTypes>, Tape>;

  template<typename Tape>
  using enableIfForwardTape = typename std::enable_if<isForwardTape<Tape>::value>::type;

  template<typename Tape>
  using isPrimalValueTape = std::integral_constant<bool,
         std::is_base_of<PrimalValueIndexTape<typename Tape::TapeTypes>, Tape>::value
      || std::is_base_of<PrimalValueTape<typename Tape::TapeTypes>, Tape>::value
    >;

  template<typename Tape>
  using enableIfPrimalValueTape = typename std::enable_if<isPrimalValueTape<Tape>::value>::type;

  template<typename Tape>
  using isJacobianTape = std::integral_constant<bool,
       std::is_base_of<JacobiIndexTape<typename Tape::TapeTypes>, Tape>::value
    || std::is_base_of<JacobiTape<typename Tape::TapeTypes>, Tape>::value
    >;

  template<typename Tape>
  using enableIfJacobianTape = typename std::enable_if<isJacobianTape<Tape>::value>::type;

  template<typename Tape>
  using isReverseTape = std::integral_constant<bool,
         isJacobianTape<Tape>::value
      || isPrimalValueTape<Tape>::value
    >;

  template<typename Tape>
  using enableIfReverseTape = typename std::enable_if<isReverseTape<Tape>::value>::type;
}
