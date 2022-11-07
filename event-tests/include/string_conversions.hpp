/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 *  - SciComp, TU Kaiserslautern:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "codi.hpp"

std::string to_string(codi::EventHints::Direction direction) {
  switch (direction) {
    case codi::EventHints::Direction::Forward:
      return "forward";
    case codi::EventHints::Direction::Reverse:
      return "reverse";
    default:
      return "unknown";
  }
}

std::string to_string(codi::EventHints::Endpoint endpoint) {
  switch (endpoint) {
    case codi::EventHints::Endpoint::Begin:
      return "begin";
    case codi::EventHints::Endpoint::End:
      return "end";
    default:
      return "unknown";
  }
}

std::string to_string(codi::EventHints::Statement statement) {
  switch (statement) {
    case codi::EventHints::Statement::Expression:
      return "expression";
    case codi::EventHints::Statement::Copy:
      return "copy";
    case codi::EventHints::Statement::Passive:
      return "passive";
    default:
      return "unknown";
  }
}

std::string to_string(codi::EventHints::Reset reset) {
  switch (reset) {
    case codi::EventHints::Reset::Full:
      return "full";
    case codi::EventHints::Reset::Hard:
      return "hard";
    case codi::EventHints::Reset::To:
      return "to";
    default:
      return "unknown";
  }
}
