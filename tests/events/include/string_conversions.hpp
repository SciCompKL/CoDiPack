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

#include "codi.hpp"

std::string to_string(codi::EventHints::EvaluationKind kind) {
  switch (kind) {
    case codi::EventHints::EvaluationKind::Primal:
      return "primal";
    case codi::EventHints::EvaluationKind::Forward:
      return "forward";
    case codi::EventHints::EvaluationKind::Reverse:
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
