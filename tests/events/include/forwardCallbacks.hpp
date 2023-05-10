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

#include <codi.hpp>

#include "string_conversions.hpp"

namespace ForwardCallbacks {

  template<typename Tape>
  void onStatementPrimal(Tape&, typename Tape::Real const& lhsValue, typename Tape::Identifier const& lhsDotValue,
                         typename Tape::Real const& newValue, codi::EventHints::Statement statement, void*) {
    std::cout << "StatementPrimal " << to_string(statement) << " lhsValue " << lhsValue << " lhsDotValue "
              << lhsDotValue << " newValue " << newValue << std::endl;
  }

  template<typename Tape>
  std::list<typename codi::EventSystem<Tape>::Handle> registerAll() {
    std::list<typename codi::EventSystem<Tape>::Handle> handles;
    handles.push_back(codi::EventSystem<Tape>::registerStatementPrimalListener(onStatementPrimal<Tape>));
    return handles;
  }
}

template<typename Tape>
void deregisterCallbacks(std::list<typename codi::EventSystem<Tape>::Handle> const& handles) {
  for (auto const& handle : handles) {
    codi::EventSystem<Tape>::deregisterListener(handle);
  }
}
