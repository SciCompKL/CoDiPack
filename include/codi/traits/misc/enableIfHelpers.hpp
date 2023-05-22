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

#include <type_traits>

#include "../../config.h"
#include "../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Enable if abbreviation for std::is_base_of
  template<typename Base, typename Impl, typename R = void>
  using enable_if_base_of = std::enable_if<std::is_base_of<Base, Impl>::value, R>;

  /// Enable if abbreviation for \c "std::is_same"
  template<typename T1, typename T2, typename R = void>
  using enable_if_same = std::enable_if<std::is_same<T1, T2>::value, R>;

  /// Enable if abbreviation for \c "!std::is_same"
  template<typename T1, typename T2, typename R = void>
  using enable_if_not_same = std::enable_if<!std::is_same<T1, T2>::value, R>;
}
