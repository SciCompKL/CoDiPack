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

/*
 * In order to include this file the user has to define the preprocessor macro ENUM. *
 *
 * The define ENUM will be undefined at the end of this template.
 */

#ifndef ENUM
  #error Please define the name of the ENUM.
#endif

// Create a correct include environment for viewing and programming in an IDE.
#ifndef ENUM
  #define PROXY

  #include "../misc/macros.hpp"
  #include "../config.h"
  #include "enumBitset.hpp"
  #define ENUM Enum

namespace codi {
#endif

  /// Return a EnumBitset structure when to enums are ored.
  CODI_INLINE EnumBitset<ENUM> operator|(ENUM a, ENUM b) {
    return EnumBitset<ENUM>(a) | b;
  }

  /// Return a boolean structure when to enums are anded.
  CODI_INLINE bool operator&(ENUM a, ENUM b) {
    return a == b;
  }

// Create a correct include environment for viewing and programming in an IDE.
#ifdef PROXY
  #undef PROXY
}
#endif

#undef ENUM
