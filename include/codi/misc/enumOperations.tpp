/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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

/*
 * This file creates:
 * - EnumBitset<ENUM> operator|(ENUM, ENUM)
 * - EnumBitset<ENUM> operator&(ENUM, ENUM)
 *
 * In order to include this file the user has to define the preprocessor macro ENUM.
 *
 * The defined ENUM will be undefined at the end of this template.
 */

#ifndef ENUM
  #error Please define the name of the ENUM.
#endif

// Create a correct include environment for viewing and programming in an IDE.
#ifndef ENUM
  #define PROXY

  #include "../config.h"
  #include "enumBitset.hpp"
  #include "enumInterface.hpp"
  #include "macros.hpp"
  #define ENUM EnumInterface

namespace codi {
#endif

  /// Return an EnumBitset structure when two enums are combined with or.
  CODI_INLINE EnumBitset<ENUM> operator|(ENUM a, ENUM b) {
    return EnumBitset<ENUM>(a) | b;
  }

  /// Return a boolean structure when two enums are combined with and.
  CODI_INLINE EnumBitset<ENUM> operator&(ENUM a, ENUM b) {
    return EnumBitset<ENUM>(a) & b;
  }

// Create a correct include environment for viewing and programming in an IDE.
#ifdef PROXY
  #undef PROXY
}
#endif

#undef ENUM
