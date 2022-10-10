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

#include "../../misc/macros.hpp"
#include "../../config.h"
#include "../data/position.hpp"
#include "fullTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Edit tapes after they have been recorded.
   *
   * These interface functions can be used to modify the tape after it has been recorded.
   *
   * @tparam T_Position  Global tape position, usually chosen as Tape::Position.
   */
  template<typename T_Position, typename T_Impl>
  struct EditingTapeInterface {
    public:

      using Position = CODI_DD(T_Position, EmptyPosition);  ///< See EditingTapeInterface.
      using Impl = CODI_DD(T_Impl, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));

      /*******************************************************************************/
      /// @name Interface definition

      /// Erase a part of the tape. It has to hold start <= end.
      void erase(Position const& start, Position const& end);

      /// Copy the specified range of the source tape and append it to the end of this tape. It has to hold
      /// start <= end.
      void append(Impl& source, Position const& start, Position const& end);
  };
}
