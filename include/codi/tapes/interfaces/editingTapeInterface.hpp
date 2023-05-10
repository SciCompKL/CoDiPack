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

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../data/position.hpp"
#include "fullTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Edit tapes after they have been recorded.
   *
   * These interface functions can be used to modify the tape after it has been recorded. Specifically, they allow to
   * erase parts of a tape, and to append a specific range of a source tape to a destination tape.
   *
   * This interface was introduced for additional flexibility when managing multiple tapes in a shared-memory parallel
   * context. The erase function, for example, can be used to remove a preliminary recording from the tape once
   * additional information is available. The append function can be used to move recordings ending up in the wrong tape
   * to the correct one. This is only required in edge cases and most AD workflows will never make use of this
   * interface, especially if they use only a single tape. Other cases might be covered by (positional) tape resets.
   *
   * Note that tapes with a linear index management strategy (see LinearIndexManager) can't implement this interface
   * because a statement's left hand side index is implicitly encoded in the statement's position on the tape. Erasing
   * parts of a tape would produce wrong subsequent left hand side indices, and appending statements from one tape to
   * another is not meaningful because the sequences of left hand side indices are tape-specific.
   *
   * @tparam T_Position  Global tape position, usually chosen as Tape::Position.
   */
  template<typename T_Position>
  struct EditingTapeInterface {
    public:

      using Position = CODI_DD(T_Position, EmptyPosition);  ///< See EditingTapeInterface.

      /*******************************************************************************/
      /// @name Interface definition

      /// Erase a part of the tape. It has to hold start <= end.
      void erase(Position const& start, Position const& end);

      /// @brief Erase a part of the tape. It has to hold start <= end.
      /// This variant of erase takes a reference to an empty helper tape. It is used as a buffer to implement erase via
      /// reset and append while avoiding the overhead of allocating a temporary tape for each erase call. Upon
      /// returning, emptyTape is guaranteed to be empty again, in the sense of a tape reset.
      void erase(Position const& start, Position const& end, EditingTapeInterface& emptyTape);

      /// Copy the specified range of the source tape and append it to the end of this tape. It has to hold
      /// start <= end.
      void append(EditingTapeInterface& source, Position const& start, Position const& end);
  };
}
