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

#include <vector>

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../misc/tapeValues.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Indices enable the mapping of primal values to their adjoint counterparts.
   *
   * In operator overloading AD, each primal variable (e.g. \f$ w \f$) needs to be mapped to the adjoint counterpart
   * (e.g. \f$ \bar w \f$. Since the adjoint cannot be stored in the primal, an identifier (usually an index) is
   * associated with each primal. This identifier is then used to access the adjoint variable.
   *
   * The interface defines the three basic operations which can be applied to a variable: assign, copy and free.
   * For each of these operations on the primal variable, the corresponding function on the identifier needs to be
   * called.
   *
   * freeIndex() only needs to be called in destructors. If a variable is overwritten, only assign needs to be called on
   * the left hand side identifier. The index manager has to decide how the old identifier is handled.
   *
   * assignUnusedIndex() provides identifiers that have not been used after the last reset. These identifiers can be
   * used for input values of the program because the adjoint will not be overwritten by intermediate variables.
   *
   * CopyNeedsStatement is a static check if the index manager implements a copy optimization, that is, if it creates
   * a new identifier for the left hand side or copies the right hand side identifier instead. Not all index management
   * approaches admit a copy optimization.
   *
   * IsLinear indicates whether the indices are coupled to the statements of a program. The tape needs to be managed
   * accordingly.
   *
   * NeedsStaticStorage indicates whether the index manager is specific to a tape type (as opposed to a specific tape
   * instance). Depending on this setting, it is stored statically or non-statically in the tape.
   *
   * Mathematical and implementational details are explained in \ref SBG2021Index.
   *
   * @tparam T_Index  Type for the identifier, usually an integer type.
   */
  template<typename T_Index>
  struct IndexManagerInterface {
    public:

      using Index = CODI_DD(T_Index, int);  ///< See IndexManagerInterface.

      /*******************************************************************************/
      /// @name Global constants

      static Index constexpr InactiveIndex = Index(0);  ///< Default inactive index for all index managers.
      static Index constexpr InvalidIndex =
          Index(-1);  ///< Default invalid index for all index mangers (max value for unsigned types).

      /*******************************************************************************/
      /// @name Identifier handling

      static bool constexpr CopyNeedsStatement =
          CODI_UNDEFINED_VALUE;  ///< True if no copy optimization is implemented. See IndexManagerInterface.
      static bool constexpr IsLinear =
          CODI_UNDEFINED_VALUE;  ///< True if identifiers are coupled to the statements. See IndexManagerInterface.

      /// True if the index manager is specific to a tape type (and not a tape instance). See IndexManagerInterface.
      static bool constexpr NeedsStaticStorage = CODI_UNDEFINED_VALUE;

      /// @brief Call on assignment of a primal value, e.g. on `w` for `w = a + b`.
      /// @return true if new indices have been generated internally.
      template<typename Tape>
      bool assignIndex(Index& index);

      /// @brief Call on registering input values.
      /// @return true if new indices have been generated internally.
      template<typename Tape>
      bool assignUnusedIndex(Index& index);

      template<typename Tape>
      void copyIndex(Index& lhs, Index const& rhs);  ///< Call on copy of a primal value, e.g. `w = a`.

      template<typename Tape>
      void freeIndex(Index& index);  ///< Call on destruction of a primal value. Usually called from the destructor.

      void reset();  ///< Reset for a new recording.

      /*******************************************************************************/
      /// @name Misc functions

      /**
       * @brief Add storage and other information to the tape values.
       * @param[inout] values  Will only create new data entries and no new section.
       */
      void addToTapeValues(TapeValues& values) const;

      /**
       * @brief Returns the largest created index.
       *
       * This is the largest entry in the adjoint vectors created by the tapes.
       */
      Index getLargestCreatedIndex() const;
  };

  // clang-format off
  template<typename Index>
  CODI_DD(Index, int) constexpr IndexManagerInterface<Index>::InactiveIndex;
  template<typename Index>
  CODI_DD(Index, int) constexpr IndexManagerInterface<Index>::InvalidIndex;
  // clang-format on
}
