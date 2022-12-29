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
#include "../interfaces/fullTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Abstracts the internal set of adjoint variables provided as part of the tape.
   *
   * This interface describes the interaction between a tape and its associated adjoint variables as well as the state
   * of the adjoint variables. Details on how the adjoint variables are implemented are abstracted away by this
   * interface.
   *
   * The adjoint variables can be read and written, resized, zeroed, and swapped. The number of adjoint variables can be
   * queried, and, if applicable, a raw pointer to an underlying array implementation can be obtained.
   *
   * The adjoint variables can be "in use" or "not in use". For example, they should be considered in use during a tape
   * evaluation, which means that no resizing should take place until the evaluation is finished. The tape is
   * responsible for setting the state of the adjoint variables accordingly.
   *
   * A tape that maintains its adjoints internally against this interface can easily exchange the adjoint
   * implementation. The principle use case of this interface is, at the moment, replacing a classical, tape-local
   * adjoint vector by a thread-safe global one for use in a shared-memory parallel setting, see LocalAdjoints and
   * ThreadSafeGlobalAdjoints.
   *
   * @tparam T_Gradient    The gradient type of a tape, usually chosen as ActiveType::Gradient.
   * @tparam T_Identifier  The adjoint/tangent identification of a tape, usually chosen as ActiveType::Identifier.
   * @tparam T_Tape        The associated tape type.
   */
  template<typename T_Gradient, typename T_Identifier, typename T_Tape>
  struct InternalAdjointsInterface {
    public:

      /// See InternalAdjointsInterface.
      using Tape = CODI_DD(T_Tape, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));
      using Gradient = CODI_DD(T_Gradient, double);   ///< See InternalAdjointsInterface.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See InternalAdjointsInterface.

      /// Constructor
      /// @param initialSize  Initial number of adjoint variables.
      InternalAdjointsInterface(size_t initialSize) {
        CODI_UNUSED(initialSize);
      }

      /// Reference access to the adjoint variable identified by identifier.
      CODI_INLINE Gradient& operator[](Identifier const& identifier);

      /// Constant reference access to the adjoint variable identified by identifier.
      CODI_INLINE Gradient const& operator[](Identifier const& identifier) const;

      /// Pointer to an underlying array implementation.
      CODI_INLINE Gradient* data();

      /// Returns the number of adjoint variables.
      CODI_INLINE size_t size() const;

      /// Ensure that identifiers up to newSize can be passed to operator[] without error.
      CODI_NO_INLINE void resize(Identifier const& newSize);

      /// Set all variables with identifiers start...end-1 to zero.
      CODI_INLINE void zero(Identifier const& start, Identifier const& end);

      /// Set all adjoint variables to Gradient().
      CODI_INLINE void zeroAll();

      /// Swap two sets of adjoint variables.
      template<typename Impl>
      CODI_INLINE void swap(Impl& other);

      /// Declare that the adjoints are in use, e.g., during a tape evaluation, and cannot be resized right now.
      CODI_INLINE void beginUse();

      /// Declare that the adjoints are no longer occupied.
      CODI_INLINE void endUse();
  };

}
