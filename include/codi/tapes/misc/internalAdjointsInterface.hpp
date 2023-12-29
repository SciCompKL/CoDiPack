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
   * The set of adjoint variables can be "in use" or "not in use". The adjoint variables are "in use" whenever there is
   * read or write access to adjoint variables, or when any general property of the set of adjoint variables such as
   * size is queried. The implementations of this interface ensure mutual exclusion between the "in use" state and
   * reallocations of the set of adjoint variables due to resizing. Resizing is only allowed if the adjoint variables
   * are "not in use". The implementations of data() and size() are expected to declare usage
   * internally, if needed. For performance reasons, the operator[]() accessors and zeroAll() don't declare usage
   * internally. Instead, the tape is responsible for this. It should declare usage by calls to beginUse() and endUse().
   * This way, multiple such calls can be safeguarded by a single usage declaration.
   *
   * The tape must not call resize() as long as it has declared usage.
   *
   * To give an example, tape evaluation involves multiple operator[]() calls. Prior to the evaluation, the tape ensures
   * that the set of adjoint variables is sufficiently large. It calls beginUse() before the evaluation and endUse()
   * after it. During the evaluation, no further resizing  of the set of adjoint variables takes place.
   *
   * See codi::DataManagementTapeInterface for a multithreading perspective on the "in use" mechanism.
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

      /// Returns the number of adjoint variables. Internally, declares usage of the adjoints.
      CODI_INLINE size_t size() const;

      /// Ensure that identifiers up to newSize can be passed to operator[] without error.
      CODI_NO_INLINE void resize(Identifier const& newSize);

      /// Set all adjoint variables to Gradient().
      CODI_INLINE void zeroAll();

      /// Swap two sets of adjoint variables. Internally, declares usage of the adjoints.
      template<typename Impl>
      CODI_INLINE void swap(Impl& other);

      /// Declare that the adjoints are in use, e.g., during a tape evaluation, and cannot be resized right now.
      CODI_INLINE void beginUse();

      /// Declare that the adjoints are no longer occupied.
      CODI_INLINE void endUse();
  };

}
