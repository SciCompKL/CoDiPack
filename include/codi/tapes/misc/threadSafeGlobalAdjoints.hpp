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

#include "../../tools/parallel/parallelToolbox.hpp"
#include "internalAdjointsInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Provides global adjoint variables owned by a tape type. Thread-safe for use in parallel taping.
   *
   * @tparam T_Gradient         The gradient type of a tape, usually chosen as ActiveType::Gradient.
   * @tparam T_Identifier       The adjoint/tangent identification of a tape, usually chosen as ActiveType::Identifier.
   * @tparam T_Tape             The associated tape type.
   * @tparam T_ParallelToolbox  The parallel toolbox used in the associated tape. See codi::ParallelToolbox.
   */
  template<typename T_Gradient, typename T_Identifier, typename T_Tape, typename T_ParallelToolbox>
  struct ThreadSafeGlobalAdjoints : public InternalAdjointsInterface<T_Gradient, T_Identifier, T_Tape> {
    public:

      /// See ThreadSafeGlobalAdjoints.
      using Tape = CODI_DD(T_Tape, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));
      using Gradient = CODI_DD(T_Gradient, double);   ///< See ThreadSafeGlobalAdjoints.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See ThreadSafeGlobalAdjoints.
      /// See ThreadSafeGlobalAdjoints.
      using ParallelToolbox = CODI_DD(T_ParallelToolbox, CODI_DEFAULT_PARALLEL_TOOLBOX);

      using ReadWriteMutex = typename ParallelToolbox::ReadWriteMutex;  ///< See ParallelToolbox.
      using LockForUse = typename ParallelToolbox::LockForRead;         ///< See ParallelToolbox.
      using LockForRealloc = typename ParallelToolbox::LockForWrite;    ///< See ParallelToolbox.

    private:

      static std::vector<Gradient> adjoints;  ///< Vector of adjoint variables.

      /// @brief Protects adjoints.
      /// Read lock locks for using the adjoint vector. Write lock locks for reallocating it.
      static ReadWriteMutex adjointsMutex;

    public:

      /// Constructor
      ThreadSafeGlobalAdjoints(size_t initialSize)
          : InternalAdjointsInterface<Gradient, Identifier, Tape>(initialSize) {}

      /// \copydoc InternalAdjointsInterface::operator[](Identifier const&) <br><br>
      /// Implementation: No locking is performed, beginUse and endUse have to be used accordingly.
      CODI_INLINE Gradient& operator[](Identifier const& identifier) {
        return adjoints[(size_t)identifier];
      }

      /// \copydoc InternalAdjointsInterface::operator[](Identifier const&) const <br><br>
      /// Implementation: No locking is performed, beginUse and endUse have to be used accordingly.
      CODI_INLINE Gradient const& operator[](Identifier const& identifier) const {
        return adjoints[(size_t)identifier];
      }

      /// \copydoc InternalAdjointsInterface::data
      CODI_INLINE Gradient* data() {
        LockForUse lock(adjointsMutex);
        return adjoints.data();
      }

      /// \copydoc InternalAdjointsInterface::size
      CODI_INLINE size_t size() const {
        LockForUse lock(adjointsMutex);
        return adjoints.size();
      }

      /// \copydoc InternalAdjointsInterface::resize
      CODI_NO_INLINE void resize(Identifier const& newSize) {
        LockForRealloc lock(adjointsMutex);
        adjoints.resize((size_t)newSize);
      }

      /// \copydoc InternalAdjointsInterface::zeroAll
      CODI_INLINE void zeroAll() {
        for (Gradient& gradient : adjoints) {
          gradient = Gradient();
        }
      }

      /// \copydoc InternalAdjointsInterface::swap
      CODI_INLINE void swap(ThreadSafeGlobalAdjoints&) {
        /* Adjoints in this implementation are a static global member. Therefore, there is no need to swap them. */
      }

      /// \copydoc InternalAdjointsInterface::beginUse <br><br>
      /// Implementation: Sets an internal lock.
      CODI_INLINE void beginUse() {
        adjointsMutex.lockRead();
      }

      /// \copydoc InternalAdjointsInterface::endUse <br><br>
      /// Implementation: Unsets an internal lock.
      CODI_INLINE void endUse() {
        adjointsMutex.unlockRead();
      }
  };

  template<typename Gradient, typename Identifier, typename Tape, typename ParallelToolbox>
  std::vector<CODI_DD(Gradient, double)>
      ThreadSafeGlobalAdjoints<Gradient, Identifier, Tape, ParallelToolbox>::adjoints(1);

  template<typename Gradient, typename Identifier, typename Tape, typename ParallelToolbox>
  typename CODI_DD(ParallelToolbox, CODI_DEFAULT_PARALLEL_TOOLBOX)::ReadWriteMutex
      ThreadSafeGlobalAdjoints<Gradient, Identifier, Tape, ParallelToolbox>::adjointsMutex;
}
