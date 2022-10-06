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

#include <algorithm>
#include <vector>

#include "interalAdjointsInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Adjoint variables owned by a tape instance.
   *
   * @tparam T_Gradient    The gradient type of a tape, usually chosen as ActiveType::Gradient.
   * @tparam T_Identifier  The adjoint/tangent identification of a tape, usually chosen as ActiveType::Identifier.
   * @tparam T_Tape        The associated tape type.
   */
  template<typename T_Gradient, typename T_Identifier, typename T_Tape>
  struct LocalAdjoints : public InternalAdjointsInterface<T_Gradient, T_Identifier, T_Tape> {
    public:

      using Tape = CODI_DD(T_Tape, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>)); ///< See LocalAdjoints.
      using Gradient = CODI_DD(T_Gradient, double);   ///< See LocalAdjoints.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See LocalAdjoints.

    private:

      bool inUse;

      std::vector<Gradient> adjoints;  ///< Vector of adjoint variables.

    public:

      /// Constructor
      LocalAdjoints(size_t initialSize) : InternalAdjointsInterface<Gradient, Identifier, Tape>(initialSize),
                                          inUse(false),
                                          adjoints(initialSize) {}

      /// See InternalAdjointsInterface.
      CODI_INLINE Gradient& operator[](Identifier const& identifier) {
        return adjoints[(size_t)identifier];
      }

      /// See InternalAdjointsInterface.
      CODI_INLINE Gradient const& operator[](Identifier const& identifier) const {
        return adjoints[(size_t)identifier];
      }

      /// See InternalAdjointsInterface.
      CODI_INLINE Gradient* data() {
        return adjoints.data();
      }

      /// See InternalAdjointsInterface.
      CODI_INLINE size_t size() const {
        return adjoints.size();
      }

      /// See InternalAdjointsInterface.
      CODI_NO_INLINE void resize(Identifier const& newSize) {
        if (inUse) {
          CODI_EXCEPTION("Cannot resize adjoints while they are in use.");
        }

        adjoints.resize((size_t)newSize);
      }

      /// See InternalAdjointsInterface.
      CODI_INLINE void zeroAll() {
        for (Gradient& gradient : adjoints) {
          gradient = Gradient();
        }
      }

      /// See InternalAdjointsInterface.
      CODI_INLINE void swap(LocalAdjoints& other) {
        std::swap(adjoints, other.adjoints);
      }

      /// See InternalAdjointsInterface.
      CODI_INLINE void beginUse() {
        inUse = true;
      }

      /// See InternalAdjointsInterface.
      CODI_INLINE void endUse() {
        inUse = false;
      }

  };
}
