/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2024 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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

/** \copydoc codi::Namespace */
namespace codi {

  template<typename Gradient, typename Identifier, typename Tape, typename ParallelToolbox>
  struct ThreadSafeGlobalAdjoints;

  /// Traits for the internal adjoint variables maintained by the tape.
  namespace InternalAdjointVectorTraits {

    /// Whether the adjoint vector is global, that is, shared among different tapes.
    template<typename InternalAdjointType>
    struct IsGlobal : std::false_type {};

#ifndef DOXYGEN_DISABLE
    template<typename Gradient, typename Identifier, typename Tape, typename ParallelToolbox>
    struct IsGlobal<ThreadSafeGlobalAdjoints<Gradient, Identifier, Tape, ParallelToolbox>> : std::true_type {};
#endif

  }

  /// General traits of adjoint vectors.
  namespace AdjointVectorTraits {

    /**
     * @brief Trait implementation to deduce the entry type from an adjoint vector type.
     *
     * Default implementation for STL containers.
     */
    template<typename AdjointVector>
    struct GradientImplementation {
      public:
        using Gradient = typename AdjointVector::value_type;  ///< Type of adjoint vector entries.
    };

#ifndef DOXYGEN_DISABLE
    /// Specialization for references.
    template<typename AdjointVector>
    struct GradientImplementation<AdjointVector&> : public GradientImplementation<AdjointVector> {};

    /// Specialization for pointers.
    template<typename T_Gradient>
    struct GradientImplementation<T_Gradient*> {
      public:
        using Gradient = T_Gradient;
    };
#endif

    /**
     * @brief Deduce the entry type from an adjoint vector type, usually identical to the gradient type of a tape.
     *
     * @tparam AdjointVector Type of the adjoint vector.
     */
    template<typename AdjointVector>
    using Gradient = typename GradientImplementation<AdjointVector>::Gradient;
  }
}
