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

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Provides basic synchronization facilities.
   */
  struct SynchronizationInterface {
    public:
      /**
       * @brief Ensures that only one among the calling threads calls the given function object.
       */
      template<typename FunctionObject>
      static CODI_INLINE void serialize(FunctionObject const& func);

      /**
       * @brief Does not return until called by all threads.
       */
      static CODI_INLINE void synchronize();
  };

  /**
   * @brief Default implementation of SynchronizationInterface for serial applications.
   */
  struct DefaultSynchronization : public SynchronizationInterface {
    public:
      /// \copydoc SynchronizationInterface::serialize
      /// <br> Implementation: does not synchronize, just calls the function object.
      template<typename FunctionObject>
      static CODI_INLINE void serialize(FunctionObject const& func) {
        func();
      }

      /// \copydoc  SynchronizationInterface::synchronize
      /// <br> Implementation: empty.
      static CODI_INLINE void synchronize() {}
  };
}
