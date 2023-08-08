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
   * @brief Provides information on threads.
   */
  struct ThreadInformationInterface {
    public:
      /**
       * @brief Provides an upper bound on the number of threads.
       */
      static CODI_INLINE int getMaxThreads();

      /**
       * @brief Returns the id of the calling thread.
       *
       * Thread ids are integers 0, 1, ..., maximum number of threads - 1.
       */
      static CODI_INLINE int getThreadId();
  };

  /**
   * @brief Default implementation of ThreadInformationInterface for serial applications.
   */
  struct DefaultThreadInformation : public ThreadInformationInterface {
    public:
      /// \copydoc ThreadInformationInterface::getMaxThreads
      static CODI_INLINE int getMaxThreads() {
        return 1;
      }

      /// \copydoc ThreadInformationInterface::getThreadId
      static CODI_INLINE int getThreadId() {
        return 0;
      }
  };
}
