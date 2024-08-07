/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2024 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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

#include <omp.h>

#include "../threadInformationInterface.hpp"
#include "openMPAtomic.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Thread information for OpenMP.
   */
  struct OpenMPThreadInformation : public ThreadInformationInterface {
    public:

      /// \copydoc ThreadInformationInterface::getMaxThreads()
      /// <br><br> Implementation: Limit applies to all threads, also those due to nesting.
      static CODI_INLINE int getMaxThreads() {
        return 512;
      }

      /// \copydoc ThreadInformationInterface::getThreadId()
      /// <br><br> Implementation: Returns custom IDs to account for nesting, in particular not omp_get_thread_num().
      static CODI_INLINE int getThreadId() {
        static OpenMPAtomic<int> nextThreadId = 0;

        static int myThreadId = -1;
        CODI_OMP_THREADPRIVATE(myThreadId)

        if (myThreadId == -1) {
          myThreadId = nextThreadId++;
        }

        codiAssert(myThreadId < getMaxThreads());

        return myThreadId;
      }
  };
}
