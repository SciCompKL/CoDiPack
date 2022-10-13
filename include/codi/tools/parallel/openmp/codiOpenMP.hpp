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

#include "../../../../codi.hpp"
#include "../../../expressions/parallelActiveType.hpp"
#include "../../../tapes/indices/parallelReuseIndexManager.hpp"
#include "../../../tapes/misc/threadSafeGlobalAdjoints.hpp"
#include "../../helpers/openMPExternalFunctionHelper.hpp"
#include "openMPAtomic.hpp"
#include "openMPMutex.hpp"
#include "openMPStaticThreadLocalPointer.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Parallel toolbox for OpenMP.
  using OpenMPToolbox = ParallelToolbox<OpenMPAtomic, OpenMPMutex, OpenMPStaticThreadLocalPointer>;

  /// Thread-safe global adjoints for OpenMP.
  template<typename Gradient, typename Identifier, typename Tape>
  using OpenMPGlobalAdjoints = ThreadSafeGlobalAdjoints<Gradient, Identifier, Tape, OpenMPToolbox>;

  /// \copydoc codi::RealReverseIndexGen <br><br>
  /// This a thread-safe implementation for use with OpenMP.
  template<typename Real, typename Gradient = OpenMPAtomic<Real>,
           typename IndexManager = ParallelReuseIndexManager<int, OpenMPToolbox>>
  using RealReverseIndexOpenMPGen = ParallelActiveType<
                                      JacobianReuseTape<
                                        JacobianTapeTypes<Real, Gradient, IndexManager, DefaultChunkedData,
                                          OpenMPGlobalAdjoints>>,
                                      OpenMPToolbox>;

  /// \copydoc codi::RealReverseIndexOpenMPGen
  using RealReverseIndexOpenMP = RealReverseIndexOpenMPGen<double>;
}
