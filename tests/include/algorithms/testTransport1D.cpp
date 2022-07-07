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

#include "applications/transport1D.hpp"

#include "../utils/fileSystem.hpp"

int main(int nargs, char** args) {
  (void)nargs;
  (void)args;

  using Real = codi::RealReverse;
  using Problem = Transport1D<Real>;

  FileSystem::makePath("testTransport1D");

  Problem app("testTransport1D/primal.out");
  app.generalSettings.outputDir = "testTransport1D/primal";
  codi::algorithms::PrimalEvaluation<Problem> pe{codi::algorithms::PrimalEvaluationSettings()};
  pe.settings.checkRelConvergence = false;
  pe.settings.absThreshold = 0.000000001;
  pe.run(app);

  app.setOutputFile("testTransport1D/revAcc.out");
  app.generalSettings.outputDir = "testTransport1D/revAcc";
  codi::algorithms::ReverseAccumulation<Problem> ra{codi::algorithms::ReverseAccumulationSettings()};
  ra.settings.checkRelConvergence = false;
  ra.settings.absThreshold = 0.000000001;
  ra.run(app);

  app.setOutputFile("testTransport1D/blackBox.out");
  app.setIteration(0);
  app.generalSettings.outputDir = "testTransport1D/blackBox";
  codi::algorithms::BlackBox<Problem> bb{codi::algorithms::BlackBoxSettings()};
  bb.settings.checkRelConvergence = false;
  bb.settings.absThreshold = 0.000000001;
  bb.run(app);

  return 0;
}
