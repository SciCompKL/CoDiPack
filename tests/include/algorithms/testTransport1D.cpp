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

#include "../utils/fileSystem.hpp"
#include "applications/transport1D.hpp"

using Real = codi::RealReverse;
using Problem = Transport1D<Real>;

void prepare(Problem& app, std::string const& folder, std::string file) {
  FileSystem::makePath(folder.c_str());
  app.setOutputFolder(folder);
  app.setOutputFile(folder + "/" + file);
}

int main(int nargs, char** args) {
  (void)nargs;
  (void)args;

  FileSystem::makePath("testTransport1D");

  Problem app;
  prepare(app, "testTransport1D/primal", "run.out");
  codi::algorithms::PrimalEvaluation<Problem> pe{codi::algorithms::PrimalEvaluationSettings()};
  pe.settings.checkRelConvergence = false;
  pe.settings.absThreshold = 0.000000001;
  pe.run(app);

  prepare(app, "testTransport1D/revAcc", "run.out");
  codi::algorithms::ReverseAccumulation<Problem> ra{codi::algorithms::ReverseAccumulationSettings()};
  ra.settings.checkRelConvergence = false;
  ra.settings.absThreshold = 0.000000001;
  ra.run(app);

  prepare(app, "testTransport1D/blackBox", "run.out");
  app.setIteration(0);
  codi::algorithms::BlackBox<Problem> bb{codi::algorithms::BlackBoxSettings()};
  bb.settings.checkRelConvergence = false;
  bb.settings.absThreshold = 0.000000001;
  bb.run(app);

  prepare(app, "testTransport1D/checkpointTest", "run.out");
  app.setIteration(0);
  app.initialize();
  app.getIOInterface()->onlyWriteFinal = false;
  codi::algorithms::CheckpointTest<Problem> ct{codi::algorithms::CheckpointTestSettings()};
  ct.run(app);

  return 0;
}
