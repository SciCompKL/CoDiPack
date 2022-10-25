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

using codi::algorithms::ApplicationHints;
using codi::algorithms::ApplicationFlags;
using codi::algorithms::ReverseAccumulationSettings;
using codi::algorithms::CheckpointHandle;

char const* const OUTPUT_DIR = "testTransport1D";
char const* const CHECKPOINT_DIR = "testTransport1D_checkpoints";

template<typename Type>
void prepare(Transport1D<Type>& app, std::string const& folder, std::string file) {
  FileSystem::makePath(folder.c_str());
  app.setOutputFolder(folder);
  app.setOutputFile(folder + "/" + file);
  app.initialize();

  std::cout << "Running '" << folder << "'" << std::endl;
}

struct AppConfig {
    std::string name;
    ApplicationHints hints;

    AppConfig(std::string const& name, ApplicationHints hints) :
      name(name),
      hints(hints)
    {}

    AppConfig() = default;
};

template<typename T_Type>
struct VectorConfig {
    using Type = T_Type;
    std::string name;
    int vectorFunctions;

    VectorConfig(std::string const& name, int vectorFunctions) :
      name(name),
      vectorFunctions(vectorFunctions)
    {}
};

size_t constexpr CONFIG_SIZE = 4;
AppConfig appConfigs[CONFIG_SIZE] = {
  {"InitRecord", ApplicationFlags::InitializationComputesP | ApplicationFlags::PStateIsAvailable | ApplicationFlags::FComputationIsAvailable},
  {"InitRecompute_PIterableYes", ApplicationFlags::PComputationIsAvailable | ApplicationFlags::PStateIsAvailable | ApplicationFlags::FComputationIsAvailable},
  {"InitRecompute_PIterableNo", ApplicationFlags::PComputationIsAvailable | ApplicationFlags::FComputationIsAvailable},
  {"InitRecord_FComputeNo", ApplicationFlags::InitializationComputesP | ApplicationFlags::PStateIsAvailable}
};


std::vector<AppConfig> selectAll(AppConfig* configs) {
  return std::vector<AppConfig>(&configs[0], &configs[CONFIG_SIZE]);
}

std::vector<AppConfig> select(AppConfig* config, std::initializer_list<int> l) {
  std::vector<AppConfig> s(l.size());

  size_t pos = 0;
  for(int i : l) {
    s[pos] = config[i];
    pos += 1;
  }

  return s;
}

template<template <typename> class Algo, typename Type, typename Settings>
void runProblem(VectorConfig<Type> const& vecConf, std::vector<AppConfig> const& appConf, Settings& settings, std::string const& prefix) {

  Transport1D<Type> app;
  app.getCheckpointInterface()->setFolder(CHECKPOINT_DIR);

  for(AppConfig const& curAppConfig: appConf) {
    std::string name = codi::StringUtil::format("%s/%s_%s_%s", OUTPUT_DIR, prefix.c_str(), curAppConfig.name.c_str(), vecConf.name.c_str());

    app.setHints(curAppConfig.hints);
    app.settings.functionalNumber = vecConf.vectorFunctions;
    prepare(app, name, "run.out");

    Algo<Transport1D<Type>> algo(settings);
    algo.run(app);
  }
}

void runForwardTests()
{
  codi::algorithms::ForwardModeSettings settings;
  settings.maxIterations = 455;

  runProblem<codi::algorithms::ForwardMode>(VectorConfig<codi::RealForward>("Vec1", 1), select(appConfigs, {1}), settings, "forward");
  runProblem<codi::algorithms::ForwardMode>(VectorConfig<codi::RealForwardVec<2>>("Vec2", 2), select(appConfigs, {1}), settings, "forward");
}


void runRATests()
{
  codi::algorithms::ReverseAccumulationSettings settings;
  settings.start = 455;
  settings.maxIterations = 550;
  settings.checkRelConvergence = false;
  settings.absThreshold = 0.000000001;

  runProblem<codi::algorithms::ReverseAccumulation>(VectorConfig<codi::RealReverse>("TapeVec1_Functional1", 1), selectAll(appConfigs), settings, "revAcc");
  runProblem<codi::algorithms::ReverseAccumulation>(VectorConfig<codi::RealReverse>("TapeVec1_Functional2", 2), select(appConfigs, {0}), settings, "revAcc");

  runProblem<codi::algorithms::ReverseAccumulation>(VectorConfig<codi::RealReverseVec<4>>("TapeVec4_Functional1", 1), select(appConfigs, {0}), settings, "revAcc");
  runProblem<codi::algorithms::ReverseAccumulation>(VectorConfig<codi::RealReverseVec<4>>("TapeVec4_Functional4", 4), select(appConfigs, {0}), settings, "revAcc");
  runProblem<codi::algorithms::ReverseAccumulation>(VectorConfig<codi::RealReverseVec<4>>("TapeVec4_Functional5", 5), select(appConfigs, {0}), settings, "revAcc");

  runProblem<codi::algorithms::ReverseAccumulation>(VectorConfig<codi::RealReverse>("CustomVec_Functional4", 4), select(appConfigs, {0}), settings, "revAcc");
  runProblem<codi::algorithms::ReverseAccumulation>(VectorConfig<codi::RealReverse>("CustomVec_Functional5", 5), select(appConfigs, {0}), settings, "revAcc");
}

void runBBTests(Problem& app)
{
  prepare(app, codi::StringUtil::format("%s/blackBox", OUTPUT_DIR), "run.out");
  codi::algorithms::BlackBox<Problem> bb{codi::algorithms::BlackBoxSettings()};
  bb.settings.checkRelConvergence = false;
  bb.settings.absThreshold = 0.000000001;
  bb.run(app);
}

void runBBWCTests()
{

  // Setup checkpoints
  Problem app;
  app.getCheckpointInterface()->setFolder(CHECKPOINT_DIR);
  prepare(app, codi::StringUtil::format("%s/blackBoxWithCheck_writeCheck", OUTPUT_DIR), "run.out");
  codi::algorithms::PrimalEvaluation<Problem> pe{};
  pe.settings.checkRelConvergence = false;
  pe.settings.absThreshold = 0.000000001;
  pe.settings.writeCheckpoints = true;
  pe.run(app);

  codi::algorithms::BlackBoxWithCheckpointsSettings settings;
  settings.start = 0;
  settings.end = 455;
  settings.verbose = true;

  runProblem<codi::algorithms::BlackBoxWithCheckpoints>(VectorConfig<codi::RealReverse>("TapeVec1_Functional1", 1), selectAll(appConfigs), settings, "blackBoxWithCheck");

  runProblem<codi::algorithms::BlackBoxWithCheckpoints>(VectorConfig<codi::RealReverse>("TapeVec1_Functional2", 2), select(appConfigs, {0}), settings, "blackBoxWithCheck");

  runProblem<codi::algorithms::BlackBoxWithCheckpoints>(VectorConfig<codi::RealReverseVec<4>>("TapeVec4_Functional1", 1), select(appConfigs, {0}), settings, "blackBoxWithCheck");
  runProblem<codi::algorithms::BlackBoxWithCheckpoints>(VectorConfig<codi::RealReverseVec<4>>("TapeVec4_Functional4", 4), select(appConfigs, {0}), settings, "blackBoxWithCheck");
  runProblem<codi::algorithms::BlackBoxWithCheckpoints>(VectorConfig<codi::RealReverseVec<4>>("TapeVec4_Functional5", 5), select(appConfigs, {0}), settings, "blackBoxWithCheck");

  runProblem<codi::algorithms::BlackBoxWithCheckpoints>(VectorConfig<codi::RealReverse>("CustomVec_Functional4", 4), select(appConfigs, {0}), settings, "blackBoxWithCheck");
  runProblem<codi::algorithms::BlackBoxWithCheckpoints>(VectorConfig<codi::RealReverse>("CustomVec_Functional5", 5), select(appConfigs, {0}), settings, "blackBoxWithCheck");
}

int main(int nargs, char** args) {
  (void)nargs;
  (void)args;

  FileSystem::makePath(OUTPUT_DIR);
  FileSystem::makePath(CHECKPOINT_DIR);

  Problem app;
  app.getCheckpointInterface()->setFolder(CHECKPOINT_DIR);
  prepare(app, codi::StringUtil::format("%s/primal", OUTPUT_DIR), "run.out");
  codi::algorithms::PrimalEvaluation<Problem> pe{codi::algorithms::PrimalEvaluationSettings()};
  pe.settings.checkRelConvergence = false;
  pe.settings.absThreshold = 0.000000001;
  pe.run(app);

  codi::algorithms::CheckpointManagerInterface* cm = app.getCheckpointInterface();
  codi::algorithms::CheckpointHandle* checkpoint = cm->create();
  cm->write(checkpoint);

  runRATests();

  runBBTests(app);

  runBBWCTests();

  runForwardTests();

  prepare(app, codi::StringUtil::format("%s/checkpointTest", OUTPUT_DIR), "run.out");
  app.setIteration(0);
  app.initialize();
  app.getIOInterface()->onlyWriteFinal = false;
  codi::algorithms::CheckpointTest<Problem> ct{codi::algorithms::CheckpointTestSettings()};
  ct.run(app);

  return 0;
}
