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

#include <codi/misc/fileSystem.hpp>

#include "applications/transport1D.hpp"

using codi::algorithms::ApplicationHints;
using codi::algorithms::ApplicationHintsFlags;
using codi::algorithms::ReverseAccumulationSettings;
using codi::algorithms::CheckpointHandle;

char const* const OUTPUT_DIR = "testTransport1D";
char const* const CHECKPOINT_DIR = "testTransport1D_checkpoints";

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
  {"InitRecord", ApplicationHintsFlags::InitializationComputesP | ApplicationHintsFlags::PStateIsAvailable | ApplicationHintsFlags::FComputationIsAvailable},
  {"InitRecompute_PIterableYes", ApplicationHintsFlags::PComputationIsAvailable | ApplicationHintsFlags::PStateIsAvailable | ApplicationHintsFlags::FComputationIsAvailable},
  {"InitRecompute_PIterableNo", ApplicationHintsFlags::PComputationIsAvailable | ApplicationHintsFlags::FComputationIsAvailable},
  {"InitRecord_FComputeNo", ApplicationHintsFlags::InitializationComputesP | ApplicationHintsFlags::PStateIsAvailable}
};

AppConfig defaultAppConfig =
    {"", ApplicationHints::NONE()};

template<typename Type>
void prepare(Transport1D<Type>& app, std::string const& folder, std::string file) {
  FileSystem::makePath(folder.c_str());
  app.setOutputFolder(folder);
  app.setOutputFile(folder + "/" + file);
  app.initialize();

  std::cout << "Running '" << folder << "'" << std::endl;
}

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

std::vector<AppConfig> selectAll(AppConfig const& config) {
  std::vector<AppConfig> s(1);
  s[0] = config;

  return s;
}

template<template <typename> class Algo, typename Type, typename Settings>
void runProblem(VectorConfig<Type> const& vecConf, std::vector<AppConfig> const& appConf, Settings& settings, std::string const& prefix, bool onlyWriteFinal = true) {

  Transport1D<Type> app;
  app.getCheckpointInterface()->setFolder(CHECKPOINT_DIR);
  app.getIOInterface()->onlyWriteFinal = onlyWriteFinal;

  for(AppConfig const& curAppConfig: appConf) {
    std::string name = codi::StringUtil::format("%s/%s", OUTPUT_DIR, prefix.c_str());
    if(0 != curAppConfig.name.size()) {
      name += "_" + curAppConfig.name;
    }
    if(0 != vecConf.name.size()) {
      name += "_" + vecConf.name;
    }

    app.setHints(curAppConfig.hints);
    app.settings.functionalNumber = vecConf.vectorFunctions;
    prepare(app, name, "run.out");

    Algo<Transport1D<Type>> algo(settings);
    algo.run(app);
  }
}

template<template <typename> class Algo, typename Type, typename Settings>
void runProblem(Settings& settings, std::string const& prefix, bool onlyWriteFinal = true) {
  runProblem<Algo>(VectorConfig<Type>("", 1), selectAll(defaultAppConfig), settings, prefix, onlyWriteFinal);
}

void createBasicCheckpoint() {
  codi::algorithms::PrimalEvaluationSettings settings;
  settings.checkRelConvergence = false;
  settings.absThreshold = 0.000000001;
  settings.writeFinalCheckpoint = true;
  runProblem<codi::algorithms::PrimalEvaluation, double>(settings, "primal");
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

void runBBTests()
{
  codi::algorithms::BlackBoxSettings settings;
  settings.checkRelConvergence = false;
  settings.absThreshold = 0.000000001;

  runProblem<codi::algorithms::BlackBox, codi::RealReverse>(settings, "blackBox");
}

void runBBWCTests()
{

  // Setup checkpoints
  codi::algorithms::PrimalEvaluationSettings checkSettings;
  checkSettings.checkRelConvergence = false;
  checkSettings.absThreshold = 0.000000001;
  checkSettings.writeCheckpoints = true;
  runProblem<codi::algorithms::PrimalEvaluation, double>(checkSettings, "blackBoxWithCheck_writeCheck");

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

void runCheckpointTest() {
  codi::algorithms::CheckpointTestSettings settings{};

  runProblem<codi::algorithms::CheckpointTest, codi::RealReverse>(settings, "checkpointTest");
}

void runFiniteDifferenceEvaluationTest()
{
  codi::algorithms::FiniteDifferenceEvaluationSettings settings{};
  settings.fullJacobian = true;
  settings.maxIterations = 455;
  settings.primalValidationThreshold = 1e-10;
  settings.relativeStepSize = true;
  settings.stepSizes = {1e-2};

  runProblem<codi::algorithms::FiniteDifferenceEvaluation, double>(settings, "finiteDifferenceEvaluation");
}

void runForwardTests()
{
  codi::algorithms::ForwardModeSettings settings;
  settings.maxIterations = 455;

  runProblem<codi::algorithms::ForwardMode>(VectorConfig<codi::RealForward>("Vec1", 1), select(appConfigs, {1}), settings, "forward");
  runProblem<codi::algorithms::ForwardMode>(VectorConfig<codi::RealForwardVec<2>>("Vec2", 2), select(appConfigs, {1}), settings, "forward");

  settings.fullJacobian = true;
  settings.primalValidationThreshold = 1e-10;
  runProblem<codi::algorithms::ForwardMode, codi::RealForward>(settings, "forwardFullJacobian");
}

int main(int nargs, char** args) {
  (void)nargs;
  (void)args;

  FileSystem::makePath(OUTPUT_DIR);
  FileSystem::makePath(CHECKPOINT_DIR);

  createBasicCheckpoint();

  runCheckpointTest();

  runRATests();

  runBBTests();

  runBBWCTests();

  runForwardTests();

  runFiniteDifferenceEvaluationTest();

  return 0;
}
