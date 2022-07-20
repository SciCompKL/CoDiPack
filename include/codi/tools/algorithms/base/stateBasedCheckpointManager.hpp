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

#include <dirent.h>  // Unix only

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <regex>
#include <string>
#include <vector>

#include "../../../config.h"
#include "../../../expressions/lhsExpressionInterface.hpp"
#include "../../../misc/fileIo.hpp"
#include "../../../misc/macros.hpp"
#include "../../../misc/stringUtil.hpp"
#include "../interfaces/applicationInterface.hpp"
#include "../interfaces/checkpointManagerInterface.hpp"
#include "../interfaces/fileIOInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  namespace algorithms {

    template<typename T_PassiveReal>
    struct StateVectorCheckpoint : public codi::algorithms::CheckpointBase {
      public:

        using PassiveReal = CODI_DD(T_PassiveReal, double);

        std::vector<PassiveReal> data;

        bool isWritten;
        bool isListed;

        StateVectorCheckpoint(int iteration)
            : codi::algorithms::CheckpointBase(iteration), data(0), isWritten(false), isListed(false) {}
    };

    /// Assumes setIteration method is available on Application.
    /// Uses app io to write checkpoint data
    template<typename T_Type, typename T_FileIO, typename T_Application>
    struct StateBasedCheckpointManager : public codi::algorithms::CheckpointManagerInterface {
      public:
        using Type = CODI_DD(T_Type, CODI_T(codi::LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
        using FileIO = CODI_DD(T_FileIO, FileIOInterface);
        using Application = CODI_DD(T_Application, CODI_T(codi::algorithms::ApplicationInterface<CODI_ANY>));

        using PassiveReal = RealTraits::PassiveReal<Type>;
        using Checkpoint = StateVectorCheckpoint<PassiveReal>;

        char const* const CHECKPOINT_NAME = "checkpoint";

      private:

        Application* app;
        FileIO* io;

      public:

        std::string folder;

      public:

        StateBasedCheckpointManager(std::string const& folder, Application* app, FileIO* io)
            : app(app), io(io), folder(folder) {}

        void setFolder(std::string const& value) {
          folder = value;
        }

        std::string getFolder() {
          return folder;
        }

        Checkpoint* create() {
          Checkpoint* cp = new Checkpoint(app->getIteration());
          cp->data.resize(app->getSizeY());

          app->iterateY([&](Type& value, size_t pos) { cp->data[pos] = codi::RealTraits::getPassiveValue(value); });

          return cp;
        }

        std::vector<codi::algorithms::CheckpointHandle*> list() {
          std::vector<codi::algorithms::CheckpointHandle*> checkList;

          std::regex checkpointRegex(StringUtil::format("%s_(\\d+).%s", CHECKPOINT_NAME, io->getFileEnding().c_str()));
          std::smatch match;

          struct dirent* drnt = NULL;

          DIR* dir = opendir(folder.c_str());
          if (dir) {
            while (NULL != (drnt = readdir(dir))) {
              if (0 != strcmp(drnt->d_name, ".") && 0 != strcmp(drnt->d_name, "..")) {
                std::string dirName = drnt->d_name;
                if (std::regex_match(dirName, match, checkpointRegex)) {
                  int iter = std::stoi(match[1].str().c_str());
                  Checkpoint* check = new Checkpoint(iter);
                  check->isListed = true;
                  check->isWritten = true;
                  checkList.push_back(check);
                }
              }
            }
            closedir(dir);
          } else {
            std::cerr << "Can not open directory '" << folder << "'" << std::endl;
          }

          return checkList;
        }

        void load(codi::algorithms::CheckpointHandle* cp) {
          Checkpoint* check = cast(cp);

          bool clearData = false;
          if (0 == check->data.size()) {
            read(check);
            clearData = true;
          }

          app->iterateY([=](Type& value, size_t pos) { value = check->data[pos]; });

          app->setIteration(check->getIteration());

          if (clearData) {
            check->data.resize(0);
          }
        }

        void remove(codi::algorithms::CheckpointHandle* cp) {
          Checkpoint* check = cast(cp);

          if (check->isWritten && !check->isListed) {
            std::remove(createFileName(check).c_str());
          }

          delete check;
        }

        void write(codi::algorithms::CheckpointHandle* cp) {
          Checkpoint* check = cast(cp);

          size_t totalSize = sizeof(PassiveReal) * check->data.size() + sizeof(size_t);

          typename FileIO::WriteHandle h = io->openWrite(createFileName(check), totalSize);
          io->write(h, check->data.size());
          io->write(h, check->data.data(), check->data.size());
          io->closeWrite(h);

          check->isWritten = true;
          check->data.resize(0);
        }

        void read(codi::algorithms::CheckpointHandle* cp) {
          Checkpoint* check = cast(cp);

          size_t size = 0;
          typename FileIO::ReadHandle h = io->openRead(createFileName(check));

          io->read(h, size);
          check->data.resize(size);

          io->read(h, check->data.data(), check->data.size());
          io->closeRead(h);
        }

      protected:

        std::string createFileName(Checkpoint* check) {
          return StringUtil::format("%s/%s_%05d.%s", folder.c_str(), CHECKPOINT_NAME, check->getIteration(),
                                    io->getFileEnding().c_str());
        }

      private:

        Checkpoint* cast(codi::algorithms::CheckpointHandle* cp) {
          return static_cast<Checkpoint*>(cp);
        }
    };
  }
}
