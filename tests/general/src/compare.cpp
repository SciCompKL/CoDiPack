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

#include <sys/stat.h>

#include <iostream>
#include <string>
#include <vector>

#include "../include/compareFiles.h"
#include "../include/output.hpp"
#include "../include/testInterface.hpp"

struct CompareOutput {
  public:

    enum struct Color {
      Red = 31,
      Green = 32,
      Yellow = 33
    };

    std::string OK;
    std::string FAILURE;
    std::string FILE_MISSING;

    size_t minFieldSize;

    double threshold;

    std::vector<std::string> drivers;
    TestNames testNames;

    bool testInHeader;

    CompareOutput()
        : OK("OK"),
          FAILURE("Failure"),
          FILE_MISSING("Missing"),
          minFieldSize(0),
          threshold(1e-16),
          drivers(),
          testNames(),
          testInHeader(true) {
      listAllNames(testNames);
      minFieldSize = std::max(OK.size(), std::max(FAILURE.size(), FILE_MISSING.size()));
    }

    bool parse(int nargs, char* args[]) {
      double allOk = true;

      std::string const THRESHOLD_OPTION("-t");
      std::string const DRIVER_OPTION("-d");
      std::string const TRANSPOSE_OPTION("--trans");

      for (int curArg = 1; curArg < nargs; curArg += 1) {
        if (TRANSPOSE_OPTION == std::string(args[curArg])) {
          testInHeader = false;
        } else if (THRESHOLD_OPTION == std::string(args[curArg])) {
          curArg += 1;
          if (curArg >= nargs) {
            std::cerr << "Error: Missing value for -t option." << std::endl;
            allOk = false;
            break;
          }
          threshold = std::stod(args[curArg]);
        } else if (DRIVER_OPTION == std::string(args[curArg])) {
          curArg += 1;
          if (curArg >= nargs) {
            std::cerr << "Error: Missing value for -d option." << std::endl;
            allOk = false;
            break;
          }
          drivers.push_back(args[curArg]);
        } else {
          std::cerr << "Error: Unknown argument: " << args[curArg] << std::endl;
          allOk = false;
        }
      }

      return allOk;
    }

    template<typename List>
    void formatHeader(size_t const maxDriverSize, List& list) {
      printf("%*s", (int)maxDriverSize, " ");
      for (std::string const& item : list) {
        int maxEntrySize = (int)std::max(item.size(), minFieldSize);
        printf(" %s", formatCenter(item, maxEntrySize, item.size()).c_str());
      }
      printf("\n");
    }

    bool getLongModeName(std::string const& driverName, std::string& modeName) {
      bool allOk = true;
      size_t modePos = driverName.find("_");
      if (std::string::npos == modePos) {
        std::cerr << "Error: could not find mode in driver name: " << driverName << std::endl;
        allOk = false;
      }

      if (allOk) {
        std::string mode = driverName.substr(0, modePos);
        if ("D0" == mode) {
          modeName = "deriv0th";
        } else if ("D1" == mode) {
          modeName = "deriv1st";
        } else if ("D2" == mode) {
          modeName = "deriv2nd";
        } else {
          std::cerr << "Error: No long mode name available for: " << mode << std::endl;
          allOk = false;
        }
      }

      return allOk;
    }

    bool formatEntry(std::string const& driver, std::string const& test, int maxCellSize, std::string& cell) {
      bool allOk = true;
      std::string modeName;
      if (!getLongModeName(driver, modeName)) {
        return false;
      }

      std::string baseFile = generateTestCompareFileName(test, modeName);
      std::string resultFile = generateDriverOutputFileName(test, driver);

      int contentSize = 0;
      std::string result = "";

      if (isTestAvail(resultFile)) {
        bool same = compareFiles(baseFile, resultFile, threshold);
        if (same) {
          contentSize = OK.size();
          result = formatColor(OK, Color::Green);
        } else {
          contentSize = FAILURE.size();
          result = formatColor(FAILURE, Color::Red);
          allOk = false;
        }
      } else {
        contentSize = FILE_MISSING.size();
        result = formatColor(FILE_MISSING, Color::Yellow);
      }

      int targetSize = std::max((int)minFieldSize, maxCellSize);
      cell = " " + formatCenter(result, targetSize, contentSize);

      return allOk;
    }

    bool run() {
      bool allOk = true;

      size_t maxDriverSize = getMaxSize(drivers);
      size_t maxTestSize = getMaxSize(testNames);

      std::string curLine;
      if (testInHeader) {
        formatHeader(maxDriverSize + 1, testNames);

        for (std::string const& curDriver : drivers) {
          curLine = "";
          curLine += format("%*s:", (int)maxDriverSize, curDriver.c_str());

          for (std::string const& curTest : testNames) {
            std::string cell;
            allOk &= formatEntry(curDriver, curTest, (int)curTest.size(), cell);
            curLine += cell;
          }
          std::cout << curLine << std::endl;
        }
      } else {
        formatHeader(maxTestSize + 1, drivers);

        for (std::string const& curTest : testNames) {
          curLine = "";
          curLine += format("%*s:", (int)maxTestSize, curTest.c_str());

          for (std::string const& curDriver : drivers) {
            std::string cell;
            allOk &= formatEntry(curDriver, curTest, (int)curDriver.size(), cell);
            curLine += cell;
          }
          std::cout << curLine << std::endl;
        }
      }

      return allOk;
    }

    std::string formatCenter(std::string const& text, int size, int contentSize) {
      int pad = size - contentSize;
      int leftPad = pad / 2 + pad % 2;
      int rightPad = pad / 2;

      return format("%*s%s%*s", leftPad, "", text.c_str(), rightPad, "");
    }

    std::string formatColor(std::string const& text, Color color) {
      return format("\033[%dm%s\033[0m", (int)color, text.c_str());
    }

    template<typename List>
    size_t getMaxSize(List const& list) {
      size_t maxSize = 0;
      for (std::string const& item : list) {
        maxSize = std::max(maxSize, item.size());
      }

      return maxSize;
    }

    std::string generateDriverOutputFileName(std::string const& test, std::string const& driver) {
      std::string file = "build/results";
      file += "/" + driver;
      file += "/" + test + ".out";

      return file;
    }

    std::string generateTestCompareFileName(std::string const& test, std::string const& mode) {
      std::string file = "results";
      file += "/" + mode;
      file += "/" + test + ".out";

      return file;
    }

    bool isTestAvail(std::string const& file) {
      struct stat sb;

      return !stat(file.c_str(), &sb);
    }
};

int main(int nargs, char* args[]) {
  CompareOutput compare;

  bool allOk = true;

  allOk = compare.parse(nargs, args);
  if (allOk) {
    allOk = compare.run();
  }

  if (allOk) {
    return 0;
  } else {
    return -1;
  }
}
