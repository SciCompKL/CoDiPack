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
#include <stdio.h>
#include <sys/stat.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct Settings {
  public:
    double threshold;

    std::vector<std::string> fileNames;
};

enum class DerivCase {
  D0,
  D1,
  D2,
  NONE
};

class ResultDiff {
  private:
    size_t fileCount;
    const Settings& settings;

    std::vector<std::ifstream> files;

    std::vector<std::string> nextLine;

    const std::string POINT_PREFIX;
    const std::string IN_PREFIX;
    const std::string OUT_PREFIX;

  public:
    ResultDiff(const Settings& settings)
        : fileCount(settings.fileNames.size()),
          settings(settings),
          files(settings.fileNames.size()),
          nextLine(settings.fileNames.size()),
          POINT_PREFIX("Point"),
          IN_PREFIX("in"),
          OUT_PREFIX("out") {}

  public:
    bool run() {
      bool allOk = true;

      if (allOk) {
        allOk = openFiles();
      }
      if (allOk) {
        allOk = compareFiles();
      }

      closeFiles();

      return allOk;
    }

  private:

    bool openFiles() {
      bool allFilesAvail = true;

      for (size_t i = 0; i < fileCount; ++i) {
        nextLine[i] = "";

        struct stat buffer;
        if (0 == stat(settings.fileNames[i].c_str(), &buffer)) {
          files[i].open(settings.fileNames[i]);
        } else {
          allFilesAvail = false;
          std::cerr << "Could not find file '" << settings.fileNames[i] << "'." << std::endl;
        }
      }

      return allFilesAvail;
    }

    void closeFiles() {
      for (size_t i = 0; i < fileCount; ++i) {
        files[i].close();
      }
    }

    void trimString(std::string& s) {
      std::string::iterator it, first, last;
      std::string temp;

      for (it = s.begin(); it < s.end(); it++) {
        if (!isspace(*it)) {
          first = it;
          break;
        }
      }

      for (it = s.end(); it > s.begin(); it--) {
        if (!isspace(*(it - 1))) {
          last = it;
          break;
        }
      }

      temp.assign(first, last);

      s = temp;
    }

    template<typename Stream>
    bool readLines(std::vector<Stream>& stream, std::vector<std::string>& lines, char delim, bool skipEmpty = true) {
      bool allEndOfFile = true;
      for (size_t i = 0; i < fileCount; ++i) {
        do {
          std::getline(stream[i], lines[i], delim);
          trimString(lines[i]);

          if (stream[i].good()) {
            allEndOfFile = false;
          }
        } while (skipEmpty && lines[i].empty() && stream[i].good());
      }

      return !allEndOfFile;
    }

    bool allStringSame(std::vector<std::string> const& strings, int& fileDiff) {
      for (size_t curFile = 1; curFile < fileCount; curFile += 1) {
        if (strings[0] != strings[curFile]) {
          fileDiff = curFile;
          return false;
        }
      }

      return true;
    }

    bool allValuesSame(int& diffPos, int& fileDiff) {
      std::vector<std::istringstream> lineStreams(fileCount);
      std::vector<std::string> tokens(fileCount);
      for (size_t curFile = 0; curFile < fileCount; curFile += 1) {
        lineStreams[curFile].str(nextLine[curFile]);
      }

      for (int curToken = 0; readLines(lineStreams, tokens, ' '); ++curToken) {
        if (0 == curToken) {
          // First token is name
          if (!allStringSame(tokens, fileDiff)) {
            diffPos = curToken;
            return false;
          }
        } else {
          // compare values
          double base = parseDouble(tokens[0]);
          for (size_t curFile = 1; curFile < fileCount; curFile += 1) {
            double curValue = parseDouble(tokens[curFile]);

            if (getDeviation(base, curValue) > settings.threshold) {
              diffPos = curToken;
              fileDiff = curFile;
              return false;
            }
          }
        }
      }

      return true;
    }

    double getDeviation(const double base, const double value) {
      double diff = std::abs(base - value);

      if (0 == diff) {
        return diff;
      } else if (0.0 == base || 0.0 == value) {
        return diff;
      } else {
        return diff / std::abs(base);
      }
    }

    bool compareFiles() {
      std::cout.precision(15);
      std::cout.setf(std::ios::scientific);
      std::cout.setf(std::ios::showpos);

      bool hasDeviation = false;
      int curPoint = 0;       // Current evaluation point. Each points resets the curDataEntry
      int curDataEntry = -1;  // -1 for header, counts up for each primal, Jacobian or Hessian entry.
      int fileDiff = 0;       // For error reporting. Which file yields the difference

      DerivCase dCase = DerivCase::NONE;

      for (int curLine = 1; readLines(files, nextLine, '\n', false) && !hasDeviation; ++curLine) {
        if (0 == nextLine[0].find(POINT_PREFIX)) {  // Check for evaluation point
          curPoint += 1;
          curDataEntry = -1;
          // Evaluation points need to be the same
          if (!allStringSame(nextLine, fileDiff)) {
            hasDeviation = true;
            std::cerr << errorFileOutput(fileDiff) << ": Evaluation point differs in line " << curLine << "."
                      << std::endl;
            break;
          }
        } else if (nextLine[0].empty()) {  // Check for empty lines. In Hessian case this indicates a new matrix.
          if (!allStringSame(nextLine, fileDiff)) {
            std::cerr << errorFileOutput(fileDiff) << ": Difference in line " << curLine << "." << std::endl;
            exit(0);
          }
          // Reset data point for D2
          if (DerivCase::D2 == dCase) {
            curDataEntry = -1;
          }
        } else {  // Data case
          if (DerivCase::NONE == dCase) {
            // Determine the type of the result file: D0, D1, D2
            if (0 == nextLine[0].find(IN_PREFIX)) {
              dCase = DerivCase::D1;
            } else if (0 == nextLine[0].find(OUT_PREFIX)) {
              if (std::string::npos != nextLine[0].find(IN_PREFIX, OUT_PREFIX.size())) {
                dCase = DerivCase::D2;
              } else {
                dCase = DerivCase::D0;
              }
            } else {
              std::cerr << "Error: Could not determine derivative case in file '" << settings.fileNames[0] << "' line "
                        << curLine << "." << std::endl;
              exit(0);
            }
          }

          // Compare header for D1 and D2
          if (-1 == curDataEntry) {
            if (DerivCase::D1 == dCase || DerivCase::D2 == dCase) {
              if (!allStringSame(nextLine, fileDiff)) {
                hasDeviation = true;
                std::cerr << errorFileOutput(fileDiff) << ": Header differs for point " << curPoint << " in line "
                          << curLine << "." << std::endl;
                break;
              }
            } else {
              // D0 has no header
              curDataEntry = 0;
            }
          }

          if (-1 != curDataEntry) {
            int diffPos = 0;
            if (!allValuesSame(diffPos, fileDiff)) {
              hasDeviation = true;
              std::cerr << errorFileOutput(fileDiff) << ": Value entry differs for point " << curPoint << " in line "
                        << curLine << "." << std::endl;
              break;
            }
          }

          curDataEntry += 1;
        }
      }

      return !hasDeviation;
    }

    std::string errorFileOutput(int filePos) {
      return settings.fileNames[0] + " " + settings.fileNames[filePos];
    }

    double parseDouble(std::string const& str) {
      try {
        return std::stod(str);
      } catch (std::out_of_range const&) {
        return 0.0;
      }
    }
};

bool compareFiles(std::string const& base, std::string const& result, double threshold) {
  Settings settings;
  settings.threshold = threshold;
  settings.fileNames.push_back(base);
  settings.fileNames.push_back(result);

  ResultDiff diff(settings);

  return diff.run();
}
