
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "../include/compareFiles.h"
#include "../include/testInterface.hpp"

struct CompareOutput {
    double threshold;

    std::vector<std::string> drivers;

    CompareOutput() :
      threshold(1e-16),
      drivers() {}

    bool parse(int nargs, char* args[]) {
      double allOk = true;

      std::string const THRESHOLD_OPTION("-t");
      std::string const DRIVER_OPTION("-d");

      for(int curArg = 1; curArg < nargs; curArg += 1) {
        if(THRESHOLD_OPTION == std::string(args[curArg])) {
          curArg += 1;
          if(curArg >= nargs) {
            std::cerr << "Error: Missing value for -t option." << std::endl;
            allOk = false;
            break;
          }
          threshold = std::stod(args[curArg]);
        } else if(DRIVER_OPTION == std::string(args[curArg])) {
          curArg += 1;
          if(curArg >= nargs) {
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

    void formatHeader(size_t const maxDriverSize, TestNames& testNames) {
      printf("%*s", (int)maxDriverSize, " ");
      for(std::string const& curTest : testNames) {
        printf(" %s", curTest.c_str());
      }
      printf("\n");
    }

    bool getLongModeName(std::string const& driverName, std::string& modeName) {
      bool allOk = true;
      size_t modePos = driverName.find("_");
      if(std::string::npos == modePos) {
        std::cerr << "Error: could not find mode in driver name: " << driverName << std::endl;
        allOk = false;
      }

      if(allOk) {
        std::string mode = driverName.substr(0, modePos);
        if("D1" == mode) {
          modeName = "deriv1st";
        } else {
          std::cerr << "Error: No long mode name available for: " << mode << std::endl;
          allOk = false;
        }
      }

      return allOk;
    }

    bool run() {
      bool allOk = true;
      TestNames testNames;
      listAllNames(testNames);

      size_t maxDriverSize = getMaxDriverSize();

      formatHeader(maxDriverSize + 1, testNames);

      for(std::string const& curDriver : drivers) {
        printf("%*s:", (int)maxDriverSize, curDriver.c_str());

        std::string modeName;
        if(!getLongModeName(curDriver, modeName)) {
          break;
        }

        for(std::string const& curTest : testNames) {
          std::string baseFile = generateTestComparetFile(curTest, modeName);
          std::string resultFile = generateDriverOutputFile(curTest, curDriver);

          if(isTestAvail(resultFile)) {
            bool same = compareFiles(baseFile, resultFile, threshold);
            if(same) {
              printf(" %*s", (int)curTest.size(), "OK");
            } else {
              printf(" %*s", (int)curTest.size(), "Failure");
              allOk = false;
            }
          } else {
            printf(" %*s", (int)curTest.size(), "n/a");
          }
        }
        printf("\n");
      }

      return allOk;
    }

    size_t getMaxDriverSize() {
      size_t maxSize = 0;
      for(std::string const& curDriver : drivers) {
        maxSize = std::max(maxSize, curDriver.size());
      }

      return maxSize;
    }

    std::string generateDriverOutputFile(std::string const& test, std::string const& driver) {
      std::string file = "build/results";
      file += "/" + driver;
      file += "/" + test + ".out";

      return file;
    }

    std::string generateTestComparetFile(std::string const& test, std::string const& mode) {
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
  if(allOk) {
    allOk =  compare.run();
  }

  if(allOk) {
    return 0;
  } else {
    return -1;
  }
}

