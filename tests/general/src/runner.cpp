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

#include <libgen.h>
#include <string.h>
#include <sys/stat.h>

#include <codi/misc/exceptions.hpp>
#include <codi/misc/macros.hpp>

#include "../include/drivers/driverInterface.hpp"
#include "../include/testInterface.hpp"

#ifndef DRIVER
  #error A driver include needs to be specified
#endif

#ifndef DRIVER_NAME
  #error A driver name needs to be specified
#endif

#include DRIVER

using Driver = CODI_DECLARE_DEFAULT(DRIVER_NAME, CODI_TEMPLATE(DriverInterface<double>));

struct Runner {
  public:

    using Number = typename Driver::Number;
    Driver driver;

    Runner() : driver() {}

    void run() {
      TestVector<Number> testInfos = driver.getTestInfos();

      for (auto& curInfo : testInfos) {
        std::string outFile = generateOutputFile(curInfo.test);

        FILE* out = fopen(outFile.c_str(), "w");
        std::cout << "Running Driver: " << driver.getName() << " Test: " << curInfo.test->getName() << std::endl;
        driver.runTest(curInfo, out);
        fclose(out);
      }
    }

  private:

    int makePath(char const* dir, mode_t mode) {
      struct stat sb;

      if (!stat(dir, &sb)) {
        return 0;
      }

      makePath(dirname(strdupa(dir)), mode);

      return mkdir(dir, mode);
    }

    std::string generateOutputFile(TestInterface* test) {
      std::string file = "build/results";
      file += "/" + driver.getName();

      makePath(file.c_str(), 0755);

      file += "/" + test->getName() + ".out";

      return file;
    }
};

int main(int nargs, char** args) {
  (void)nargs;
  (void)args;

  Runner runner;

  runner.run();

  return 0;
}
