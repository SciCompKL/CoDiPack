
#include <libgen.h>
#include <string.h>
#include <sys/stat.h>

#include <codi/aux/exceptions.hpp>
#include <codi/aux/macros.h>


#include "../include/testInterface.hpp"
#include "../include/drivers/driverInterface.hpp"

#ifndef DRIVER
  #error A driver include needs to be specified
#endif

#ifndef DRIVER_NAME
  #error A driver name needs to be specified
#endif

#include DRIVER

using Driver = DECLARE_DEFAULT(DRIVER_NAME, TEMPLATE(DriverInterface<double>));


struct Runner {

    using Number = typename Driver::Number;
    Driver driver;

    std::map<DriverOrder, std::string> orderName;

    Runner() : driver(), orderName() {
      orderName[DriverOrder::Primal] = "primal";
      orderName[DriverOrder::Deriv1st] = "deriv1st";
      orderName[DriverOrder::Deriv2nd] = "deriv2nd";
    }

    void run() {
      TestVector<Number> testInfos = driver.getTestInfos();

      for(auto& curInfo : testInfos) {
        std::string outFile = generateOutputFile(curInfo.test);

        FILE* out = fopen(outFile.c_str(), "w");
        std::cout << "Running Driver: " << driver.getName() << " Test: " << curInfo.test->getName() << std::endl;
        driver.runTest(curInfo, out);
        fclose(out);
      }
    }

  private:

    int makePath(char const* dir, mode_t mode)
    {
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

    std::string getOrderName(DriverOrder const order) {
      auto pos = orderName.find(order);

      if(pos == orderName.end()) {
        codi::CODI_EXCEPTION("Missing name for driver order.");
      }

      return pos->second;
    }
};

int main(int nargs, char** args) {
  (void)nargs;
  (void)args;

  Runner runner;

  runner.run();

  return 0;
}
