#include <string>
#include <iostream>
#include <sstream>

#include <sys/types.h>
#include <unistd.h>

#include "../../testInterface.hpp"

struct TestIO : public TestInterface {
  public:
    NAME("IO")
    IN(1)
    OUT(1)
    POINTS(1) = {{1.0}};

    template<typename Number>
    static void func(Number* x, Number* y) {
      y[0] = x[0];

    #if REVERSE_TAPE
      auto& tape = Number::getGlobalTape();
      std::stringstream filename;
      filename << "test" << getpid() << ".tape";

      tape.writeToFile(filename.str());
      tape.deleteData();
      tape.readFromFile(filename.str());

      unlink(filename.str().c_str());
    #endif
    }
};
