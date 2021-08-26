#pragma once

#include "../tests/allTests.hpp"
#include "driverInterface.hpp"

template<typename T_Number>
struct DriverBase : public DriverInterface<T_Number> {
  public:

    using Number = CODI_DECLARE_DEFAULT(T_Number, double);

    virtual void createAllTests(TestVector<Number>& tests) = 0;

  private:

    std::string name;

  public:

    DriverBase(std::string const& name) : name(name) {}

    std::string getName() {
      return name;
    }

    TestVector<Number> getTestInfos() {
      TestVector<Number> testInfos;

      createAllTests(testInfos);

      return testInfos;
    }

  protected:
    template<typename Number>
    void createTests(TestVector<Number>& tests) {
      (void)tests;
    }

    template<typename Number, typename Test, typename... Args>
    void createTests(TestVector<Number>& tests) {
      tests.push_back(TestInfo<Number>(new Test(), Test::template func<Number>));

      createTests<Number, Args...>(tests);
    }

    void prepare(Number* x, Number* y, int curPoint, TestInterface* test, FILE* out) {
      fprintf(out, "Point %d : {", curPoint);

      for (int i = 0; i < test->getInputCount(); ++i) {
        if (i != 0) {
          fprintf(out, ", ");
        }
        double val = test->getEvalPoint(curPoint, i);
        fprintf(out, "%f", val);

        x[i] = (Number)(val);
      }
      fprintf(out, "}\n");

      for (int i = 0; i < test->getOutputCount(); ++i) {
        y[i] = 0.0;
      }
    }
};
