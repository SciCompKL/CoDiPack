#pragma once

#include <map>
#include <set>
#include <string>
#include <vector>

#include "testMacros.hpp"

struct TestInterface {
  public:

    virtual ~TestInterface() {}

    virtual double getEvalPoint(int point, int col) = 0;
    virtual int getEvalPointsCount() = 0;
    virtual int getInputCount() = 0;
    virtual std::string getName() = 0;
    virtual int getOutputCount() = 0;

    template<typename Number>
    static void func(Number* x, Number* y);
};

template<typename Number>
using TestFunc = void(*)(Number* x, Number* u);

template<typename Number>
struct TestInfo {

    TestInterface* test;
    TestFunc<Number> func;

    TestInfo() = default;
    TestInfo(TestInterface* test, TestFunc<Number> func) :
      test(test),
      func(func) {}
};

template<typename Number> using TestVector = std::vector<TestInfo<Number>>;
template<typename Number> using TestMap    = std::map<std::string, TestInfo<Number>>;
                          using TestNames  = std::set<std::string>;

void listAllNames(TestNames& names);
