#pragma once

#include <map>
#include <set>
#include <string>
#include <vector>

#include "testMacros.hpp"

template<typename Number>
struct TestInterface {
  public:

    virtual ~TestInterface() {}

    virtual double getEvalPoint(int point, int col) = 0;
    virtual int getEvalPointsCount() = 0;
    virtual int getInputCount() = 0;
    virtual std::string getName() = 0;
    virtual int getOutputCount() = 0;

    virtual void func(Number* x, Number* y) = 0;

};

template<typename Number> using TestVector = std::vector<TestInterface<Number>*>;
template<typename Number> using TestMap    = std::map<std::string, TestInterface<Number>*>;
                          using TestNames  = std::set<std::string>;

void listAllNames(TestNames& names);
