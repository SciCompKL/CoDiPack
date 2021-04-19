#pragma once

#include <stdio.h>

#include "../testInterface.hpp"

enum class DriverOrder
{
  Deriv0th,
  Deriv1st,
  Deriv2nd
};

template<typename Number>
struct DriverInterface {
  public:

    virtual ~DriverInterface() {}

    virtual std::string getName() = 0;
    virtual TestVector<Number> getTestInfos() = 0;

    virtual void runTest(TestInfo<Number>& test, FILE* out) = 0;
};
