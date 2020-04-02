#pragma once

#include "../../testInterface.hpp"

#include "testCopy.hpp"
#include "testExpr.hpp"

template<typename Number>
void listTestsBasic(TestVector<Number>& tests) {
  tests.push_back(new TestCopy<Number>());
  tests.push_back(new TestExpr<Number>());
}
