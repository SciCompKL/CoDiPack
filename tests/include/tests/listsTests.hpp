#pragma once

#include "../testInterface.hpp"

#include "basic/listBasic.hpp"

template<typename Number>
void listTestAll(TestVector<Number>& tests) {
  listTestsBasic(tests);
}
