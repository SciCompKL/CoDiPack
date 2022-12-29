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
using TestFunc = void (*)(Number* x, Number* y);

template<typename Number>
struct TestInfo {
    TestInterface* test;
    TestFunc<Number> func;

    TestInfo() = default;
    TestInfo(TestInterface* test, TestFunc<Number> func) : test(test), func(func) {}
};

template<typename Number>
using TestVector = std::vector<TestInfo<Number>>;
template<typename Number>
using TestMap = std::map<std::string, TestInfo<Number>>;
using TestNames = std::set<std::string>;

void listAllNames(TestNames& names);
