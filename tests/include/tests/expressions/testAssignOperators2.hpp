/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 *  - SciComp, TU Kaiserslautern:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#include "../../testInterface.hpp"

struct TestAssignOperators2 : public TestInterface {
  public:
    NAME("AssignOperators2")
    IN(2)
    OUT(3)
    POINTS(18) =  // clang-format off
    {
      {-10.0,   -10},
      {-10.0,    -5},
      {-10.0,     5},
      {-10.0,    10},
      { -5.0,   -10},
      { -5.0,    -5},
      { -5.0,     5},
      { -5.0,    10},
      {  0.0,     5},
      {  0.0,    10},
      {  5.0,   -10},
      {  5.0,    -5},
      {  5.0,     5},
      {  5.0,    10},
      { 10.0,   -10},
      { 10.0,    -5},
      { 10.0,     5},
      { 10.0,    10}
    };  // clang-format on

    template<typename Number>
    static void func(Number* x, Number* y) {
      y[0] = x[0];
      y[0] /= x[1];  // R x (R \ {0})
      y[1] = 5.00;
      y[1] /= x[1];  // R x (R \ {0})
      y[2] = x[0];
      y[2] /= 5.00;  // R x (R \ {0})
    }
};
