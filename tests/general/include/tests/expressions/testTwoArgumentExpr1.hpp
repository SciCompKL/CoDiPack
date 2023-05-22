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
#include "../../testInterface.hpp"

struct TestTwoArgumentExpr1 : public TestInterface {
  public:
    NAME("TwoArgumentExpr1")
    IN(2)
    OUT(32)
    POINTS(25) =  // clang-format off
    {
      {-10.0,   -10},
      {-10.0,    -5},
      {-10.0,     0},
      {-10.0,     5},
      {-10.0,    10},
      { -5.0,   -10},
      { -5.0,    -5},
      { -5.0,     0},
      { -5.0,     5},
      { -5.0,    10},
      {  0.0,   -10},
      {  0.0,    -5},
      {  0.0,     0},
      {  0.0,     5},
      {  0.0,    10},
      {  5.0,   -10},
      {  5.0,    -5},
      {  5.0,     0},
      {  5.0,     5},
      {  5.0,    10},
      { 10.0,   -10},
      { 10.0,    -5},
      { 10.0,     0},
      { 10.0,     5},
      { 10.0,    10}
    };  // clang-format on

    template<typename Number>
    static void func(Number* x, Number* y) {
      y[0] = x[0] + x[1];        // R x R
      y[1] = 5.00 + x[1];        // R x R
      y[2] = x[0] + 5.00;        // R x R
      y[3] = x[0] - x[1];        // R x R
      y[4] = 5.00 - x[1];        // R x R
      y[5] = x[0] - 5.00;        // R x R
      y[6] = x[0] * x[1];        // R x R
      y[7] = 5.00 * x[1];        // R x R
      y[8] = x[0] * 5.00;        // R x R
      y[9] = min(x[0], x[1]);    // R x R
      y[10] = min(5.00, x[1]);   // R x R
      y[11] = min(x[0], 5.00);   // R x R
      y[12] = max(x[0], x[1]);   // R x R
      y[13] = max(5.00, x[1]);   // R x R
      y[14] = max(x[0], 5.00);   // R x R
      y[15] = fmin(x[0], x[1]);  // R x R the results for the points with the same values are reversed here because the
                                 // other one uses the standard template.
      y[16] = fmin(5.00, x[1]);  // R x R
      y[17] = fmin(x[0], 5.00);  // R x R
      y[18] = fmax(x[0], x[1]);  // R x R the results for the points with the same values are reversed here because the
                                 // other one uses the standard template.
      y[19] = fmax(5.00, x[1]);  // R x R
      y[20] = fmax(x[0], 5.00);  // R x R
      y[21] = pow(x[0], x[1]);   // R x R
      y[22] = pow(5.00, x[1]);   // R x R
      y[23] = pow(x[0], 5.00);   // R x R
      y[24] = copysign(x[0], x[1]);   // R x R
      y[25] = copysign(5.00, x[1]);   // R x R
      y[26] = copysign(x[0], 5.00);   // R x R
      y[27] = copysignf(x[0], x[1]);  // R x R
      y[28] = copysignf(5.00, x[1]);  // R x R
      y[29] = copysignf(x[0], 5.00);  // R x R

      y[30] = x[0];
      y[31] = x[1];
      swap(y[30], y[31]);  // R x R
    }
};
