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

struct TestOneArgumentExpr2 : public TestInterface {
  public:
    NAME("OneArgumentExpr2")
    IN(1)
    OUT(4)
    POINTS(20) =  // clang-format off
    {
      {  0.5000},
      {  1.0000},
      {  1.5000},
      {  2.0000},
      {  2.5000},
      {  3.0000},
      {  3.5000},
      {  4.0000},
      {  4.5000},
      {  5.0000},
      {  5.5000},
      {  6.0000},
      {  6.5000},
      {  7.0000},
      {  7.5000},
      {  8.0000},
      {  8.5000},
      {  9.0000},
      {  9.5000},
      { 10.0000}
    };  // clang-format on

    template<typename Number>
    static void func(Number* x, Number* y) {
      y[0] = log(x[0]);     // (0, inf)
      y[1] = log10(x[0]);   // (0, inf)
      y[2] = sqrt(x[0]);    // [0, inf)
      y[3] = tgamma(x[0]);  // R currently only defined for positive arguments
    }
};
