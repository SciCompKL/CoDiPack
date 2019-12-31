/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2020 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
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
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *     Max Sagebaum
 *     Tim Albring
 *     Johannes Blühdorn
 */

#include <toolDefines.h>

IN(2)
OUT(9)
POINTS(25) =
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
};

void func(NUMBER* x, NUMBER* y) {
  y[0] = x[0];
  y[0] += x[1];  // R x R
  y[1] = 5.00;
  y[1] += x[1];  // R x R
  y[2] = x[0];
  y[2] += 5.00;  // R x R
  y[3] = x[0];
  y[3] -= x[1];  // R x R
  y[4] = 5.00;
  y[4] -= x[1];  // R x R
  y[5] = x[0];
  y[5] -= 5.00;  // R x R
  y[6] = x[0];
  y[6] *= x[1];  // R x R
  y[7] = 5.00;
  y[7] *= x[1];  // R x R
  y[8] = x[0];
  y[8] *= 5.00;  // R x R
}
