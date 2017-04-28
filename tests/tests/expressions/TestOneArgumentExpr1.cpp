/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2017 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#include <toolDefines.h>

IN(1)
OUT(15)
POINTS(41) =
{
  {-10.0000},
  { -9.5000},
  { -9.0000},
  { -8.5000},
  { -8.0000},
  { -7.5000},
  { -7.0000},
  { -6.5000},
  { -6.0000},
  { -5.5000},
  { -5.0000},
  { -4.5000},
  { -4.0000},
  { -3.5000},
  { -3.0000},
  { -2.5000},
  { -2.0000},
  { -1.5000},
  { -1.0000},
  { -0.5000},
  {  0.0000},
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
};

void func(NUMBER* x, NUMBER* y) {
  y[ 0] =      +x[0];   // R
  y[ 1] =      -x[0];   // R
  y[ 2] =   exp(x[0]);  // R
  y[ 3] =   sin(x[0]);  // R
  y[ 4] =   cos(x[0]);  // R
  y[ 5] =  atan(x[0]);  // R
  y[ 6] =  sinh(x[0]);  // R
  y[ 7] =  cosh(x[0]);  // R
  y[ 8] =   abs(x[0]);  // R
  y[ 9] =  fabs(x[0]);  // R
  y[10] = floor(x[0]);  // R
  y[11] =  tanh(x[0]);  // R
  y[12] =   tan(x[0]);  // R \ {(1/2 + i) * PI | for all i in Z}
  y[13] =   erf(x[0]);  // R
  y[14] =  erfc(x[0]);  // R
}
