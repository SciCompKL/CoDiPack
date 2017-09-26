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
OUT(3)
POINTS(39) =
{
  { -0.9500},
  { -0.9000},
  { -0.8500},
  { -0.8000},
  { -0.7500},
  { -0.7000},
  { -0.6500},
  { -0.6000},
  { -0.5500},
  { -0.5000},
  { -0.4500},
  { -0.4000},
  { -0.3500},
  { -0.3000},
  { -0.2500},
  { -0.2000},
  { -0.1500},
  { -0.1000},
  { -0.0500},
  {  0.0000},
  {  0.0500},
  {  0.1000},
  {  0.1500},
  {  0.2000},
  {  0.2500},
  {  0.3000},
  {  0.3500},
  {  0.4000},
  {  0.4500},
  {  0.5000},
  {  0.5500},
  {  0.6000},
  {  0.6500},
  {  0.7000},
  {  0.7500},
  {  0.8000},
  {  0.8500},
  {  0.9000},
  {  0.9500}
};

void func(NUMBER* x, NUMBER* y) {
  y[0] = atanh(x[0]);  // (-1, 1)
  y[1] =  asin(x[0]);  // [-1, 1]
  y[2] =  acos(x[0]);  // [-1, 1]
}
