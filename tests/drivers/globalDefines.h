/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

#define POINTS(number) \
  int getEvalPointsCount() {return number;} \
  extern double points[number][in_count]; \
  double getEvalPoint(int point, int col) { return points[point][col]; } \
  double points[number][in_count]

#define IN(number) \
  const int in_count = number; \
  int getInputCount() {return number;}

#define OUT(number) \
  int getOutputCount() {return number;}

#ifndef SECOND_ORDER
# define SECOND_ORDER 0
#endif

#ifndef PRIMAL
# define PRIMAL 0
#endif

#ifndef REVERSE_TAPE
# define REVERSE_TAPE 0
#endif

#if !defined(DOUBLE)
  using Real = NUMBER::Real;
  using Gradient = NUMBER::GradientValue;
#else
  using Real = double;
  using Gradient = double;
#endif

int getEvalPointsCount();
double getEvalPoint(int point, int col);
int getInputCount();
int getOutputCount();

void func(NUMBER* x, NUMBER* y);
