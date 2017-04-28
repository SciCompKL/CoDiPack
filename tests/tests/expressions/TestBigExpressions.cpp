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

IN(12)
OUT(1)
POINTS(1) =
{
  {1.25, 2.5, 3.25, 4.5, 5.75, 6.25, 7.5, 8.5, 9.25, 10.25, 11.75, 12.5}

};

void func(NUMBER* x, NUMBER* y) {
  NUMBER* val_coord_Edge_CG = &x[0];
  NUMBER* val_coord_FaceElem_CG = &x[3];
  NUMBER* val_coord_Elem_CG = &x[6];
  NUMBER* val_coord_Point = &x[9];

  auto vec_a0 =  val_coord_Edge_CG[0]-val_coord_Point[0];
  auto vec_a1 =  val_coord_Edge_CG[1]-val_coord_Point[1];
  auto vec_a2 =  val_coord_Edge_CG[2]-val_coord_Point[2];

  auto vec_b0 =  val_coord_FaceElem_CG[0]-val_coord_Point[0];
  auto vec_b1 =  val_coord_FaceElem_CG[1]-val_coord_Point[1];
  auto vec_b2 =  val_coord_FaceElem_CG[2]-val_coord_Point[2];

  auto vec_c0 =  val_coord_Elem_CG[0]-val_coord_Point[0];
  auto vec_c1 =  val_coord_Elem_CG[1]-val_coord_Point[1];
  auto vec_c2 =  val_coord_Elem_CG[2]-val_coord_Point[2];

  auto vec_d0 = vec_a1*vec_b2-vec_a2*vec_b1;
  auto vec_d1 = -(vec_a0*vec_b2-vec_a2*vec_b0);
  auto vec_d2 = vec_a0*vec_b1-vec_a1*vec_b0;
//  NUMBER vec_a0 =  val_coord_Edge_CG[0]-val_coord_Point[0];
//  NUMBER vec_a1 =  val_coord_Edge_CG[1]-val_coord_Point[1];
//  NUMBER vec_a2 =  val_coord_Edge_CG[2]-val_coord_Point[2];
//
//  NUMBER vec_b0 =  val_coord_FaceElem_CG[0]-val_coord_Point[0];
//  NUMBER vec_b1 =  val_coord_FaceElem_CG[1]-val_coord_Point[1];
//  NUMBER vec_b2 =  val_coord_FaceElem_CG[2]-val_coord_Point[2];
//
//  NUMBER vec_c0 =  val_coord_Elem_CG[0]-val_coord_Point[0];
//  NUMBER vec_c1 =  val_coord_Elem_CG[1]-val_coord_Point[1];
//  NUMBER vec_c2 =  val_coord_Elem_CG[2]-val_coord_Point[2];
//
//  NUMBER vec_d0 = vec_a1*vec_b2-vec_a2*vec_b1;
//  NUMBER vec_d1 = -(vec_a0*vec_b2-vec_a2*vec_b0);
//  NUMBER vec_d2 = vec_a0*vec_b1-vec_a1*vec_b0;

  y[0] = fabs(vec_c0*vec_d0 + vec_c1*vec_d1 + vec_c2*vec_d2)/6.0;

}
