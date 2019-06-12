/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2018 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

#pragma once

#include <codi.hpp>

  const size_t DefaultVectorDim = 5;

#if defined(FWD)
  typedef codi::RealForward NUMBER;

#elif defined(FWD2nd)
  #define SECOND_ORDER 1
  typedef codi::RealForwardGen<codi::RealForward> NUMBER;

#elif defined(FWD_Vec)
  typedef codi::RealForwardVec<DefaultVectorDim> NUMBER;

#elif defined(RWS_Chunk)
  typedef codi::RealReverse NUMBER;
  #define REVERSE_TAPE 1

#elif defined(RWS_ChunkVec)
  typedef codi::RealReverseVec<DefaultVectorDim> NUMBER;
  #define REVERSE_TAPE 1

#elif defined(RWS_Unch)
  typedef codi::RealReverseUnchecked NUMBER;
  #define REVERSE_TAPE 1

#elif defined(RWS_ChunkInd)
  typedef codi::RealReverseIndex NUMBER;
  #define REVERSE_TAPE 1

#elif defined(RWS_ChunkIndVec)
  typedef codi::RealReverseIndexVec<DefaultVectorDim> NUMBER;
  #define REVERSE_TAPE 1

#elif defined(RWS_UnchInd)
  typedef codi::RealReverseIndexUnchecked NUMBER;
  #define REVERSE_TAPE 1

#elif defined(RWS_Prim)
  typedef codi::RealReversePrimal NUMBER;
  #define REVERSE_TAPE 1

#elif defined(RWS_PrimVec)
  typedef codi::RealReversePrimalVec<DefaultVectorDim> NUMBER;
  #define REVERSE_TAPE 1

#elif defined(RWS_PrimIndex)
  typedef codi::RealReversePrimalIndex NUMBER;
  #define REVERSE_TAPE 1

#elif defined(RWS_PrimUnch)
  typedef codi::RealReversePrimalUnchecked NUMBER;
  #define REVERSE_TAPE 1

#elif defined(RWS2nd)
  #define SECOND_ORDER 1
  typedef codi::RealReverseGen<codi::RealForward> NUMBER;
  #define REVERSE_TAPE 1

#elif defined(RWS2nd_Prim)
  #define SECOND_ORDER 1
  typedef codi::RealReversePrimalGen<codi::RealForward> NUMBER;
  #define REVERSE_TAPE 1

#elif defined(DOUBLE)
  using namespace std;
  typedef double NUMBER;

#else
# error "No CoDi type defined"

#endif

#if !defined(DOUBLE)
  using Real = NUMBER::Real;
  using Gradient = NUMBER::GradientValue;
#else
  using Real = double;
  using Gradient = double;
#endif

#include "globalDefines.h"
