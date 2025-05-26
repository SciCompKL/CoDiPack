/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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

#include <codi/expressions/complex/stdComplex.hpp>
#include <complex>

template<typename C, typename R>
static void assignToComplex(C* c, R const* r, int complexCount) {
  for (int i = 0; i < complexCount; i += 1) {
    c[i] = C(r[i * 2], r[i * 2 + 1]);
  }
}

template<typename R, typename C>
static void assignToReal(R* r, C const* c, int complexCount) {
  for (int i = 0; i < complexCount; i += 1) {
    r[i * 2] = std::real(c[i]);
    r[i * 2 + 1] = std::imag(c[i]);
  }
}

template<typename Arg>
static codi::RealTraits::PassiveReal<Arg> passive(Arg const& arg) {
  return codi::RealTraits::getPassiveValue(arg);
}

template<typename Arg>
static std::complex<codi::RealTraits::PassiveReal<Arg>> passive(std::complex<Arg> const& arg) {
  return {codi::RealTraits::getPassiveValue(std::real(arg)), codi::RealTraits::getPassiveValue(std::imag(arg))};
}

template<typename Arg>
static std::complex<codi::RealTraits::PassiveReal<Arg>> passive(codi::ActiveComplex<Arg> const& arg) {
  return {codi::RealTraits::getPassiveValue(std::real(arg)), codi::RealTraits::getPassiveValue(std::imag(arg))};
}

#if CODI_SpecializeStdComplex
template<typename T>
using TestComplex = std::complex<T>;
#else
template<typename T>
using TestComplex = codi::ActiveComplex<T>;
#endif
