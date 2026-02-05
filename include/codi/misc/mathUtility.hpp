/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), RPTU University Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, RPTU University Kaiserslautern-Landau)
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
 *  - SciComp, RPTU University Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */

#pragma once

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Helper function for overallocation in multiples of a given chunk size.
   * @tparam IntegralType   An integral type, usually size_t.
   * @param targetSize      Minimum returned size.
   * @param chunkSize       Returned size is a multiple of chunkSize.
   * @return Next multiple of chunkSize larger or equal to targetSize.
   */
  template<typename IntegralType>
  IntegralType getNextMultiple(IntegralType const& targetSize, IntegralType const& chunkSize) {
    IntegralType chunkCount = (targetSize + chunkSize - 1) / chunkSize;
    return chunkCount * chunkSize;
  }
}
