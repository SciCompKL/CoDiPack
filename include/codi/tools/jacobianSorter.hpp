/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
#include <array>

#include "../configure.h"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Sorts Jacobian entries on Jacobian tapes.
   *
   * The sorter buffers the pushed entries for each statement and
   * add Jacobian values for arguments that have the same identifier.
   *
   * @tparam         Real  The floating point value used by the tapes.
   * @tparam GradientData  The type of the Jacobian values.
   */
  template<typename Real, typename GradientData>
  struct JacobianSorter {

      std::array<GradientData, MaxStatementIntSize> indices; /**< Array of the identifiers */
      std::array<Real, MaxStatementIntSize> jacobies; /**< Array of the Jacobian values */
      StatementInt size; /**< Current number of arguments for the expression */

      JacobianSorter() = default; /**< Default constructor */

      /**
       * @brief Wrapper method that buffers the arguments for the statement.
       *
       * The method checks first if the identifier is already in the buffer.
       * If not not a new entry is created.
       *
       * @param[in] jacobi  The Jacobian value for the argument.
       * @param[in]  index  The identifier for the argument.
       */
      CODI_INLINE void setDataAndMove(const Real& jacobi, const GradientData& index) {
        bool found = false;
        StatementInt pos;
        for(pos = 0; pos < size; pos += 1) {
          if(indices[pos] == index) {
            found = true;
            break;
          }
        }

        if(!found) {
          size += 1;
          indices[pos] = index;
          jacobies[pos] = jacobi;
        } else {
          jacobies[pos] += jacobi;
        }

      }

      /**
       * @brief Adds the buffered data to the vector.
       *
       * @param[in,out] vec  The vector on which the data is added.
       *
       * @tparam Vec  Requires a 'setDataAndMove' method as in the codi::ChunkVector interface.
       */
      template<typename Vec>
      CODI_INLINE void storeData(Vec& vec) {
        for(StatementInt pos = 0; pos < size; pos += 1) {
          vec.setDataAndMove(jacobies[pos], indices[pos]);
        }

        // Reset the data for the next statement
        size = 0;
      }
  };
}
