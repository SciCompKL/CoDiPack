/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 *  - SciComp, TU Kaiserslautern:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <cstddef>

#include "../../config.h"
#include "../../misc/fileIo.hpp"
#include "../../misc/macros.hpp"
#include "chunk.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @ brief Chunk with one entry per item.
   *
   * @tparam Chunk1  Any type.
   */
  struct ByteChunk final : public Chunk1<char> {
    public:

      using Base = Chunk1;  ///< Abbreviation for the base class type.

    public:

      // Use constructors.
      using Base::Base;

      /// \copydoc ChunkBase::pushData
      template<typename Type>
      CODI_INLINE void pushData(Type const& value1) {

        codiAssert(getUnusedSize() != 0);

        Type* typedData = (Type*)(&data1[usedSize]);
        typedData[0] = value1;
        usedSize += sizeof(Type);
      }
  };
}
