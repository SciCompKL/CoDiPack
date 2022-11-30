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

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "../../../config.h"
#include "../../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  namespace algorithms {

    struct FileIOInterface {
        using WriteHandle = void*;
        using ReadHandle = void*;

        WriteHandle openWrite(std::string const& filename, size_t totalSize);
        ReadHandle openRead(std::string const& filename);

        void closeWrite(WriteHandle handle);
        void closeRead(ReadHandle handle);

        std::string getFileEnding();

        template<typename T>
        void write(WriteHandle handle, T const* data, size_t elements);
        template<typename T>
        void write(WriteHandle handle, T const& data);

        template<typename T>
        void read(ReadHandle handle, T* data, size_t elements);
        template<typename T>
        void read(ReadHandle handle, T& data);
    };
  }
}
