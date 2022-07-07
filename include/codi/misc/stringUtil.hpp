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

#include <stdarg.h>
#include <string>

#include "../config.h"
#include "macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  namespace StringUtil {

    /**
     * String format function which uses printf formatting specifiers.
     *
     * @param format        The format specifier like printf
     * @param list          The variable argument list for the format string
     *
     * @return  The output with the formated values.
     */
    inline std::string vformat(const char* format, va_list list) {
        const int bufferSize = 200;
        char buffer[bufferSize];

        // copy the list if we need to iterate through the variables again
        va_list listCpy;
        va_copy(listCpy, list);


        int outSize = vsnprintf(buffer, bufferSize, format, list);

        std::string result;
        if(outSize + 1 > bufferSize) {
            char* newBuffer = new char[outSize + 1];

            outSize = vsnprintf(newBuffer, outSize + 1, format, listCpy);

            result = newBuffer;

            delete [] newBuffer;
        } else {
            result = buffer;
        }

        // cleanup the copied list
        va_end (listCpy);

        return result;
    }

    /**
     * String format function which uses printf formatting specifiers.
     *
     * @param format        The format specifier like printf
     * @param ...           The values for the format string
     *
     * @return  The output with the formated values.
     */
    inline std::string format(const char* format, ...) {
        va_list list;
        va_start(list, format);
        std::string output = vformat(format, list);
        va_end(list);

        return output;
    }
  }
}
