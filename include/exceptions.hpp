/**
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing, TU Kaiserslautern
 *
 * This file is part of CoDiPack.
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
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
 * Authors: TODO
 */
#pragma once

#include <cstdio>
#include <cstdarg>


namespace codi {
  /**
   * @brief Generates an exception.
   *
   * @param ...  Arguments for a printf like output and format.
   */
  #define CODI_EXCEPTION(...) outputException( __func__, __FILE__, __LINE__, __VA_ARGS__)

  /**
   * @brief Prints the positions and the message of the exception.
   *
   * The position and function where the exceptions occurred is printed. The message will be handled
   * and formatted by a printf like function.
   *
   * @param[in] function  Name of the function from which the exception is generated.
   * @param[in]     file  File where the exception was generated.
   * @param[in]     line  Line inside the file where the exception was generated.
   * @param[in]  message  The exception message and the arguments for the formatting in the message.
   */
  inline void outputException(const char function[], const char file[], const int line, const char* message, ...) {
    fprintf(stderr, "Error in function %s (%s:%d)\nThe message is: ", function, file, line);

    va_list vl;
    va_start(vl, message);
    vfprintf(stderr, message, vl);
    va_end(vl);

    fprintf(stderr, "\n");
    exit(-1);
  }
}

