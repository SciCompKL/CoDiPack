/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Checks the assert statement and aborts the program if the statement is false.
   *
   * The method prints an error message that contains the failed expression, the function name, the file, and the line.
   * This method is usually accessed via the codiAssert macro, see codi::Config::EnableAssert.
   *
   * @param[in]       condition  The evaluated value of the condition.
   * @param[in] conditionString  The condition as a string value.
   * @param[in]        function  The name of the function that caused the assert.
   * @param[in]            file  The file were the function is defined.
   * @param[in]            line  The line in the file were the assert is defined.
   */
  inline void checkAndOutputAssert(bool const condition, char const* conditionString, char const* function,
                                   char const* file, int line) {
    if (!condition) {
      std::cerr << "codiAssertion failed: " << conditionString << " in function " << function << " at " << file << ":"
                << line << std::endl;
      abort();
    }
  }

/**
 * @brief Generates an exception.
 *
 * @param ...  Arguments for a printf like output and format.
 */
#define CODI_EXCEPTION(...) outputException(__func__, __FILE__, __LINE__, __VA_ARGS__)

  /**
   * @brief Prints the positions and the message of the exception.
   *
   * The position and function where the exceptions occurred are printed. The message will be handled
   * and formatted by a printf like function.
   *
   * @param[in] function  Name of the function from which the exception is generated.
   * @param[in]     file  File where the exception was generated.
   * @param[in]     line  Line inside the file where the exception was generated.
   * @param[in]  message  The exception message and the arguments for the formatting in the message.
   */
  inline void outputException(char const function[], char const file[], int const line, char const* message, ...) {
    fprintf(stderr, "Error in function %s (%s:%d)\nThe message is: ", function, file, line);

    va_list vl;
    va_start(vl, message);
    vfprintf(stderr, message, vl);
    va_end(vl);

    fprintf(stderr, "\n");
    exit(-1);
  }

#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
  #define DEPRECATE(foo, msg) foo __attribute__((deprecated(msg)))
#elif defined(_MSC_VER)
  #define DEPRECATE(foo, msg) __declspec(deprecated(msg)) foo
#else
  #error This compiler is not supported.
#endif

  /// Presents presenting compile time warnings to the user.
  /// The warning is presented as a deprecated note.
  struct Warning {
    public:
      /// Show a warning about an implicit cast of an active real type.
      template<bool v>
      static void implicitCast() {
        implicitCastStatic(::std::integral_constant<bool, v>());
      }

      /// Implementation that displayes the warning.
      DEPRECATE(static void implicitCastStatic(::std::false_type const&),
                "static_warning: Implicit conversion of CoDiPack type to real.") {}

      /// Implementation that ignores the warning.
      static void implicitCastStatic(::std::true_type const&) {}
  };
}
