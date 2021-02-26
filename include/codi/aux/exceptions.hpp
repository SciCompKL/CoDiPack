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
   * The method prints an error message that contains the failed expression, the function name, the file and the line.
   *
   * @param[in]       condition  The evaluated value of the condition.
   * @param[in] conditionString  The condition as a string value.
   * @param[in]        function  The name of the function that caused the assert.
   * @param[in]            file  The file were the function is defined.
   * @param[in]            line  The line in the file were the assert is defined.
   */
  inline void checkAndOutputAssert(bool const condition, char const* conditionString, char const* function, char const* file, int line) {
    if (!condition) {
      std::cerr << "codiAssertion failed: " << conditionString << " in function " << function << " at " << file << ":" << line << std::endl;
      abort();
    }
  }

  /**
   * @brief Generates an exception.
   *
   * @param...  Arguments for a printf like output and format.
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
  inline void outputException(char const function[], char const file[], int const line, char const* message,...) {
    fprintf(stderr, "Error in function %s (%s:%d)\nThe message is: ", function, file, line);

    va_list vl;
    va_start(vl, message);
    vfprintf(stderr, message, vl);
    va_end(vl);

    fprintf(stderr, "\n");
    exit(-1);
  }
}

