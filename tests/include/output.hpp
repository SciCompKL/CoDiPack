#pragma once

#include <stdarg.h>
#include <stdio.h>
#include <string>
#include <vector>

#include "../../include/codi/tools/data/jacobian.hpp"
//#include ../../include/codi/tools/data/hessian.hpp"

char const *const HEADER_FORMAT = "%6s_%03zd";
char const *const VALUE_FORMAT  = "%10g";
char const *const COL_SEPERATOR = " ";
char const *const LINE_END      = "\n";
char const *const BLANK         = "          ";


inline std::string vformat(char const* format, va_list list) {
  int const bufferSize = 200;
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

inline std::string format(char const* format, ...) {
  va_list list;
  va_start(list, format);
  std::string result = vformat(format, list);
  va_end(list);

  return result;
}

template<typename T>
void writeOutputPrimal(FILE* out, std::vector<T> const& primal) {

  for(size_t curOut = 0; curOut < primal.size(); curOut += 1) {
    fprintf(out, HEADER_FORMAT, "out", curOut);
    fprintf(out, COL_SEPERATOR);
    fprintf(out, VALUE_FORMAT, primal[curOut]);
    fprintf(out, LINE_END);
  }
}

template<typename Vec>
void writeOutputJacobian(FILE* out, codi::Jacobian<Vec> const& jac) {

  // print header
  fprintf(out, BLANK);
  for(size_t curIn = 0; curIn < jac.getN(); curIn += 1) {
    fprintf(out, COL_SEPERATOR);
    fprintf(out, HEADER_FORMAT, "in", curIn);
  }
  fprintf(out, LINE_END);

  for(size_t curOut = 0; curOut < jac.getM(); curOut += 1) {
    fprintf(out, HEADER_FORMAT, "out", curOut);
    for(size_t curIn1st = 0; curIn1st < jac.getN(); curIn1st += 1) {
      fprintf(out, COL_SEPERATOR);
      fprintf(out, VALUE_FORMAT, jac(curOut, curIn1st));
    }

    fprintf(out, LINE_END);
  }
}

//template<typename Vec>
//void writeOutputHessian(codi::Hessian<Vec> const& hes) {

//  for(size_t curOut = 0; curOut < hes.getM(); curOut += 1) {

//    // print header
//    fprintf(out, HEADER_FORMAT, "out", curOut);
//    for(size_t curIn = 0; curIn < hes.getN(); curIn += 1) {
//      fprintf(out, COL_SEPERATOR);
//      fprintf(out, HEADER_FORMAT, "in", curIn);
//    }
//    fprintf(out, LINE_END);

//    for(size_t curIn1st = 0; curIn1st < hes.getN(); curIn1st += 1) {
//      fprintf(out, HEADER_FORMAT, "in", curIn1st);

//      for(size_t curIn2nd = 0; curIn2nd < hes.getN(); curIn2nd += 1) {
//        fprintf(out, COL_SEPERATOR);
//        fprintf(out, VALUE_FORMAT, hes(curOut, curIn1st, curIn2nd));
//      }

//      fprintf(out, LINE_END);
//    }

//    fprintf(out, LINE_END);
//  }
//}
