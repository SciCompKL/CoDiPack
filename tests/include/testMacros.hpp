#pragma once

#define IN(number)                        \
  static int constexpr in_count = number; \
  int getInputCount() {                   \
    return number;                        \
  }

#define NAME(name)        \
  std::string getName() { \
    return name;          \
  }

#define OUT(number)                        \
  static int constexpr out_count = number; \
  int getOutputCount() {                   \
    return number;                         \
  }

#define POINTS(number)                      \
  int getEvalPointsCount() {                \
    return number;                          \
  }                                         \
  double getEvalPoint(int point, int col) { \
    return points[point][col];              \
  }                                         \
  double points[number][in_count]
