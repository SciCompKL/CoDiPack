#include <toolDefines.h>

IN(2)
OUT(1)
POINTS(1) = {{1.0, 2.0}};

void func(NUMBER* x, NUMBER* y) {
  y[0] = x[0] + x[1];
}
