#include <toolDefines.h>

#include <iostream>

int main(int nargs, char** args) {
  int evalPoints = getEvalPointsCount();
  int inputs = getInputCount();
  int outputs = getOutputCount();
  NUMBER* x = new NUMBER[inputs];
  NUMBER* y = new NUMBER[outputs];

  for(int curPoint = 0; curPoint < evalPoints; ++curPoint) {
    std::cout << "Point " << curPoint << " : {";

    for(int i = 0; i < inputs; ++i) {
      if(i != 0) {
        std::cout << ", ";
      }
      double val = getEvalPoint(curPoint, i);
      std::cout << val;

      x[i] = (NUMBER)(val);
    }
    std::cout << "}\n";

    for(int i = 0; i < outputs; ++i) {
      y[i] = 0.0;
    }

    for(int curIn = 0; curIn < inputs; ++curIn) {
      x[curIn].set_gradient(1.0);
      for(int i = 0; i < outputs; ++i) {
        y[i].set_gradient(0.0);
      }

      func(x, y);

      for(int curOut = 0; curOut < outputs; ++curOut) {
        std::cout << curIn << " " << curOut << " " << y[curOut].get_gradient() << std::endl;
      }

      x[curIn].set_gradient(0.0);
    }
  }
}
