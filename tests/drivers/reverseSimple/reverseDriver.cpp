#include <toolDefines.h>

#include <iostream>
#include <vector>

int main(int nargs, char** args) {
  int evalPoints = getEvalPointsCount();
  int inputs = getInputCount();
  int outputs = getOutputCount();
  NUMBER* x = new NUMBER[inputs];
  NUMBER* y = new NUMBER[outputs];

  codi::SimpleTape<double, int>& tape = codi::RealReverseSimple::globalTape;
  tape.resize(1000, 1000);

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

    std::vector<std::vector<double> > jac(outputs);
    for(int curOut = 0; curOut < outputs; ++curOut) {
      for(int i = 0; i < inputs; ++i) {
        tape.registerInput(x[i]);
      }

      func(x, y);

      y[curOut].setGradient(1.0);
      tape.evaluate();

      for(int curIn = 0; curIn < inputs; ++curIn) {
        jac[curOut].push_back(x[curIn].getGradient());
      }

      tape.reset();
    }

    for(int curIn = 0; curIn < inputs; ++curIn) {
      for(int curOut = 0; curOut < outputs; ++curOut) {
        std::cout << curIn << " " << curOut << " " << jac[curOut][curIn] << std::endl;
      }
    }
  }
}
