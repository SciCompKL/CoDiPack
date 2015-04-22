#include <toolDefines.h>

#include <iostream>

int main(int nargs, char** args) {
  int evalPoints = getEvalPointsCount();
  int inputs = getInputCount();
  int outputs = getOutputCount();
  NUMBER* x = new NUMBER[inputs];
  NUMBER* y = new NUMBER[outputs];

  adept::Stack tape(true);
  tape.pause_recording();

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
      tape.continue_recording();
      for(int i = 0; i < inputs; ++i) {
        x[i].register_gradient();
      }

      func(x, y);
      tape.pause_recording();

      y[curOut].set_gradient(1.0);
      tape.compute_adjoint();

      for(int curIn = 0; curIn < inputs; ++curIn) {
        jac[curOut].push_back(x[curIn].get_gradient());
      }

      tape.clear_gradients();
    }

    for(int curIn = 0; curIn < inputs; ++curIn) {
      for(int curOut = 0; curOut < outputs; ++curOut) {
        std::cout << curIn << " " << curOut << " " << jac[curOut][curIn] << std::endl;
      }
    }
  }
}
