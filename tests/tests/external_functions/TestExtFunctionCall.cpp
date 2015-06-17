#include <toolDefines.h>

#include <tools/dataStore.hpp>

#include <iostream>

IN(2)
OUT(1)
POINTS(1) = {{2.0, 3.0}};

void func_forward(NUMBER& z, const NUMBER& w, const NUMBER& v){
  z = w*v;
}

#ifdef REVERSE_TAPE
static void extFunc(void* checkpoint){
  codi::DataStore *check = static_cast<codi::DataStore*>(checkpoint);
  NUMBER *x, w;
  check->getData(x);
  check->getData(w);
  x[0].gradient() += x[1].getValue()*w.getGradient();
  x[1].gradient() += x[0].getValue()*w.getGradient();
}

static void delFunc(void* checkpoint){
  codi::DataStore *check = static_cast<codi::DataStore*>(checkpoint);
  delete check;
  std::cout << "Delete" << std::endl;
}

void func(NUMBER* x, NUMBER* y) {
  NUMBER::TapeType& tape = NUMBER::getGlobalTape();
  codi::DataStore *checkpoint = new codi::DataStore;
  NUMBER w;
  tape.setPassive();
  func_forward(w,x[0],x[1]);
  tape.setActive();

  tape.registerInput(w);
  checkpoint->addData(x);
  checkpoint->addData(w);
  tape.pushExternalFunctionHandle(&extFunc, checkpoint, delFunc);


  y[0] = w*w;
}
#else
void func(NUMBER*x, NUMBER* y){
  NUMBER w;
  func_forward(w, x[0], x[1]);
  y[0] = w*w;
}
#endif

