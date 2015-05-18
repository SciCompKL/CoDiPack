#include <toolDefines.h>

#include <tools/DataStore.hpp>

#include <iostream>

IN(2)
OUT(1)
POINTS(1) = {{2.0, 3.0}};

void func_forward(NUMBER& z, const NUMBER& w, const NUMBER& v){
  z = w*v;
}

const int ITER = 5;

#ifdef REVERSE_TAPE
static void extFunc(void* checkpoint){
  DataStore *check = static_cast<DataStore*>(checkpoint);
  NUMBER *x, w0, w1;
  check->getData(x);
  check->getData(w0);
  check->getData(w1);
  w0.gradient() += x[1].getValue()*w1.getGradient();
  x[1].gradient() += w0.getValue()*w1.getGradient();
}

static void delFunc(void* checkpoint){
  DataStore *check = static_cast<DataStore*>(checkpoint);
  delete check;
  std::cout << "Delete" << std::endl;
}

void func(NUMBER* x, NUMBER* y) {
  NUMBER::TapeType& tape = NUMBER::globalTape;
  NUMBER w[ITER];

  w[0] = x[0];
  for(int i = 1; i < ITER; ++i) {
    tape.setPassive();
    func_forward(w[i],w[i - 1],x[1]);
    tape.setActive();

    DataStore *checkpoint = new DataStore();
    tape.registerInput(w[i]);
    checkpoint->addData(x);
    checkpoint->addData(w[i-1]);
    checkpoint->addData(w[i]);
    tape.pushExternalFunctionHandle(&extFunc, checkpoint, delFunc);
  }


  y[0] = w[ITER - 1]*w[ITER - 1];
}
#else
void func(NUMBER*x, NUMBER* y){
  NUMBER w[ITER];

  w[0] = x[0];
  for(int i = 1; i < ITER; ++i) {
    func_forward(w[i],w[i - 1],x[1]);
  }


  y[0] = w[ITER - 1]*w[ITER - 1];
}
#endif

