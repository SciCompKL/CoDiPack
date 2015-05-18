#include <toolDefines.h>

IN(2)
OUT(1)
POINTS(1) = {{2.0, 3.0}};

void func_forward(NUMBER& z, const NUMBER& w, const NUMBER& v){
  z = w*v;
}
#ifdef CHUNK_TAPE
static void ext_func(void* checkpoint){
  DataStore *check = static_cast<DataStore*>(checkpoint);
  NUMBER *x, w;
  check->getData(x);
  check->getData(w); 
  x[0].gradient() += x[1].getValue()*w.getGradient();
  x[1].gradient() += x[0].getValue()*w.getGradient(); 
}
void func(NUMBER* x, NUMBER* y) {
  DataStore *Checkpoint = new DataStore;
  NUMBER w;
  codi::RealReverse::globalTape.setPassive();
  func_forward(w,x[0],x[1]);
  codi::RealReverse::globalTape.setActive();

  codi::RealReverse::globalTape.registerInput(w);
  Checkpoint->addData(x);
  Checkpoint->addData(w);
  codi::RealReverse::globalTape.pushExternalFunction(&ext_func, Checkpoint, NULL);


  y[0] = w*w;
}
#else
void func(NUMBER*x, NUMBER* y){
  NUMBER w;
  func_forward(w, x[0], x[1]);
  y[0] = w*w;
}
#endif

