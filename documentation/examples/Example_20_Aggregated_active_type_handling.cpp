//! [Example 20 - Aggregated active type handling]

#include <codi.hpp>
#include <iostream>

using Real = codi::RealReverse;
using Tape = typename Real::Tape;
using Identifier = typename Real::Identifier;
using RealBase = typename Real::Real;

//! [Function]
template<typename Type>
Type func(const Type& x) {
  return x * x;
}
//! [Function]


//! [Typed external function]
template<typename Type>
void extFunc_rev(Tape* t, void* d, codi::VectorAccessInterface<RealBase, Identifier>* va) {
  codi::ExternalFunctionUserData* data = (codi::ExternalFunctionUserData*)d;

  // Step 3: Create a wrapped vector access interface.
  using Factory = codi::AggregatedTypeVectorAccessWrapperFactory<Type>;
  using VectorWrapper = typename Factory::RType;
  VectorWrapper* vaType = Factory::create(va);

  using TypeIdentifier = typename VectorWrapper::Identifier;
  using TypeReal = typename VectorWrapper::Real;

  // Step 4: Get the external function data
  TypeReal x_v = data->getData<TypeReal>();
  TypeIdentifier x_i = data->getData<TypeIdentifier>();
  TypeIdentifier w_i = data->getData<TypeIdentifier>();

  // Step 5: Use the wrapped vector access interface and perform the adjoint operation
  TypeReal w_b = vaType->getAdjoint(w_i, 0);
  TypeReal t_b = 2.0 * codi::ComputationTraits::transpose(x_v) * w_b;

  vaType->updateAdjoint(x_i, 0, t_b);
  vaType->resetAdjoint(w_i, 0);

  // Step 6: Delete the created wrapper.
  Factory::destroy(vaType);
}

void extFunc_del(Tape* t, void* d) {
  codi::ExternalFunctionUserData* data = (codi::ExternalFunctionUserData*)d;

  delete data;

  std::cout << " Reset: data is deleted." <<  std::endl;
}

template<typename Type>
Type addExternalFunc(Type const& x) {

  Tape& tape = Real::getTape();

  // Step 1: Perform the passive function evaluation.
  tape.setPassive();
  Type w = func(x);
  tape.setActive();
  codi::RealTraits::registerExternalFunctionOutput(w);

  // Step 2: Use the general access routines on the values to extract the primal and identifier data.
  codi::ExternalFunctionUserData* data = new codi::ExternalFunctionUserData();
  data->addData(codi::RealTraits::getValue(x));
  data->addData(codi::RealTraits::getIdentifier(x));
  data->addData(codi::RealTraits::getIdentifier(w));

  tape.pushExternalFunction(codi::ExternalFunction<Tape>::create(extFunc_rev<Type>, data, extFunc_del));

  return w;
}
//! [Typed external function]

int main(int nargs, char** args) {
  Real x = 3.0;

  Tape& tape = Real::getTape();
  tape.setActive();

  tape.registerInput(x);
  Real t1 = addExternalFunc(x);

  std::complex<Real> c(t1, -t1);
  std::complex<Real> t2 = addExternalFunc(c);

  Real y = std::abs(t2);
  tape.registerOutput(y);

  tape.setPassive();
  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "x = " << x << std::endl;
  std::cout << "y = " << y << std::endl;
  std::cout << "dy/dx = " << x.getGradient() << std::endl;

  tape.reset();

  return 0;
}
//! [Example 20 - Aggregated active type handling]
