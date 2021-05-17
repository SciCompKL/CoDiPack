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

  using Factory = codi::AggregatedTypeVectorAccessWrapperFactory<Type>;
  using VectorWrapper = typename Factory::RType;                    // Step 2: Create a wrapped vector access interface.

  VectorWrapper* vaType = Factory::create(va);

  using TypeIdentifier = typename VectorWrapper::Identifier;
  using TypeReal = typename VectorWrapper::Real;

  TypeIdentifier t_i = data->getData<TypeIdentifier>();

  TypeReal t_b = vaType->getAdjoint(t_i, 0);                        // Step 3: Use the wrapped vector access interface.

  std::cout << " Reverse: t_b = " << t_b  << std::endl;

  Factory::destroy(vaType);                                         // Step 4: Delete the created wrapper.
}

void extFunc_del(Tape* t, void* d) {
  codi::ExternalFunctionUserData* data = (codi::ExternalFunctionUserData*)d;

  delete data;

  std::cout << " Reset: data is deleted." <<  std::endl;
}

template<typename Type>
void addExternalFunc(Type const& v) {

  Tape& tape = Real::getGlobalTape();

  codi::ExternalFunctionUserData* data = new codi::ExternalFunctionUserData();
  data->addData(codi::RealTraits::getIdentifier(v));                // Step 1: Use the general access routines an values of the type.

  tape.pushExternalFunction(codi::ExternalFunction<Tape>::create(extFunc_rev<Type>, data, extFunc_del));
}
//! [Typed external function]

int main(int nargs, char** args) {
  Real x = 4.0;

  Tape& tape = Real::getGlobalTape();
  tape.setActive();

  tape.registerInput(x);
  Real t = func(x);

  addExternalFunc(t);

  std::complex<Real> c(t, -t);
  std::complex<Real> w = func(c);

  addExternalFunc(w);

  Real y = func(std::abs(w));
  tape.registerOutput(y);

  tape.setPassive();
  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "f(4.0) = " << y << std::endl;
  std::cout << "df/dx(4.0) = " << x.getGradient() << std::endl;

  tape.reset();

  return 0;
}
//! [Example 20 - Aggregated active type handling]
