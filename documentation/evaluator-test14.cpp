//
// Test codi::Evaluator with c++14 generic lambdas
//

#include <iostream>
#include <cstdio>
#include <codi.hpp>



int main() 
{

  auto myfunc=[](const auto  &x, auto &y)
    {
      y[0]=sin(x[0]);
    };
  
  auto evaluator= codi::Evaluator<decltype(myfunc)>(1,1,myfunc);

  std::vector<double> xx(1);
  std::vector<double> yy(1);

  for (double x=0.0; x<10.0; x+=1.0)
  {
    xx[0]=x;
    evaluator.call(xx);
    printf("f(x)=%g f'(x)=%g cos(x)=%g\n", evaluator.result(0),evaluator.jacobian(0,0),cos(x));
  }


  auto myfunc2=[](const auto &x, auto &y)
  {
    y[0]=x[0]*11+x[1]*12;
    y[1]=x[0]*21+x[1]*22;
  };

  auto evaluator2=codi::Evaluator<decltype(myfunc2)>(2,2,myfunc2);
  std::vector<double> xx2(2);
  std::vector<double> yy2(2);

  
  xx2[0]=1;
  xx2[1]=1;
  evaluator2.call(xx2);
  printf("f:  %g %g\n",evaluator2.result(0), evaluator2.result(1));
  printf("f': %g %g\n",evaluator2.jacobian(0,0),evaluator2.jacobian(0,1));
  printf("    %g %g\n",evaluator2.jacobian(1,0),evaluator2.jacobian(1,1));
  return 0;
}
