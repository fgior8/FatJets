#include "OtherFunctions.h"

double generate_test_func()
{
  MTRand_int32 irand(123);

  //max of function
  double max=2.36866988642;
  
  double gen_x, gen_fx;

  do {
    gen_x=irand.rand_double();
    gen_fx=irand.rand_double()*max;
  }  while(test_func(gen_x) < gen_fx);

  return gen_x;
}


double test_func(double x)
{
  //function is (1-x)^2 * x^4 / (2-x)
  // x from 0 to 1, normalized to one
  double normalization = 142.418582618187;
  
  return normalization*(1-x)*(1-x)*x*x*x*x/(2-x);
}
