
// Floating point debugging
#include <fenv.h>

#include "../poisson_test_problem.h"

using namespace oomph;

// ============================================================
// Generic poisson test
// ============================================================
int generic_poisson_test()
{
 GenericPoissonForTests test_poisson;
 return test_poisson.run_test();
}


int main()
{
 // Enable some floating point error checkers
 feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

 int result = generic_poisson_test();

 if(result == 0)
  {
   oomph_info << "***" <<std::endl;
   oomph_info << "*** generic_poisson_test() passed" << std::endl;
   oomph_info << "***" <<std::endl;
  }

 return result;
}
