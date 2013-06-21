
// Floating point debugging
#include <fenv.h>

#include "../../micromag.h"

using namespace oomph;


int main(int argc, char *argv[])
{
  // Start MPI if necessary
  MPI_Helpers::init(argc,argv);

  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);




  // Shut down oomph-lib's MPI
  MPI_Helpers::finalize();

  return 0;
}
