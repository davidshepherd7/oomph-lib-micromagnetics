// parameters for M = [cos(2*pi*x)*cos(t), 0, 0]

// Note that MxH ends up as zero - so not much is tested here...
// see 26/10/11 for more details

#include "math.h"
using namespace std;
using namespace MathematicalConstants;

namespace TwoDMicromagSetup
{

  using namespace TwoDMicromagSetup;

  /// Compute all components of the exact solution at time t, position x
  void exact_solution(const double& t, const Vector<double>& x,
		      Vector<double>& solution)
  {
    // phi
    solution[0] = 2*sin(2*Pi*x[0])*cos(t);

    // x,y,z components of M respectively
    solution[1] = cos(2*Pi*x[0])*cos(t);
    solution[2] = 0.0;
    solution[3] = 0.0;
  }

  void llg_source_function(const double& t, const Vector<double>& x,
			   Vector<double>& source)
  {
    // Source function to exactly cancel out contributions from not being a real exact solution
    source[0] = -cos(Pi*x[0]*2.0)*sin(t);
    source[1] = 0.0;
    source[2] = 0.0;
  }

}; // End of namespace
