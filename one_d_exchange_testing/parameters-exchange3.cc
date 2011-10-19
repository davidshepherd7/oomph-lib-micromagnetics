// parameters for M = cos(t) [cos(2*pi*x), cos(2*pi*x),0]
// see 17/10/11

#include "math.h"
using namespace std;
using namespace MathematicalConstants;

namespace OneDMicromagSetup
{

  using namespace OneDMicromagSetup;

  // Time stepping parameters
  double t_max = 5;
  double dt = 0.01;

  // Number of elements
  double n_x_elements = 40;

  /// Compute all components of the exact solution at time t, position x
  void exact_solution(const double& t, const Vector<double>& x,
		      Vector<double>& solution)
  {
    // phi
    solution[0] = 2.0*cos(t)*sin(2*Pi*x[0]);

    // x,y,z components of M respectively
    solution[1] = cos(2*Pi*x[0])*cos(t);
    solution[1] = cos(2*Pi*x[0])*cos(t);
    solution[2] = 0.0;
  }

  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {
    // Source function to exactly cancel out contributions from not being a real exact solution
    source[0] = -cos(Pi*x[0]*2.0)*sin(t)-cos(Pi*x[0]*2.0)*cos(t)*((Pi*Pi)*pow(cos(Pi*x[0]*2.0),2.0)*pow(cos(t),2.0)*4.0-cos(Pi*x[0]*2.0)*cos(t)*(Pi*cos(Pi*x[0]*2.0)*cos(t)*4.0+(Pi*Pi)*cos(Pi*x[0]*2.0)*cos(t)*4.0));
    source[1] = -cos(Pi*x[0]*2.0)*sin(t)+cos(Pi*x[0]*2.0)*cos(t)*((Pi*Pi)*pow(cos(Pi*x[0]*2.0),2.0)*pow(cos(t),2.0)*4.0-cos(Pi*x[0]*2.0)*cos(t)*(Pi*cos(Pi*x[0]*2.0)*cos(t)*4.0+(Pi*Pi)*cos(Pi*x[0]*2.0)*cos(t)*4.0));
    source[2] = (Pi*Pi)*pow(cos(Pi*x[0]*2.0),2.0)*pow(cos(t),2.0)*-4.0+cos(Pi*x[0]*2.0)*cos(t)*(Pi*cos(Pi*x[0]*2.0)*cos(t)*4.0+(Pi*Pi)*cos(Pi*x[0]*2.0)*cos(t)*4.0);
  }

}; // End of namespace
