// parameters for M = t^2*[1,1,1]

#include "math.h"
using namespace std;
using namespace MathematicalConstants;

namespace OneDMicromagSetup
{

  using namespace OneDMicromagSetup;

  // Time stepping parameters
  double t_max = 3;
  double dt = 0.05;

  // Number of elements
  unsigned n_x_elements = 3;


  // Prototypes for coeff functions
  double llg_damping_coeff(const double& t, const Vector<double>& x);
  double llg_precession_coeff(const double& t, const Vector<double>& x);
  double exchange_coeff(const double& t, const Vector<double>& x);

  /// Compute all components of the exact solution at time t, position x
  void exact_solution(const double& t, const Vector<double>& x,
		      Vector<double>& solution)
  {
    // phi
    solution[0] = -x[0];

    // x,y,z components of M respectively
    solution[1] = 0.0;
    solution[2] = t*t;
    solution[3] = t*t;
  }

  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {

    source[0] = t*2.0-llg_damping_coeff(t,x)*(t*t*t*t)*2.0;
    source[1] = t*2.0+llg_damping_coeff(t,x)*(t*t*t*t)+llg_precession_coeff(t,x)*(t*t);
    source[2] = t*2.0+llg_damping_coeff(t,x)*(t*t*t*t)-llg_precession_coeff(t,x)*(t*t);

  }

}; // End of namespace
