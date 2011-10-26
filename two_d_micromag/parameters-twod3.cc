// parameters for M = cos(t) *[cos(2*pi*x)*sin(2*pi*y), cos(2*pi*y)*sin(2*pi*x),0]
// see 26/10/11

#include "math.h"
using namespace std;
using namespace MathematicalConstants;

namespace TwoDMicromagSetup
{

  using namespace TwoDMicromagSetup;

  // Time stepping parameters
  double t_max = 5;
  double dt = 0.05;

  // Number of elements
  unsigned n_x = 10;
  unsigned n_y = 10;



  // Prototypes for coeff functions
  double llg_damping_coeff(const double& t, const Vector<double>& x);
  double llg_precession_coeff(const double& t, const Vector<double>& x);
  double exchange_coeff(const double& t, const Vector<double>& x);

  /// Compute all components of the exact solution at time t, position x
  void exact_solution(const double& t, const Vector<double>& x,
		      Vector<double>& solution)
  {
    // phi
    solution[0] =

    // x,y,z components of M respectively
    solution[1] = cos(2*Pi*x[0])*sin(2*Pi*x[1])*cos(t);
    solution[2] = cos(2*Pi*x[1])*sin(2*Pi*x[0])*cos(t);
    solution[3] = 0.0;
  }

  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {
source[0] = -cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(t)-llg_damping_coeff(t,x)*cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*cos(t)*(cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*(Pi*cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*cos(t)*8.0)-cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*cos(t)*(Pi*cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*8.0));
  source[1] = -cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(t)+llg_damping_coeff(t,x)*cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*(cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*(Pi*cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*cos(t)*8.0)-cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*cos(t)*(Pi*cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*8.0));
  source[2] = -llg_precession_coeff(t,x)*(cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*(Pi*cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*cos(t)*8.0)-cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*cos(t)*(Pi*cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*8.0));

  }

}; // End of namespace
