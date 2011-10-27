// parameters for M = cos(t) *[cos(2*pi*x)*sin(2*pi*y)*sin(2*pi*z),
//                             sin(2*pi*x)*cos(2*pi*y)*sin(2*pi*z),
//                             sin(2*pi*x)*sin(2*pi*y)*cos(2*pi*z)]

// see 27/10/11

#include "math.h"
using namespace std;
using namespace MathematicalConstants;

namespace ThreeDMicromagSetup
{

  using namespace ThreeDMicromagSetup;

  // Time stepping parameters
  double t_max = 5;
  double dt = 0.1;

  // Number of elements
  unsigned n_x = 6;
  unsigned n_y = 6;
  unsigned n_z = 6;

  // Prototypes for coeff functions
  double llg_damping_coeff(const double& t, const Vector<double>& x);
  double llg_precession_coeff(const double& t, const Vector<double>& x);
  double exchange_coeff(const double& t, const Vector<double>& x);

  /// Compute all components of the exact solution at time t, position x
  void exact_solution(const double& t, const Vector<double>& x,
		      Vector<double>& solution)
  {
    // phi
    solution[0] = 2*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t);

    // M
    solution[1] = cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t);
    solution[2] = cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t);
    solution[3] = cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t);
  }

  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {
source[0] = -llg_precession_coeff(t,x)*(cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*(Pi*cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*1.2E1)-cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*(Pi*cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*1.2E1))-llg_damping_coeff(t,x)*(cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*(cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*(Pi*cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*1.2E1)-cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*(Pi*cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*1.2E1))+cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*(cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*(Pi*cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*1.2E1)-cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*(Pi*cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*1.2E1)))-cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*sin(t);
  source[1] = llg_precession_coeff(t,x)*(cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*(Pi*cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*1.2E1)-cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*(Pi*cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*1.2E1))+llg_damping_coeff(t,x)*(cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*(cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*(Pi*cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*1.2E1)-cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*(Pi*cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*1.2E1))-cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*(cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*(Pi*cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*1.2E1)-cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*(Pi*cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*1.2E1)))-cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*sin(t);
  source[2] = -llg_precession_coeff(t,x)*(cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*(Pi*cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*1.2E1)-cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*(Pi*cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*1.2E1))+llg_damping_coeff(t,x)*(cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*(cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*(Pi*cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*1.2E1)-cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*(Pi*cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*1.2E1))+cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*(cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*(Pi*cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*1.2E1)-cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*cos(t)*(Pi*cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[1]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[2]*2.0)*cos(t)*1.2E1)))-cos(Pi*x[2]*2.0)*sin(Pi*x[0]*2.0)*sin(Pi*x[1]*2.0)*sin(t);

  }

}; // End of namespace
