// parameters for M = cos(2*pi*x)* [cos(t), sin(t), 0]

#include "math.h"
using namespace std;
using namespace MathematicalConstants;

namespace OneDMicromagSetup
{

  using namespace OneDMicromagSetup;

  // Time stepping parameters
  double t_max = 5;
  double dt = 0.05;

  // Number of elements
  unsigned n_x_elements = 40;


  // Prototypes for coeff functions
  double llg_damping_coeff(const double& t, const Vector<double>& x);
  double llg_precession_coeff(const double& t, const Vector<double>& x);
  double exchange_coeff(const double& t, const Vector<double>& x);

  /// Compute all components of the exact solution at time t, position x
  void exact_solution(const double& t, const Vector<double>& x,
		      Vector<double>& solution)
  {
    // phi
    solution[0] = 2*sin(2*Pi*x[0])*cos(t);

    // x,y,z components of M respectively
    solution[1] = cos(2*Pi*x[0])*cos(t);
    solution[2] = cos(2*Pi*x[0])*sin(t);
    solution[3] = 0.0;

    // x,y,z components of H_exchange respectively
    solution[4] = (Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[0]*2.0)*cos(t)*-4.0;
    solution[5] = (Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[0]*2.0)*sin(t)*-4.0;
    solution[6] = 0.0;
  }

  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {

    source[0] = -cos(Pi*x[0]*2.0)*sin(t)+llg_damping_coeff(t,x)*cos(Pi*x[0]*2.0)*sin(t)*(cos(Pi*x[0]*2.0)*sin(t)*(Pi*cos(Pi*x[0]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[0]*2.0)*cos(t)*4.0)-(Pi*Pi)*exchange_coeff(t,x)*pow(cos(Pi*x[0]*2.0),2.0)*cos(t)*sin(t)*4.0);
    source[1] = cos(Pi*x[0]*2.0)*cos(t)-llg_damping_coeff(t,x)*cos(Pi*x[0]*2.0)*cos(t)*(cos(Pi*x[0]*2.0)*sin(t)*(Pi*cos(Pi*x[0]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[0]*2.0)*cos(t)*4.0)-(Pi*Pi)*exchange_coeff(t,x)*pow(cos(Pi*x[0]*2.0),2.0)*cos(t)*sin(t)*4.0);
    source[2] = llg_precession_coeff(t,x)*(cos(Pi*x[0]*2.0)*sin(t)*(Pi*cos(Pi*x[0]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[0]*2.0)*cos(t)*4.0)-(Pi*Pi)*exchange_coeff(t,x)*pow(cos(Pi*x[0]*2.0),2.0)*cos(t)*sin(t)*4.0);


  }

}; // End of namespace
