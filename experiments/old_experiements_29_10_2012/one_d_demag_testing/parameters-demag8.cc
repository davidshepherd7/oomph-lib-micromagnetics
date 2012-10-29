// parameters for the solution M = cos(2*pi*x)* [cos(t), sin(t), 0]

// see ??ds for details

# include "math.h"
using namespace std;

namespace OneDMicromagSetup
{
  using namespace MathematicalConstants;
  using namespace OneDMicromagSetup;

  // Time stepping parameters
  double t_max = 2;
  double dt = 0.0001;

  // Number of elements
  unsigned n_x_elements = 5;

  // Prototypes for coeff functions
  double llg_damping_coeff(const double& t, const Vector<double>& x);
  double llg_precession_coeff(const double& t, const Vector<double>& x);
  double exchange_coeff(const double& t, const Vector<double>& x);

  void exact_solution(const double& t, const Vector<double>& x, Vector<double>& solution)
  {
    // Phi exact solution
    //??ds for piecewise constant llg_precession_coeff
    solution[0] = llg_precession_coeff(x,t)*2*sin(2*Pi*x[0])*cos(t);

    // x, y, z components of M respectively
    solution[1] = llg_precession_coeff(x,t)*cos(2*Pi*x[0])*cos(t);
    solution[2] = llg_precession_coeff(x,t)*cos(2*Pi*x[0])*sin(t);
    solution[3] = 0.0;
  }

  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {
    // Source function to exactly cancel out contributions from solution not being a real solution

source[0] = -llg_precession_coeff(t,x)*cos(Pi*x[0]*2.0)*sin(t)+llg_damping_coeff(t,x)*llg_precession_coeff(t,x)*cos(Pi*x[0]*2.0)*sin(t)*(llg_precession_coeff(t,x)*cos(Pi*x[0]*2.0)*sin(t)*(Pi*llg_precession_coeff(t,x)*cos(Pi*x[0]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*llg_precession_coeff(t,x)*cos(Pi*x[0]*2.0)*cos(t)*4.0)-(Pi*Pi)*exchange_coeff(t,x)*(llg_precession_coeff(t,x)*llg_precession_coeff(t,x))*pow(cos(Pi*x[0]*2.0),2.0)*cos(t)*sin(t)*4.0);
  source[1] = llg_precession_coeff(t,x)*cos(Pi*x[0]*2.0)*cos(t)-llg_damping_coeff(t,x)*llg_precession_coeff(t,x)*cos(Pi*x[0]*2.0)*cos(t)*(llg_precession_coeff(t,x)*cos(Pi*x[0]*2.0)*sin(t)*(Pi*llg_precession_coeff(t,x)*cos(Pi*x[0]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*llg_precession_coeff(t,x)*cos(Pi*x[0]*2.0)*cos(t)*4.0)-(Pi*Pi)*exchange_coeff(t,x)*(llg_precession_coeff(t,x)*llg_precession_coeff(t,x))*pow(cos(Pi*x[0]*2.0),2.0)*cos(t)*sin(t)*4.0);
  source[2] = llg_precession_coeff(t,x)*(llg_precession_coeff(t,x)*cos(Pi*x[0]*2.0)*sin(t)*(Pi*llg_precession_coeff(t,x)*cos(Pi*x[0]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*llg_precession_coeff(t,x)*cos(Pi*x[0]*2.0)*cos(t)*4.0)-(Pi*Pi)*exchange_coeff(t,x)*(llg_precession_coeff(t,x)*llg_precession_coeff(t,x))*pow(cos(Pi*x[0]*2.0),2.0)*cos(t)*sin(t)*4.0);

  }

}; // end of namespace
