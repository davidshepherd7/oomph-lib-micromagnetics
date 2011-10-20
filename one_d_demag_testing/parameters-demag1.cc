// parameters for the solution M = cos(2*pi*x)* [cos(t), sin(t), 0]

// see ??ds for details

# include "math.h"
using namespace std;

namespace OneDMicromagSetup
{
  using namespace MathematicalConstants;
  using namespace OneDMicromagSetup;

  // Time stepping parameters
  double t_max = 4;
  double dt = 0.05;

  // Number of elements
  double n_x_elements = 40;

  // Prototypes for coeff functions
  double llg_damping_coeff(const double& t, const Vector<double>& x);
  double llg_precession_coeff(const double& t, const Vector<double>& x);
  double exchange_coeff(const double& t, const Vector<double>& x);

  void exact_solution(const double& t, const Vector<double>& x, Vector<double>& solution)
  {
    // Phi exact solution
    solution[0] = 2*sin(2*Pi*x[0])*cos(t);

    // x, y, z components of M respectively
    solution[1] = cos(2*Pi*x[0])*cos(t);
    solution[2] = cos(2*Pi*x[0])*sin(t);
    solution[3] = 0.0;
  }

  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {
    // Source function to exactly cancel out contributions from solution not being a real solution

    source[0] = -cos(Pi*x[0]*2.0)*sin(t)+llg_damping_coeff(t,x)*pow(cos(Pi*x[0]*2.0),2.0)*pow(sin(t),2.0)*(Pi*cos(Pi*x[0]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[0]*2.0)*cos(t)*4.0);
    source[1] = cos(Pi*x[0]*2.0)*cos(t)-llg_damping_coeff(t,x)*pow(cos(Pi*x[0]*2.0),2.0)*cos(t)*sin(t)*(Pi*cos(Pi*x[0]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[0]*2.0)*cos(t)*4.0);
    source[2] = llg_precession_coeff(t,x)*cos(Pi*x[0]*2.0)*sin(t)*(Pi*cos(Pi*x[0]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*cos(Pi*x[0]*2.0)*cos(t)*4.0);

  }

}; // end of namespace
