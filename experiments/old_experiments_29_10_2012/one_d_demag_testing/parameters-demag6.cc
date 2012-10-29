// Parameters for M = sin(2*pi*x) * [cos(t), sin(t), 0] (M=0 at ends)

// See 7/10/11 (2) for details

using namespace std;
using namespace MathematicalConstants;

namespace OneDMicromagSetup
{
  using namespace OneDMicromagSetup;

  // Time stepping parameters
  double t_max = 3;
  double dt = 0.05;

  // Number of elements
  unsigned n_x_elements = 41;

  // Prototypes for coeff functions
  double llg_damping_coeff(const double& t, const Vector<double>& x);
  double llg_precession_coeff(const double& t, const Vector<double>& x);
  double exchange_coeff(const double& t, const Vector<double>& x);

  void exact_solution(const double& t, const Vector<double>& x, Vector<double>& exact)
  {
    // phi
    exact[0] = -2.0*cos(2*Pi*x[0])*cos(t);

    // x,y,z components of M
    exact[1] = sin(2*Pi*x[0])*cos(t);
    exact[2] = sin(2*Pi*x[0])*sin(t);
    exact[3] = 0.0;
  }

  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {
    // Source function to exactly cancel out contributions from solution not being a real solution

source[0] = -sin(Pi*x[0]*2.0)*sin(t)+llg_damping_coeff(t,x)*pow(sin(Pi*x[0]*2.0),2.0)*pow(sin(t),2.0)*(Pi*sin(Pi*x[0]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*sin(Pi*x[0]*2.0)*cos(t)*4.0);
  source[1] = sin(Pi*x[0]*2.0)*cos(t)-llg_damping_coeff(t,x)*pow(sin(Pi*x[0]*2.0),2.0)*cos(t)*sin(t)*(Pi*sin(Pi*x[0]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*sin(Pi*x[0]*2.0)*cos(t)*4.0);
  source[2] = llg_precession_coeff(t,x)*sin(Pi*x[0]*2.0)*sin(t)*(Pi*sin(Pi*x[0]*2.0)*cos(t)*4.0+(Pi*Pi)*exchange_coeff(t,x)*sin(Pi*x[0]*2.0)*cos(t)*4.0);

  }

}; // end of parameters
