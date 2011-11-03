// parameters for the solution with M = t*t*(x - x*x)[1,1,0]
//
// See 27/9/11 (2) for details

namespace OneDMicromagSetup
{
  using namespace MathematicalConstants;
  using namespace OneDMicromagSetup;

  // Time stepping parameters
  double t_max = 6;
  double dt = 0.05;

  // Number of elements
  unsigned n_x_elements = 10;

  // Prototypes for coeff functions
  double llg_damping_coeff(const double& t, const Vector<double>& x);
  double llg_precession_coeff(const double& t, const Vector<double>& x);
  double exchange_coeff(const double& t, const Vector<double>& x);

  void exact_solution(const double& t, const Vector<double>& x,
		      Vector<double>& exact)
  {
    // phi
    exact[0] = (-2)*Pi*t*t*(x[0]*x[0] - (2*x[0]*x[0]*x[0])/3);

    // x, y, z components of M respectivelx[1]
    exact[1] = -(t*t)*(x[0]-x[0]*x[0]);
    exact[2] = -(t*t)*(x[0]-x[0]*x[0]);
    exact[3] = 0.0;

    // x,y,z components of H_exchange respectively
    exact[4] = exchange_coeff(t,x)*2*t*t;
    exact[5] = exchange_coeff(t,x)*2*t*t;
    exact[6] = 0.0;
  }

  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {

  source[0] = t*x[0]*(x[0]-1.0)*(Pi*llg_damping_coeff(t,x)*(t*t*t*t*t)*(x[0]*x[0])*2.0-Pi*llg_damping_coeff(t,x)*(t*t*t*t*t)*(x[0]*x[0]*x[0])*4.0+Pi*llg_damping_coeff(t,x)*(t*t*t*t*t)*(x[0]*x[0]*x[0]*x[0])*2.0+1.0)*2.0;
  source[1] = t*x[0]*(x[0]-1.0)*(Pi*llg_damping_coeff(t,x)*(t*t*t*t*t)*(x[0]*x[0])*2.0-Pi*llg_damping_coeff(t,x)*(t*t*t*t*t)*(x[0]*x[0]*x[0])*4.0+Pi*llg_damping_coeff(t,x)*(t*t*t*t*t)*(x[0]*x[0]*x[0]*x[0])*2.0-1.0)*-2.0;
  source[2] = Pi*llg_precession_coeff(t,x)*(t*t*t*t)*(x[0]*x[0])*pow(x[0]-1.0,2.0)*4.0;


  }
}; // End of namespace
