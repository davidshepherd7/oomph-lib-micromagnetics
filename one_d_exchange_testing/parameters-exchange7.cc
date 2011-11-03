// parameters for the solution with M = t*(x*x*x/3 - x*x/2)*[1,1,0]
//
// See 3/11/11 (3) for details

// note that dMx/dx = 0 at boundaries

namespace OneDMicromagSetup
{
  using namespace MathematicalConstants;
  using namespace OneDMicromagSetup;

  // Time stepping parameters
  double t_max = 6;
  double dt = 0.05;

  // Number of elements
  unsigned n_x_elements = 40;

  // Prototypes for coeff functions
  double llg_damping_coeff(const double& t, const Vector<double>& x);
  double llg_precession_coeff(const double& t, const Vector<double>& x);
  double exchange_coeff(const double& t, const Vector<double>& x);

  void exact_solution(const double& t, const Vector<double>& x,
		      Vector<double>& exact)
  {
    // phi
    exact[0] = (1.0/3.0)*Pi*t* x[0]*x[0]*x[0]* (x[0] - 2);

    // x[0], x[1], x[2] components of M respectivelx[1]
    exact[1] = -t*( (x[0]*x[0])/2.0 - (x[0]*x[0]*x[0])/3.0 );
    exact[2] = -t*( (x[0]*x[0])/2.0 - (x[0]*x[0]*x[0])/3.0 );
    exact[3] = 0.0;

    // x,y,z components of H_exchange respectively
    exact[4] = exchange_coeff(t,x)*t*(2*x[0] - 1);
    exact[5] = exchange_coeff(t,x)*t*(2*x[0] - 1);
    exact[6] = 0.0;
  }

  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {

  source[0] = (x[0]*x[0])*(x[0]*2.0-3.0)*(Pi*llg_damping_coeff(t,x)*(t*t*t)*(x[0]*x[0]*x[0]*x[0])*9.0-Pi*llg_damping_coeff(t,x)*(t*t*t)*(x[0]*x[0]*x[0]*x[0]*x[0])*1.2E1+Pi*llg_damping_coeff(t,x)*(t*t*t)*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0])*4.0+9.0)*(1.0/5.4E1);
  source[1] = (x[0]*x[0])*(x[0]*2.0-3.0)*(Pi*llg_damping_coeff(t,x)*(t*t*t)*(x[0]*x[0]*x[0]*x[0])*9.0-Pi*llg_damping_coeff(t,x)*(t*t*t)*(x[0]*x[0]*x[0]*x[0]*x[0])*1.2E1+Pi*llg_damping_coeff(t,x)*(t*t*t)*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0])*4.0-9.0)*(-1.0/5.4E1);
  source[2] = Pi*llg_precession_coeff(t,x)*(t*t)*(x[0]*x[0]*x[0]*x[0])*pow(x[0]*2.0-3.0,2.0)*(1.0/9.0);

  }
}; // End of namespace
