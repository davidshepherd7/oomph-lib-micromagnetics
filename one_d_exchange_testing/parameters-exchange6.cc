// parameters for the solution with M = t*t*[x, x, x]
//
// See 27/9/11 (2) for details

namespace OneDMicromagSetup
{
  using namespace MathematicalConstants;
  using namespace OneDMicromagSetup;

  // Time stepping parameters
  double t_max = 3;
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
    exact[0] = (4/3)*x[0]*x[0]*x[0]*t*t*Pi;

    // x, y, z components of M respectively
    exact[1] = x[0]*x[0]*t*t;
    exact[2] = x[0]*x[0]*t*t;
    exact[3] = x[0]*x[0]*t*t;

    // x,y,z components of H_exchange respectively
    exact[4] = exchange_coeff(t,x)*2*t*t;
    exact[5] = exchange_coeff(t,x)*2*t*t;
    exact[6] = exchange_coeff(t,x)*2*t*t;
  }

  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {

  source[0] = t*(x[0]*x[0])*2.0-llg_damping_coeff(t,x)*(t*t)*(x[0]*x[0])*((t*t)*(x[0]*x[0])*(exchange_coeff(t,x)*(t*t)*2.0-Pi*(t*t)*(x[0]*x[0])*4.0)-exchange_coeff(t,x)*(t*t*t*t)*(x[0]*x[0])*2.0)*2.0;
  source[1] = t*(x[0]*x[0])*2.0+llg_precession_coeff(t,x)*((t*t)*(x[0]*x[0])*(exchange_coeff(t,x)*(t*t)*2.0-Pi*(t*t)*(x[0]*x[0])*4.0)-exchange_coeff(t,x)*(t*t*t*t)*(x[0]*x[0])*2.0)+llg_damping_coeff(t,x)*(t*t)*(x[0]*x[0])*((t*t)*(x[0]*x[0])*(exchange_coeff(t,x)*(t*t)*2.0-Pi*(t*t)*(x[0]*x[0])*4.0)-exchange_coeff(t,x)*(t*t*t*t)*(x[0]*x[0])*2.0);
  source[2] = t*(x[0]*x[0])*2.0-llg_precession_coeff(t,x)*((t*t)*(x[0]*x[0])*(exchange_coeff(t,x)*(t*t)*2.0-Pi*(t*t)*(x[0]*x[0])*4.0)-exchange_coeff(t,x)*(t*t*t*t)*(x[0]*x[0])*2.0)+llg_damping_coeff(t,x)*(t*t)*(x[0]*x[0])*((t*t)*(x[0]*x[0])*(exchange_coeff(t,x)*(t*t)*2.0-Pi*(t*t)*(x[0]*x[0])*4.0)-exchange_coeff(t,x)*(t*t*t*t)*(x[0]*x[0])*2.0);

  }
}; // End of namespace
