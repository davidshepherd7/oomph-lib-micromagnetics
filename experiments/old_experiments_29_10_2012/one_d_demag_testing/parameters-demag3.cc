// parameters for the solution with M = [x*cos(t), (1-x)*cos(t), 0]
//
// See 27/9/11 (2) for details

namespace OneDMicromagSetup
{
  using namespace MathematicalConstants;
  using namespace OneDMicromagSetup;

  // Time stepping parameters
  double t_max = 10;
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
    exact[0] = 2*Pi*cos(t)*x[0]*x[0];

    // x, y, z components of M respectively
    exact[1] = x[0]*cos(t);
    exact[2] = (1-x[0])*cos(t);
    exact[3] = 0.0;
  }

  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {
    source[0] = -x[0]*sin(t)+Pi*llg_damping_coeff(t,x)*x[0]*pow(cos(t),3.0)*pow(x[0]-1.0,2.0)*4.0;
    source[1] = sin(t)*(x[0]-1.0)+Pi*llg_damping_coeff(t,x)*(x[0]*x[0])*pow(cos(t),3.0)*(x[0]-1.0)*4.0;
    source[2] = Pi*llg_precession_coeff(t,x)*x[0]*pow(cos(t),2.0)*(x[0]-1.0)*-4.0;


  }
}; // End of namespace
