// parameters for the solution with M = [0, x*cos(t), 0 ]
//
// See 27/11/11 (4) for details

// note that we do not use the simplest solution (phi==0), instead we choose the bc s.t. phi=x.

using namespace std;

namespace OneDMicromagSetup
{
  using namespace MathematicalConstants;
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

  void exact_solution(const double& t, const Vector<double>& x, Vector<double>& solution)
  {
    //phi
    solution[0] = x[0];

    // x, y, z components of M
    solution[1] = 0.0;
    solution[2] = x[0]*cos(t);
    solution[3] = 0.0;
  }

  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {
    // Source function to exactly cancel out contributions from not being a real exact solution

  source[0] = llg_damping_coeff(t,x)*(x[0]*x[0])*pow(cos(t),2.0);
  source[1] = -x[0]*sin(t);
  source[2] = llg_precession_coeff(t,x)*x[0]*cos(t);

  }

}; // end of namespace
