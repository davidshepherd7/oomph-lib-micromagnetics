// parameters for the solution with M constant in space: M = [sin(t), cos(t), 0]
//
// See 5/10/11 for details

# include "math.h"

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

  void exact_solution(const double& t, const Vector<double>& x, Vector<double>& exact)
  {
    // phi
    exact[0] = x[0];

    // x,y,z of M
    exact[1] = sin(t);
    exact[2] = cos(t);
    exact[3] = 0.0;
  }

  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {
    source[0] = cos(t)+llg_damping_coeff(t,x)*pow(cos(t),2.0);
    source[1] = -sin(t)-llg_damping_coeff(t,x)*cos(t)*sin(t);
    source[2] = llg_precession_coeff(t,x)*cos(t);
  }

}; // End of namespace
