// parameters for the solution with M constant in space: M = [ sin(omega*t), cos(omega*t), 0]
// 
// See 5/10/11 for details

namespace OneDMicromagSetup
{
  using namespace OneDMicromagSetup;

  // Time stepping parameters
  double t_max = 5;
  double dt = 0.05;

  // Number of elements
  double n_x_elements = 40;
  
  // Time parameter in solutions
  double omega = 1;

  void exact_M_solution(const double& t, const Vector<double>& x, Vector<double>& M)
  {
    M[0] = sin(omega*t);
    M[1] = cos(omega*t);
    M[2] = 0.0;
  }

  double exact_phi_solution(const double& t, const Vector<double>& x)
  {
    return x[0];
  }
  
  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {
    source[0] = omega*cos(omega*t) + cos(omega*t)*cos(omega*t);
    source[1] = -omega*sin(omega*t) - sin(omega*t)*cos(omega*t);
    source[2] = cos(omega*t);
  }

}; // End of namespace

#include "baseparameters.cc"
