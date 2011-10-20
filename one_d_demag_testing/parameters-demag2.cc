// parameters for the solution with M = [x*cos(t), 0, 0]
//
// See 27/11/11 (1) for details

# include "math.h"
using namespace std;

namespace OneDMicromagSetup
{
  using namespace MathematicalConstants;
  using namespace OneDMicromagSetup;

  // Time stepping parameters
  double t_max = 5;
  double dt = 0.05;

  // Number of elements
  double n_x_elements = 40;

  void exact_solution(const double& t, const Vector<double>& x,
		      Vector<double>& solution)
  {
    // Phi exact solution
    solution[0] = 2*Pi*cos(t)*x[0]*x[0];

    // x, y, z components of M respectively
    solution[1] = x[0]*cos(t);
    solution[2] = 0.0;
    solution[3] = 0.0;
  }

  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {
    // Source function to exactly cancel out contributions from not being a real exact solution
    source[0] = -x[0]*sin(t);
    source[1] = 0.0;
    source[2] = 0.0;
  }

}; // end of namespace
