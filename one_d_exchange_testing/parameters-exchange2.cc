// parameters for M = [1,0,0]

// see ??ds for more detail

#include "math.h"
using namespace std;
using namespace MathematicalConstants;

namespace OneDMicromagSetup
{

  using namespace OneDMicromagSetup;

  // Time stepping parameters
  double t_max = 5;
  double dt = 0.05;

  // Number of elements
  unsigned n_x_elements = 1;

  /// Compute all components of the exact solution at time t, position x
  void exact_solution(const double& t, const Vector<double>& x,
			       Vector<double>& solution)
  {
    // phi
    solution[0] = 0.0;

    // x,y,z components of M respectively
    solution[1] = 1.0;
    solution[2] = 0.0;
    solution[3] = 0.0;
  }

  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {
    // Source function to exactly cancel out contributions from not being a real exact solution
    source[0] = 0.0;
    source[1] = 0.0;
    source[2] = 0.0;
  }
}; // End of namespace
