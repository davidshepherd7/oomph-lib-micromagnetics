// parameters for the solution M_x = cos(2*pi*x)* [cos(t), sin(t), 0]

// see ??ds for details

# include "math.h"
using namespace std;

namespace OneDMicromagSetup
{
  using namespace MathematicalConstants;
  using namespace OneDMicromagSetup;

  // Time stepping parameters
  double t_max = 4;
  double dt = 0.05;

  // Number of elements
  double n_x_elements = 40;

  void exact_solution(const double& t, const Vector<double>& x, Vector<double>& solution)
  {
    // Phi exact solution
    solution[0] = 2*sin(2*Pi*x[0])*cos(t);

    // x, y, z components of M respectively
    solution[1] = cos(2*Pi*x[0])*cos(t);
    solution[2] = cos(2*Pi*x[0])*sin(t);
    solution[3] = 0.0;
  }

  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {
    // Source function to exactly cancel out contributions from solution not being a real solution

    source[0] = 4*Pi*cos(2*Pi*x[0])*cos(2*Pi*x[0])*cos(2*Pi*x[0])*cos(t)*sin(t)*sin(t) 
      - cos(2*Pi*x[0])*sin(t);
    source[1] = cos(2*Pi*x[0])*cos(t) 
      - 4*Pi*cos(2*Pi*x[0])*cos(2*Pi*x[0])*cos(2*Pi*x[0])*cos(t)*cos(t)*sin(t);
    source[2] = 4*Pi*cos(2*Pi*x[0])*cos(2*Pi*x[0])*cos(t)*sin(t);
  }

}; // end of namespace
