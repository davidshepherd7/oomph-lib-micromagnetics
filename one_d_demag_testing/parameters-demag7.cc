// Parameters for M = (x^2 - x) * [cos(t), sin(t), 0] (M=0 at ends)
// See 7/10/11 (3) for details

using namespace std;
using namespace MathematicalConstants;

namespace OneDMicromagSetup
{
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

  void exact_solution(const double& t, const Vector<double>& x, Vector<double>& solution)
  {
    // Phi exact solution
    solution[0] =  4.0*Pi*(x[0]*x[0]*x[0]/3 - x[0]*x[0]/2)*cos(t);

    // x, y, z components of M respectively
    solution[1] = (x[0]*x[0] - x[0])*cos(t);
    solution[2] = (x[0]*x[0] - x[0])*sin(t);
    solution[3] = 0.0;
  }

  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {

    source[0] = sin(t)*(x[0]-x[0]*x[0])-llg_damping_coeff(t,x)*pow(sin(t),2.0)*pow(x[0]-x[0]*x[0],2.0)*(exchange_coeff(t,x)*cos(t)*2.0-Pi*cos(t)*pow(x[0]*2.0-1.0,2.0));
  source[1] = -cos(t)*(x[0]-x[0]*x[0])+llg_damping_coeff(t,x)*cos(t)*sin(t)*pow(x[0]-x[0]*x[0],2.0)*(exchange_coeff(t,x)*cos(t)*2.0-Pi*cos(t)*pow(x[0]*2.0-1.0,2.0));
  source[2] = llg_precession_coeff(t,x)*sin(t)*(x[0]-x[0]*x[0])*(exchange_coeff(t,x)*cos(t)*2.0-Pi*cos(t)*pow(x[0]*2.0-1.0,2.0));


  }

}; // end of namespace
