// parameters for the solution with [M = 0, x cos(omega*t), 0 ]
// 
// See 27/11/11 (4) for details

using namespace std;
using namespace MathematicalConstants;

namespace OneDMicromagSetup
{

  using namespace OneDMicromagSetup;

  // Time stepping parameters
  double t_max = 5;
  double dt = 0.05;

  // Number of elements
  double n_x_elements = 40;

  // double alpha = 0.5;   // Gibert damping constant
  // double gamma = 0.221;   // Electromagnetic ratio
  double omega = 1; // Time parameter in solution
  double damping_const = 1;
  double precession_const = 0;
 
  void exact_M_solution(const double& t, const Vector<double>& x, Vector<double>& M)
  {
    M[0] = 0.0;
    M[1] = x[0]*cos(omega*t);
    M[2] = 0.0;
  }

  double exact_phi_solution(const double& t, const Vector<double>& x)
  {
    return x[0];
  }
  
  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {
    // Source function to exactly cancel out contributions from not being a real exact solution
    source[0] = damping_const* -cos(omega*t)*x[0]*x[0];
    source[1] = -omega*sin(omega*t)*x[0];
    source[2] = precession_const * -cos(omega*t)*x[0];
  }

  /////////////////////////////////////////////////////////////////////////////////////////

  // Get the saturisation magnetisation at position x
  double sat_mag(const double& t, const Vector<double>& x)
  {
    Vector<double> M(3,0.0);
    exact_M_solution(t,x,M);
    return M[0]*M[0] + M[1]*M[1] + M[2]*M[2];
  }

  void cryst_anis_field(const double& t, const Vector<double>& x, const Vector<double>& m, Vector<double>& H_cryst_anis)
  {
    H_cryst_anis[0] = 0.0;
    H_cryst_anis[1] = 0.0;
    H_cryst_anis[2] = 0.0;
  }

  void initial_M(const Vector<double>& x, Vector<double>& M)
  {
    exact_M_solution(0.0,x,M);
  }

  void applied_field(const double& t, const Vector<double>& x, Vector<double>& H_applied)
  {
    H_applied[0] = 0.0;
    H_applied[1] = 0.0;
    H_applied[2] = 0.0;
  }

  double boundary_phi(const double& t, const Vector<double>& x)
  {
    return exact_phi_solution(t,x);
  }

  void source_function(const double& t, const Vector<double>& x, double& source)
  {
    // Just set to zero unless needed for testing - poisson source
    source = 0.0;
  }

  // The coefficient of the damping term (Mx(MxH)) of the Landau-Lifschitz-Gilbert equation
  double llg_damping_coeff(const Vector<double>& x)
  {   
    //??ds ifdef error checking then check ms is nonzero
    //    return (1/(1+alpha*alpha)) * (gamma/sat_mag(t,x));
    // return (1/(1+alpha*alpha)) * (gamma);
    return damping_const;
  }

  // The coefficient of the precession term (MxH) of the Landau-Lifschitz-Gilbert equation
  double llg_precession_coeff(const Vector<double>& x)
  {
    return precession_const; ///(1+alpha*alpha); 
  }

}; // End of namespace
