// Base parameters for all one_d_demag_testing parameters files

using namespace std;
using namespace MathematicalConstants;

namespace OneDMicromagSetup
{
  using namespace OneDMicromagSetup;

  // Derived quantities:
  //==================================================

  double alpha = 0.7;   // Gibert damping constant
  double gamma = 0.221E-8;   // Electromagnetic ratio

  // The coefficient of the damping term of the Landau-Lifschitz-Gilbert equation
  //??ds ifdef error checking then check ms is nonzero?
  double llg_damping_coeff(const double& t, const Vector<double>& x)    
  {return (1/(1+alpha*alpha))*(alpha*gamma/1.0);}

  // The coefficient of the precession term of the Landau-Lifschitz-Gilbert equation
  double llg_precession_coeff(const Vector<double>& x)
  {return (1+alpha*alpha);}

  //Get the saturisation magnetisation at position x
  double sat_mag(const double& t, const Vector<double>& x)
  {
    Vector<double> M(3,0.0);
    exact_M_solution(t,x,M);
    return sqrt(M[0]*M[0] + M[1]*M[1] + M[2]*M[2]);
  }

  void cryst_anis_field(const double& t, const Vector<double>& x, const Vector<double>& m, Vector<double>& H_cryst_anis)
  {
    H_cryst_anis[0] = 0.0;    
    H_cryst_anis[1] = 0.0;
    H_cryst_anis[2] = 0.0;
  }

  void initial_M(const double& t, const Vector<double>& x, Vector<double>& M)
  {
    exact_M_solution(t,x,M);
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
    // Just set to zero unless needed for testing
    source = 0.0;
  }

}; // End of namespace
