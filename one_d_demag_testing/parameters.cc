
using namespace std;
using namespace MathematicalConstants;

//==start_of_namespace================================================
/// Namespace for various functions
//====================================================================
namespace OneDMicromagSetup
{
  using namespace OneDMicromagSetup;

  // Time parameter
  double omega = 1;

  // Gibert damping constant
  double alpha = 0.5;

  // Electromagnetic ratio
  double gamma = 0.221;

  
  // Get the saturisation magnetisation at position x
  double sat_mag(const double& t, const Vector<double>& x)
  {
    return cos(2*Pi*x[0]);
  }

  void cryst_anis_field(const double& t, const Vector<double>& x, const Vector<double>& m, Vector<double>& H_cryst_anis)
  {
    H_cryst_anis[0] = 0;
    
    H_cryst_anis[1] = 0;

    H_cryst_anis[2] = 0;
  }

  void initial_M(const Vector<double>& x, Vector<double>& M)
  {
    //??ds assumes initial t=0
    M[0] = cos(2*Pi*x[0]);

    // Initialise y and z components of M to zero
    M[1] = 0.0;
    M[2] = 0.0;
  }

  void applied_field(const double& t, const Vector<double>& x, Vector<double>& H_applied)
  {
    double a= 4*Pi*cos(2*Pi*x[0]);
 
    H_applied[0] = a*cos(omega*t);
    H_applied[1] = 0.0;
    H_applied[2] = 0.0;

  }

  double boundary_phi(const double& t, const Vector<double>& x)
  {
    // Set all boundaries to zero for now
    //return 2*sin(2*Pi*x[0])*cos(omega*t);
    return 0.0;
  }

  
  void source_function(const double& t, const Vector<double>& x, double& source)
  {
    // Just set to zero unless needed for testing
    source = 0.0;
  }

  // The coefficient of the damping term of the Landau-Lifschitz-Gilbert equation
  double llg_damping_coeff(const double& t, const Vector<double>& x)
  {   
    //??ds ifdef error checking then check ms is nonzero
    return (1/(1+alpha*alpha)) * (gamma/sat_mag(t,x));
  }

  // The coefficient of the precession term of the Landau-Lifschitz-Gilbert equation
  double llg_precession_coeff(const Vector<double>& x)
  {
    return 1/(1+alpha*alpha); 
  }

  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {
    // Source function to exactly cancel out contributions 
    double p = 4*Pi*cos(2*Pi*x[0])*sin(omega*t)*cos(omega*t);

    source[0] = p*llg_damping_coeff(t,x)*cos(2*Pi*x[0])*sin(omega*t);
    source[1] = p*llg_damping_coeff(t,x)*cos(2*Pi*x[0])*cos(omega*t);
    source[2] = p*llg_precession_coeff(x);
  }

  void exact_M_solution(const double& t, const Vector<double>& x, Vector<double>& M)
  {
    double p = cos(2*Pi*x[0]);
    
    M[0] = p*cos(omega*t);
    M[1] = p*sin(omega*t);
    M[2] = 0.0;
  }

  double exact_phi_solution(const double& t, const Vector<double>& x)
  {
    return 2*sin(2*Pi*x[0])*cos(omega*t);
  }

}; // End of namespace
