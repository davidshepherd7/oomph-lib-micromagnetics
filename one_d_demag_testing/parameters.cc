// parameters for the solution M_x = 1

using namespace std;
using namespace MathematicalConstants;

namespace OneDMicromagSetup
{
  using namespace OneDMicromagSetup;

  // Time stepping parameters
  double t_max = 3;
  double dt = 0.05;

  // Number of elements
  double n_x_elements = 41;

  double alpha = 0.7;   // Gibert damping constant
  double gamma = 0.221E-8;   // Electromagnetic ratio
  double omega = 1; // Time parameter in solution

  // The coefficient of the damping term of the Landau-Lifschitz-Gilbert equation
  //??ds ifdef error checking then check ms is nonzero?
  //??ds just set sat_mag = 1 for now since coeffs shouldn't matter at this point
  double llg_damping_coeff(const double& t, const Vector<double>& x)    
  {return 1.0; //(1/(1+alpha*alpha))*(alpha*gamma/1.0);
  }

  // The coefficient of the precession term of the Landau-Lifschitz-Gilbert equation
  double llg_precession_coeff(const Vector<double>& x)
  {return 1.0; //(1+alpha*alpha);
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
  
  void llg_source_function(const double& t, const Vector<double>& x, Vector<double>& source)
  {
    // Source function to exactly cancel out contributions from solution not being a real solution

    source[0] = 4*Pi*cos(2*Pi*x[0])*cos(2*Pi*x[0])*cos(2*Pi*x[0])*cos(t)*sin(t)*sin(t) 
      - cos(2*Pi*x[0])*sin(t);
    source[1] = cos(2*Pi*x[0])*cos(t) 
      - 4*Pi*cos(2*Pi*x[0])*cos(2*Pi*x[0])*cos(2*Pi*x[0])*cos(t)*cos(t)*sin(t);
    source[2] = 4*Pi*cos(2*Pi*x[0])*cos(2*Pi*x[0])*cos(t)*sin(t);
  }


  // Derived quantities:
  //==================================================

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
