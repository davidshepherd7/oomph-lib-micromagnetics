// parameters for

#include "math.h"
using namespace std;
using namespace MathematicalConstants;

namespace OneDMicromagSetup
{

  using namespace OneDMicromagSetup;

  // Prototypes
  // ==========================================================
  void exact_solution(const double& t, const Vector<double>& x,
		      Vector<double>& exact);


  // Constants (can be redefined in individual parameter files)
  //===========================================================
  double alpha = 0.7;   // Gibert damping constant
  double gamma = 0.221E-8;   // Electromagnetic ratio

  // Overloadable functions
  //===========================================================

  // The coefficient of the damping term of the Landau-Lifschitz-Gilbert equation
  double llg_damping_coeff(const double& t, const Vector<double>& x)    
  {//??ds ifdef error checking then check ms is nonzero?
    return 1;}

  // The coefficient of the precession term of the Landau-Lifschitz-Gilbert equation
  double llg_precession_coeff(const double& t, const Vector<double>& x)
  {return 1;}

  double exchange_coeff(const double& t, const Vector<double>& x)
  {return 1;}

  // Crystalline anisotropy field - set to zero if unused
  void cryst_anis_field(const double& t, const Vector<double>& x, const Vector<double>& m, Vector<double>& H_cryst_anis)
  {
    H_cryst_anis[0] = 0.0;
    H_cryst_anis[1] = 0.0;
    H_cryst_anis[2] = 0.0;
  }

  // Applied field - set to zero if unused
  void applied_field(const double& t, const Vector<double>& x, Vector<double>& H_applied)
  {
    H_applied[0] = 0.0;
    H_applied[1] = 0.0;
    H_applied[2] = 0.0;
  }
  
  // Poisson source, set to zero unless needed for testing
  void source_function(const double& t, const Vector<double>& x, double& source)
  {
    source = 0.0;
  }

  // Non-overloadable functions 
  // ??ds probably never change, could add into micromagnetic elements header
  // but I need sat_mag for llg_coeff calcs...
  //===========================================================
  
  // Get the saturisation magnetisation at position x
  double sat_mag(const double& t, const Vector<double>& x)
  {
    Vector<double> exact(4,0.0);
    exact_solution(t,x,exact);
    return exact[1]*exact[1] + exact[2]*exact[2] + exact[3]*exact[3];
  }

  //??ds "legacy" functions - try to not use
  //==============================================================
  double boundary_phi(const double& t, const Vector<double>& x)
  {
    Vector<double> exact(4,0.0);
    exact_solution(t,x,exact);
    return exact[0];
  }

  void initial_M(const double& t, const Vector<double>& x, Vector<double>& M)
  {
    Vector<double> exact(4,0.0);
    exact_solution(t,x,exact);
    for(unsigned i=0; i<3; i++){M[i] = exact[i+1];}
  }

  void boundary_m(const double& t, const Vector<double>& x, Vector<double>& M)
  {
    Vector<double> exact(4,0.0);
    exact_solution(t,x,exact);
    for(unsigned i=0; i<3; i++){M[i] = exact[i+1];}
  }

  double exact_phi_solution(const double& t, const Vector<double>& x)
  {
    Vector<double> exact(4,0.0);
    exact_solution(t,x,exact);
    return exact[0];
  }

  void exact_M_solution(const double& t, const Vector<double>& x,
			Vector<double>& M)
  {
    Vector<double> exact(4,0.0);
    exact_solution(t,x,exact);
    for(unsigned i=0; i<3; i++){M[i] = exact[i+1];}
  }
  

}; // End of namespace
