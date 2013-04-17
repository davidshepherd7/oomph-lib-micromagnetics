

// Floating point debugging
#include <fenv.h>

#include<iomanip>

#include "../../old_midpoint_method.h"
#include "generic.h"
#include "./unsteady_heat_example.h"
#include "../../my_assert.h"



using namespace oomph;
using namespace MathematicalConstants;

//======start_of_UnsteadyHeatStepSolution=============================
/// Namespace for an exact solution for UnsteadyHeat equation
//====================================================================
namespace UnsteadyHeatStepSolution
{

  /// Factor controlling the rate of change
  double Gamma=10.0;

  /// Wavenumber
  double K=3.0;

  /// Angle of bump
  double Phi=1.0;

  /// Exact solution as a Vector
  void get_exact_u(const double& time, const Vector<double>& x,
                   Vector<double>& u)
  {
    u.assign(1,0.0);

    double zeta = cos(Phi)*x[0] + sin(Phi)*x[1];
    u[0] = sin(K*zeta) * 0.5 * (1.0 + tanh(Gamma*cos(2.0*Pi*time)));
  }

  /// Source function to make it an exact solution
  void get_source(const double& time, const Vector<double>& x, double& source)
  {
    source=
      -0.5*sin(K*(cos(Phi)*x[0]+sin(Phi)*x[1]))*K*K*pow(cos(Phi),2.0)*(
                                                                       0.1E1+tanh(Gamma*cos(0.2E1*0.3141592653589793E1*time)))-
      0.5*sin(K*(cos(Phi)*x[0]+sin(Phi)*x[1]))*K*K*pow(sin(Phi),2.0)*
      (0.1E1+tanh(Gamma*cos(0.2E1*0.3141592653589793E1*time)))+
      0.1E1*sin(K*(cos(Phi)*x[0]+sin(Phi)*x[1]))*
      (1.0-pow(tanh(Gamma*cos(0.2E1*0.3141592653589793E1*time)),2.0))*
      Gamma*sin(0.2E1*0.3141592653589793E1*time)*0.3141592653589793E1;
  }

} // end of UnsteadyHeatStepSolution


void main_bdf(double tmax, double tol, double dt,
              const std::string& label)
{
  BDF<1> my_bdf1(true);
  BDF<2> my_bdf2(true);

  // Build problem
  UnsteadyHeatProblem<QUnsteadyHeatElement<2,3> >
    problem(&UnsteadyHeatStepSolution::get_source,
            &UnsteadyHeatStepSolution::get_exact_u,
            &my_bdf2);

  // Fiddle with scaling parameters so midpoint works nicely ??ds
  problem.DTSF_max_increase = 4.0;
  problem.DTSF_min_decrease = 0.75;

  // Setup labels for output
  DocInfo doc_info("results"+label);

  // Open a trace file
  ofstream trace_file((doc_info.directory() + "/trace.dat").c_str());
  trace_file << "\"time\", \"u FE\", \"u exact\", \"norm of error\", \"norm of solution\", \"dt\""
             << std::endl;

  // Initialise timestep -- also sets the weights for all timesteppers
  // in the problem.
  problem.initialise_dt(dt);
  problem.set_initial_condition();

  //Output initial condition
  problem.doc_solution(doc_info, trace_file);

  //Increment counter for solutions
  doc_info.number()++;

  // Timestepping loop
  while(problem.time() < tmax)
    {
      cout << "Timestep " << doc_info.number() << std::endl;
      std::cout << "dt = " << dt << std::endl;

      // Take timestep
      double dt_next = problem.adaptive_unsteady_newton_solve(dt, tol);
      dt = dt_next;

      //Output solution
      problem.doc_solution(doc_info,trace_file);

      //Increment counter for solutions
      doc_info.number()++;
    }

  // Close trace file
  trace_file.close();

} // end of main_bdf


 void main_midpoint(double tmax, double tol, double dt,
                    const std::string& label)
 {
   OldMidpointMethod adaptive_midpoint(true, 2);
   adaptive_midpoint.Fudge_factor = 0.1;

   // Build problem
   UnsteadyHeatProblem<QUnsteadyHeatElement<2,3> >
     problem(&UnsteadyHeatStepSolution::get_source,
             &UnsteadyHeatStepSolution::get_exact_u,
             &adaptive_midpoint);

   // Nice output parameters
   // oomph_info << std::scientific << std::setprecision(2);
   // oomph_info.precision(8);

   // Fiddle with scaling parameters so midpoint works nicely ??ds
   problem.DTSF_max_increase = 4.0;
   problem.DTSF_min_decrease = 0.75;

   // Setup labels for output
   DocInfo doc_info("results" + label);

   // Open a trace file
   ofstream trace_file((doc_info.directory() + "/trace.dat").c_str());
   trace_file << "\"time\", \"u FE\", \"u exact\", \"norm of error\", \"norm of solution\", \"dt\""
              << std::endl;

   // Initialise timestep -- also sets the weights for all timesteppers
   // in the problem.
   problem.initialise_dt(dt);
   problem.set_initial_condition();

   //Output initial condition
   problem.doc_solution(doc_info, trace_file);

   //Increment counter for solutions
   doc_info.number()++;

   // Timestepping loop
   while(problem.time() < tmax)
     {
       cout << "Timestep " << doc_info.number() << std::endl;
       std::cout << "dt = " << dt << std::endl;

       // Take timestep
       double dt_next = problem.adaptive_unsteady_newton_solve(dt, tol);
       dt = dt_next;

       //Output solution
       problem.doc_solution(doc_info,trace_file);

       //Increment counter for solutions
       doc_info.number()++;
     }

   // Close trace file
   trace_file.close();

 } // end of main_midpoint




int main(int argc, char *argv[])
{
  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  CommandLineArgs::setup(argc,argv);

  double tmax = 1.0;
  double tol = 1e-3;
  double dt = 1e-6;

  // main_bdf(tmax, tol, dt, "_bdf2");
  main_midpoint(tmax, tol, dt, "_midpoint");

  return 0;
}
