

// Floating point debugging
#include <fenv.h>

#include<iomanip>

#include "../../midpoint_method.h"
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


int test_timestepper(TimeStepper* ts_pt,
                  double tmax, double tol, double dt,
                  const std::string& label, double max_global_error)
{
  // Build problem
  UnsteadyHeatProblem<QUnsteadyHeatElement<2,3> >
    problem(&UnsteadyHeatStepSolution::get_source,
            &UnsteadyHeatStepSolution::get_exact_u,
            ts_pt);

  // Nice output parameters
  // oomph_info << std::scientific << std::setprecision(2);
  // oomph_info.precision(8);

  // // Fiddle with scaling parameters so midpoint works nicely ??ds
  // problem.DTSF_max_increase = 4.0;
  // problem.DTSF_min_decrease = 0.75;

  // Setup labels for output
  DocInfo doc_info("results");

  // Open a trace file
  std::ofstream trace_file((doc_info.directory() + "/trace.dat").c_str());
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
      std::cout << "Timestep " << doc_info.number() << std::endl;
      std::cout << "dt = " << dt << std::endl;

      // Take timestep
      double dt_next = problem.adaptive_unsteady_newton_solve(dt, tol);
      dt = dt_next;

      std::cout << "Global error norm is: "
                << problem.get_error_norm() << std::endl;
      if(problem.get_error_norm() > max_global_error)
        {
          std::string error_msg = "Global error is too large! Test failed";
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      //Output solution
      problem.doc_solution(doc_info,trace_file);

      //Increment counter for solutions
      doc_info.number()++;
    }

  // Close trace file
  trace_file.close();

  return 0;

} // end of test_timestepper


int main(int argc, char *argv[])
{
  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  CommandLineArgs::setup(argc,argv);

  double tmax = 1.0;
  double tol = 1e-3;
  double dt = 1e-6;

  double max_global_error_norm = 0.02;

  test_timestepper(new MidpointMethod(true, 2, 0.1),
                   tmax, tol, dt, "_midpoint", max_global_error_norm);
  // test_timestepper(new BDF2(true), tmax, tol, dt, "_bdf2", max_global_error_norm);


  return 0;
}
