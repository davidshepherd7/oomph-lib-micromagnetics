/*
  description of file goes here
*/

#include "generic.h"
#include "micromag.h"

// Floating point error checks
#include <fenv.h>

using namespace oomph;
using namespace MathematicalConstants;
using namespace StringConversion;



int main(int argc, char *argv[])
{
  // Start MPI if necessary
  MPI_Helpers::init(argc,argv);

#ifdef PARANOID
  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
#endif

  // Read and process command line arguments
  LLGArgs args;
  args.parse(argc, argv);

  // Create problem
  LLGProblem problem;

  // Tell it if we want renormalisation or not
  problem.renormalise_each_time_step() = args.renormalise_flag();

  // Assign time stepper
  if(args.explicit_flag() && args.use_implicit_ms)
    {
      std::string err = "Cannot do explicit but fully coupled time stepping ";
      err += "when magnetostatics is included, use the decoupled driver instead";
      throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                          OOMPH_CURRENT_FUNCTION);
    }
  else
    {
      // Assign time steppers. Only one of these two will be a real
      // timestepper, the other will be null or a dummy instead.
      problem.add_time_stepper_pt(args.time_stepper_pt);
      problem.set_explicit_time_stepper_pt(args.explicit_time_stepper_pt);
    }

  if(args.use_implicit_ms && args.time_stepper_name == "midpoint")
    {
      std::string err = "Fully implicit is not set up for use with standard midpoint";
      err += " use midpoint-bdf instead.";
      throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                          OOMPH_CURRENT_FUNCTION);
    }

  problem.Use_time_adaptive_newton = args.adaptive_flag();
  problem.linear_solver_pt() = args.solver_pt;
  problem.set_mag_parameters_pt(args.magnetic_parameters_pt);
  problem.newton_solver_tolerance() = args.newton_tol;
  problem.Use_fd_jacobian = args.use_fd_jacobian;
  problem.Pin_boundary_m = args.pin_boundary_m;
  problem.Use_implicit_ms = args.use_implicit_ms;
  problem.Residual_calculator_pt = args.residual_calculator_pt;


  // Set applied field
  problem.applied_field_fct_pt() = args.h_app_fct_pt;

  // Set exact solution if we have one
  if((args.h_app_name == "minus_z")
     && (args.initial_m_name == "z")
     && (problem.mag_parameters_pt()->gilbert_damping() != 0.0)
     && (args.use_implicit_ms == false))
    {
      problem.Compare_with_mallinson = true;
    }

  // Finished setup, now we can build the problem
  problem.build(args.mesh_pts);

  // Initialise problem and output
  problem.initialise_dt(args.dt);
  problem.set_initial_condition(args.initial_m_fct_pt);

  problem.Doc_info.set_directory(args.outdir);
  problem.Doc_info.Args_pt = &args;
  problem.Doc_info.output_jacobian = args.output_jacobian;
  problem.set_doc_times(args.doc_times);

  problem.initial_doc();


  // All ready: step until completion
  double dt = args.dt;
  unsigned time_step_number = 0;
  while(problem.time() < args.tmax)
    {
      time_step_number++;

      std::cout
        << std::endl
        << std::endl
        << "Time step " << time_step_number << std::endl
        << "==========================" << std::endl
        << "time = " << problem.time()
        << ", dt = " << dt
        << ", |m| error = " << 1 - problem.mean_nodal_magnetisation_length()
        << ", energy = " << problem.micromagnetic_energy() << std::endl
        << ", previous step effective damping = " << problem.Effective_damping_constant
        << ", alt damping calc = " << problem.Alt_eff_damp
        << std::endl;

      // Do the newton solve (different ones depending flags set)
      dt = problem.smart_newton_solve(dt, args.tol);

      // Output
      problem.doc_solution();
    }


  problem.final_doc();

  // Shut down oomph-lib's MPI
  MPI_Helpers::finalize();

  return 0;
}
