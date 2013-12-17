//Generic routines
#include "generic.h"

#include "micromag.h"

// Floating point error checks
#include <fenv.h>

using namespace oomph;


int main(int argc, char* argv[])
{

#ifdef PARANOID
  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
#endif

  // Parse args
  ODECliArgs args;
  args.parse(argc, argv);

  // Create problem
  ODEProblem problem;

  // Assign from args
  problem.linear_solver_pt() = args.solver_pt;
  problem.newton_solver_tolerance() = args.newton_tol;
  // problem.Use_fd_jacobian = args.use_fd_jacobian;
  problem.Use_time_adaptive_newton = args.adaptive_flag();

  // Assign time steppers. Only one of the timesteppers will be a real
  // timestepper, the other will be null or a dummy instead.
  problem.add_time_stepper_pt(args.time_stepper_pt);
  problem.set_explicit_time_stepper_pt(args.explicit_time_stepper_pt);

  problem.build(args.mesh_pts);

  // Initialise problem and output
  problem.initialise_dt(args.dt);
  problem.set_initial_condition();

  problem.Doc_info.set_directory(args.outdir);
  problem.Doc_info.copy_args_string(&args);
  problem.Doc_info.output_jacobian = args.output_jacobian;
  problem.set_doc_times(args.doc_times);
  problem.initial_doc();

  // All ready: step until completion
  double dt = args.dt;
  unsigned time_step_number = 0;
  while (problem.time() < args.tmax)
    {
      time_step_number++;

      std::cout
        << std::endl
        << std::endl
        << "Time step " << time_step_number << std::endl
        << "==========================" << std::endl
        << "time = " << problem.time()
        << ", dt = " << dt
        << std::endl;

      // Do the newton solve (different ones depending flags set)
      dt = problem.smart_newton_solve(dt, args.tol);

      // Output
      problem.doc_solution();
    }


} // end_of_main
