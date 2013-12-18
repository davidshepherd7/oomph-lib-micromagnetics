/*
  description of file goes here
*/

#include "generic.h"

#include "../../my_general_header.h"
#include "../../semi_implicit_problem.h"

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

  // Create objects
  // ============================================================

  // Read and process command line arguments
  SemiImplicitMMArgs args;
  args.parse(argc, argv);

  // Get a pointer to the problem class
  DecoupledLLGProblem* problem_pt = args.problem_pt;

  // Assign time steppers. Only one of these two will be a real
  // timestepper, the other will be null or a dummy instead.
  problem_pt->add_time_stepper_pt(args.time_stepper_pt);
  problem_pt->set_explicit_time_stepper_pt(args.explicit_time_stepper_pt);

  if(args.adaptive_flag())
    {
      problem_pt->Use_time_adaptive_newton = true;
    }

  // Assign other values
  problem_pt->llg_sub_problem_pt()->applied_field_fct_pt() = args.h_app_fct_pt;
  problem_pt->llg_sub_problem_pt()->linear_solver_pt() = args.solver_pt;
  problem_pt->llg_sub_problem_pt()->renormalise_each_time_step() = args.renormalise_flag();
  problem_pt->llg_sub_problem_pt()->set_mag_parameters_pt(args.magnetic_parameters_pt);
  problem_pt->llg_sub_problem_pt()->newton_solver_tolerance() = args.newton_tol;
  problem_pt->llg_sub_problem_pt()->Use_fd_jacobian = args.use_fd_jacobian;

  problem_pt->Doc_info.copy_args_string(&args);

  // Create and set phi_1 sub problem
  GenericPoissonProblem phi_1_problem;
  phi_1_problem.set_flux_mesh_factory(args.phi_1_flux_mesh_factory_fct_pt);
  phi_1_problem.newton_solver_tolerance() = args.newton_tol;
  problem_pt->set_phi_1_problem_pt(&phi_1_problem);


  // Create and set phi sub problem
  GenericPoissonProblem phi_problem;
  phi_problem.newton_solver_tolerance() = args.newton_tol;
  problem_pt->set_phi_problem_pt(&phi_problem);

  // Create and set the BEM handler
  BoundaryElementHandler bem_handler;
  problem_pt->bem_handler_pt() = &bem_handler;
  problem_pt->bem_handler_pt()->Bem_element_factory = args.bem_element_factory_fct_pt;
  problem_pt->bem_handler_pt()->Use_numerical_integration = args.use_numerical_integration;


  // Do the rest of the set up
  // ============================================================

  // Set exact solution if we have one
  if((args.h_app_name == "minus_z")
     && (args.initial_m_name == "z")
     && problem_pt->mag_parameters_pt()->gilbert_damping() != 0.0)
    {
      problem_pt->Compare_with_mallinson = true;
    }

  // Customise details of the output
  problem_pt->Doc_info.set_directory(args.outdir);
  problem_pt->Doc_info.output_jacobian = args.output_jacobian;

  // Finished customising the problem, now we can build it.
  problem_pt->build(args.llg_mesh_pts,
                args.phi_mesh_pts,
                args.phi_1_mesh_pts);


  // Initialise problem and output initial conditions
  // ============================================================
  problem_pt->initialise_dt(args.dt);
  problem_pt->set_initial_condition(args.initial_condition_fpt);
  problem_pt->set_doc_times(args.doc_times);

  problem_pt->initial_doc();


  // All ready: step until completion
  // ============================================================
  double dt = args.dt;
  unsigned time_step_number = 0;
  while(problem_pt->time() < args.tmax)
    {
      std::cout
        << std::endl
        << std::endl
        << "Time step " << time_step_number << std::endl
        << "==========================" << std::endl
        << "time = " << problem_pt->time()
        << ", dt = " << dt
        << ", |m| error = " << 1 - problem_pt->mean_nodal_magnetisation_length()
        << ", energy = " << problem_pt->micromagnetic_energy() << std::endl
        << std::setprecision(8)
        << ", previous step effective damping = " << problem_pt->Effective_damping_constant
        << ", alt damping calc = " << problem_pt->Alt_eff_damp
        << std::setprecision(4)
        << std::endl;

      // Take a step
      dt = problem_pt->do_step(dt, args.tol);

      // Output
      problem_pt->doc_solution();

      time_step_number++;
    }


  problem_pt->final_doc();

  // Shut down oomph-lib's MPI
  MPI_Helpers::finalize();

  return 0;
}
