/*
  description of file goes here
*/

#include "generic.h"

#include "../../my_general_header.h"
#include "../../semi_implicit_problem.h"

#include "micromag.h"

// Floating point error checks
#include <fenv.h>

// Signal handling
#include <signal.h>

using namespace oomph;
using namespace MathematicalConstants;
using namespace StringConversion;


/// Exception for interrupt signals
class MyInterrupt : public std::exception
{
  public:
  MyInterrupt(int s) : S(s) {}
  int S;
};

/// Handle a signal by raising an exception
void sig_to_exception_handler(int s)
{
  throw MyInterrupt(s);
}



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

  // Create main semi implicit problem
  SemiImplicitHybridMicromagneticsProblem problem;
  problem.add_time_stepper_pt(args.time_stepper_pt);
  problem.llg_sub_problem_pt()->applied_field_fct_pt() = args.h_app_fct_pt;
  problem.llg_sub_problem_pt()->linear_solver_pt() = args.solver_pt;
  problem.llg_sub_problem_pt()->renormalise_each_time_step() = args.renormalise_flag();
  problem.llg_sub_problem_pt()->set_mag_parameters_pt(args.magnetic_parameters_pt);
  problem.llg_sub_problem_pt()->newton_solver_tolerance() = args.newton_tol;
  problem.llg_sub_problem_pt()->Use_fd_jacobian = args.use_fd_jacobian;
  problem.llg_sub_problem_pt()->Residual_calculator_pt = args.residual_calculator_pt;

  problem.Doc_info.Args_pt = &args;

  // Create and set phi_1 sub problem
  GenericPoissonProblem phi_1_problem;
  phi_1_problem.set_flux_mesh_factory(args.phi_1_flux_mesh_factory_fct_pt);
  phi_1_problem.newton_solver_tolerance() = args.newton_tol;
  problem.set_phi_1_problem_pt(&phi_1_problem);


  // Create and set phi sub problem
  GenericPoissonProblem phi_problem;
  phi_problem.newton_solver_tolerance() = args.newton_tol;
  problem.set_phi_problem_pt(&phi_problem);

  // Create and set the BEM handler
  BoundaryElementHandler bem_handler;
  problem.bem_handler_pt() = &bem_handler;
  problem.bem_handler_pt()->Bem_element_factory = args.bem_element_factory_fct_pt;
  problem.bem_handler_pt()->Use_numerical_integration = args.use_numerical_integration;


  // Do the rest of the set up
  // ============================================================

  // Set exact solution if we have one
  if((args.h_app_name == "minus_z")
     && (args.initial_m_name == "z")
     && problem.mag_parameters_pt()->gilbert_damping() != 0.0)
    {
      problem.Compare_with_mallinson = true;
    }

  // Customise details of the output
  problem.Doc_info.set_directory(args.outdir);
  problem.Doc_info.output_jacobian = args.output_jacobian;

  // Finished customising the problem, now we can build it.
  problem.build(args.llg_mesh_pts,
                args.phi_mesh_pts,
                args.phi_1_mesh_pts);


  // Initialise problem and output initial conditions
  // ============================================================
  problem.initialise_dt(args.dt);
  problem.set_initial_condition(args.initial_m_fct_pt);
  problem.set_doc_times(args.doc_times);

  problem.initial_doc();


  // Set up signal handler to catch (interrupt) signals and throw an
  // exception instead. Use try to catch this exception so that we can
  // write the end of the xml .pvd file (which lets us view results in
  // paraview even if the run didn't finish).
  struct sigaction sigIntHandler;
  sigIntHandler.sa_handler = sig_to_exception_handler;
  sigemptyset(&sigIntHandler.sa_mask);
  sigIntHandler.sa_flags = 0;
  sigaction(SIGINT, &sigIntHandler, NULL);
  try
    {
      // All ready: step until completion
      // ============================================================
      double dt = args.dt;
      unsigned time_step_number = 0;
      while(problem.time() < args.tmax)
        {
          std::cout
            << std::endl
            << std::endl
            << "Time step " << time_step_number << std::endl
            << "==========================" << std::endl
            << "time = " << problem.time()
            << ", dt = " << dt
            << ", |m| error = " << 1 - problem.mean_nodal_magnetisation_length()
            << ", energy = " << problem.micromagnetic_energy() << std::endl
            << std::setprecision(8)
            << ", previous step effective damping = " << problem.Effective_damping_constant
            << ", alt damping calc = " << problem.Alt_eff_damp
            << std::setprecision(4)
            << std::endl;

          // Take a step (automatically adaptive if args.tol != 0.0)
          dt = problem.semi_implicit_step(dt, args.tol);

          // Output
          problem.doc_solution();

          time_step_number++;
        }
    }
  catch(MyInterrupt& e)
    {
      std::cout << "caught signal " << e.S << std::endl;
      problem.final_doc();
      exit(1);
    }

  problem.final_doc();

  // Shut down oomph-lib's MPI
  MPI_Helpers::finalize();

  return 0;
}
