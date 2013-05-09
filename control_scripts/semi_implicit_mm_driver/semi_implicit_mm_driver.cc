/*
  description of file goes here
*/

#include "generic.h"

#include "../../my_general_header.h"
#include "../../semi_implicit_problem.h"
#include "../../midpoint_method.h"
#include "../../magnetics_helpers.h"

#include "../../micromag.h"

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

  // Create main semi implicit problem
  SemiImplicitHybridMicromagneticsProblem problem;
  problem.add_time_stepper_pt(args.time_stepper_pt);
  problem.llg_sub_problem_pt()->set_bulk_mesh_pt(args.llg_mesh_pt);
  problem.llg_sub_problem_pt()->applied_field_fct_pt() = args.h_app_fct_pt;
  problem.llg_sub_problem_pt()->linear_solver_pt() = args.solver_pt;
  problem.Doc_info.Args_pt = &args;

  // Create and set phi_1 sub problem
  GenericPoissonProblem phi_1_problem;
  phi_1_problem.set_bulk_mesh(args.phi_1_mesh_pt);
  phi_1_problem.set_flux_mesh_factory(args.phi_1_flux_mesh_factory_fct_pt);
  // (CG w/ Ilu0 prec works well enough to hard code it).
  phi_1_problem.linear_solver_pt() = Factories::linear_solver_factory("cg-ilu");
  problem.set_phi_1_problem_pt(&phi_1_problem);

  // Create and set phi sub problem
  GenericPoissonProblem phi_problem;
  phi_problem.set_bulk_mesh(args.phi_mesh_pt);
  // (CG w/ Ilu0 prec works well enough to hard code it).
  phi_problem.linear_solver_pt() = Factories::linear_solver_factory("cg-ilu");
  problem.set_phi_problem_pt(&phi_problem);

  // Create and set the BEM handler
  BoundaryElementHandler bem_handler;
  problem.bem_handler_pt() = &bem_handler;
  problem.bem_handler_pt()->Bem_element_factory = args.bem_element_factory_fct_pt;


  // Do the rest of the set up
  // ============================================================

  // Set up the magnetic parameters
  problem.mag_parameters_pt()->set_simple_llg_parameters();
  // problem.mag_parameters_pt()->magnetostatic_debug_coeff() = 0;

  // Set exact solution if we have one
  if((args.h_app_name == "minus_z") && (args.initial_m_name == "z"))
    {
      problem.Compare_with_mallinson = true;
    }

  // Customise details of the output
  problem.Doc_info.set_directory(args.outdir);
  problem.Doc_info.output_jacobian = args.output_jacobian;

  // Finished customising the problem, now we can build it.
  problem.build();


  // Initialise problem and output initial conditions
  // ============================================================
  problem.initialise_dt(args.dt);
  problem.set_initial_condition(args.initial_m_fct_pt);

  problem.initial_doc();


  // All ready: step until completion
  // ============================================================
  double dt = args.dt;
  while(problem.time() < args.tmax)
    {
      std::cout << "step number = " << problem.Doc_info.number()
                << ", time = " << problem.time()
                << ", dt = " << dt
                << ", |m| error = " << 1 - problem.mean_nodal_magnetisation_length()
                << ", error norm = " << problem.llg_sub_problem_pt()->get_error_norm()
                << std::endl;

      // Take a step (adaptive if args.tol != 0.0)
      dt = problem.semi_implicit_step(dt, args.tol);

      // Output
      problem.doc_solution();
    }


  problem.final_doc();

  // Shut down oomph-lib's MPI
  MPI_Helpers::finalize();

  return 0;
}
