/*
  description of file goes here
*/

#include "generic.h"

#include "../../my_general_header.h"
#include "../../implicit_llg_problem.h"
#include "../../midpoint_method.h"
#include "../../magnetics_helpers.h"

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
  MyCliArgs args;
  args.parse(argc, argv);

  // Create problem
  ImplicitLLGProblem problem;

  // Assign timestepper and mesh from input arguments.
  problem.add_time_stepper_pt(args.time_stepper_pt);
  problem.set_bulk_mesh_pt(args.mesh_pt);
  problem.bulk_mesh_pt()->setup_boundary_element_info();

  // Set magnetic parameters
  double gilbert_damping = 0.5, hk = 0.0; // (normalised hk)
  problem.mag_parameters_pt()->set_simple_llg_parameters();
  problem.mag_parameters_pt()->gilbert_damping() = gilbert_damping;
  // problem.mag_parameters_pt()->exchange_constant() = 0.0;
  problem.mag_parameters_pt()->magnetostatic_debug_coeff() = 0.0;
  problem.mag_parameters_pt()->k1() = hk * mag_parameters::mu0/2;

  // Set applied field
  problem.applied_field_fct_pt() = args.h_app_fct_pt;

  // Set exact solution if we have one
  if((args.h_app_name == "minus_z") && (args.initial_m_name == "z"))
    {
      problem.Compare_with_mallinson = true;
    }

  // Finished setup, now we can build the problem
  problem.build();

  // Initialise problem and output
  problem.initialise_dt(args.dt);
  problem.set_initial_condition(args.initial_m_fct_pt);

  problem.Doc_info.set_directory(args.outdir);
  problem.Doc_info.Args_pt = &args;
  problem.Doc_info.output_jacobian = args.output_jacobian;

  problem.initial_doc();

  // All ready: step until completion
  double dt = args.dt;
  while(problem.time() < args.tmax)
    {
      std::cout << "step number = " << problem.Doc_info.number()
                << ", time = " << problem.time()
                << ", dt = " << dt
                << ", |m| error = " << 1 - problem.mean_nodal_magnetisation_length()
                << std::endl;

      // The Newton step itself, adaptive if requested
      if(args.adaptive_flag())
        {
          dt = problem.adaptive_unsteady_newton_solve(dt, args.tol);
        }
      else
        {
          problem.unsteady_newton_solve(dt);
        }

      // Output
      problem.doc_solution();
    }


  problem.final_doc();

  // Shut down oomph-lib's MPI
  MPI_Helpers::finalize();

  return 0;
}
