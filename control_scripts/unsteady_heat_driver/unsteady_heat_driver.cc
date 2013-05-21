


#include "generic.h"

#include "../../my_general_header.h"
#include "../../unsteady_heat_problem.h"


// Mesh
#include "meshes/rectangular_quadmesh.h"


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
  UnsteadyHeatArgs args;
  args.parse(argc, argv);

  // Create problem
  UnsteadyHeatProblem problem;

  // Assign timestepper and mesh from input arguments.
  problem.add_time_stepper_pt(args.time_stepper_pt);
  problem.mesh_pt() = new RectangularQuadMesh<QUnsteadyHeatElement<2,3> >(10, 10, 1.0,1.0,
                                                                          args.time_stepper_pt);
  problem.mesh_pt()->setup_boundary_element_info();
  problem.linear_solver_pt() = args.solver_pt;

  // Set source function
  problem.Source_fct_pt = args.source_fct_pt;

  // Set exact solution
  problem.Exact_solution_fct_pt = args.exact_fct_pt;

  problem.build();

  // Initialise problem and output
  problem.initialise_dt(args.dt);
  problem.set_initial_condition(args.exact_fct_pt);

  problem.Doc_info.set_directory(args.outdir);
  problem.Doc_info.Args_pt = &args;
  problem.Doc_info.output_jacobian = args.output_jacobian;

  problem.initial_doc();

  // All ready: step until completion
  double dt = args.dt;
  while(problem.time() < args.tmax)
    {
      std::cout
        << std::endl
        << std::endl
        << "Time step " << problem.Doc_info.number() << std::endl
        << "==========================" << std::endl
        << "time = " << problem.time()
        << ", dt = " << dt
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
