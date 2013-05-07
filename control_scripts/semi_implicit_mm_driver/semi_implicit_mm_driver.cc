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

MicromagBEMElementEquations* temp_ele_factory(FiniteElement* const ele,
                                              const int& face)
{
  return new QMicromagBEMElement<2,2>(ele, face);
}

int main(int argc, char *argv[])
{
  // Start MPI if necessary
  MPI_Helpers::init(argc,argv);

#ifdef PARANOID
  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
#endif

  // Read and process command line arguments
  SemiImplicitMMArgs args;
  args.parse(argc, argv);

  // Create problem
  SemiImplicitHybridMicromagneticsProblem problem;

  problem.add_time_stepper_pt(args.time_stepper_pt);
  problem.llg_sub_problem_pt()->set_bulk_mesh_pt(args.llg_mesh_pt);

  GenericPoissonProblem* phi_1_problem_pt = new GenericPoissonProblem;
  phi_1_problem_pt->set_bulk_mesh(args.phi_1_mesh_pt);
  phi_1_problem_pt->Flux_mesh_factory =
    &Factories::surface_mesh_factory<QMagnetostaticFieldFluxElement<2,2> >;
  problem.set_phi_1_problem_pt(phi_1_problem_pt);

  GenericPoissonProblem* phi_problem_pt = new GenericPoissonProblem;
  phi_problem_pt->set_bulk_mesh(args.phi_mesh_pt);
  problem.set_phi_problem_pt(phi_problem_pt);

  problem.bem_handler_pt() = new BoundaryElementHandler;
  problem.bem_handler_pt()->Bem_element_factory = &temp_ele_factory;


  // Do the rest (mag parameters, phi etc.)
  // ============================================================

  // Set magnetic parameters
  double gilbert_damping = 0.5; // (normalised hk)
  problem.mag_parameters_pt()->set_simple_llg_parameters();
  problem.mag_parameters_pt()->gilbert_damping() = gilbert_damping;
  problem.mag_parameters_pt()->magnetostatic_debug_coeff() = 0;


  // Set applied field
  problem.llg_sub_problem_pt()->applied_field_fct_pt() = args.h_app_fct_pt;

  // Set exact solution if we have one
  if((args.h_app_name == "minus_z") && (args.initial_m_name == "z"))
    {
      problem.Compare_with_mallinson = true;
    }

  // Finished inputing stuff, now we can build the problem
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
