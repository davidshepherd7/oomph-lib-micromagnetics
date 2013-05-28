
// Floating point debugging
#include <fenv.h>

#include "generic.h"
#include "meshes/tetgen_mesh.h"

#include "../../semi_implicit_problem.h"
#include "../../my_assert.h"
#include "../../boundary_element_handler.h"
#include "../../generic_poisson_problem.h"
#include "../../micromagnetics_boundary_element.h"

using namespace oomph;

namespace Inputs
{
  Vector<double> exact_M(const double& t, const Vector<double> &x)
  {
    Vector<double> M(3,0.0);
    M[0] = 1;
    return M;
  }
}

int main(int argc, char *argv[])
{
  // Start MPI if necessary
  MPI_Helpers::init(argc,argv);

  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  char *fake_argv[] = {"*dummy*", "-ref", "2", "-mesh", "ut_sphere", "-happ", "zero"};
  int fake_argc = 7;

  // Read and process the fake command line arguments
  SemiImplicitMMArgs fake_args;
  fake_args.parse(fake_argc, fake_argv);

  // Create main semi implicit problem
  SemiImplicitHybridMicromagneticsProblem problem;
  problem.add_time_stepper_pt(fake_args.time_stepper_pt);
  problem.llg_sub_problem_pt()->set_bulk_mesh_pt(fake_args.llg_mesh_pt);
  problem.llg_sub_problem_pt()->applied_field_fct_pt() = fake_args.h_app_fct_pt;
  problem.Doc_info.Args_pt = &fake_args;

  // Create and set phi_1 sub problem
  GenericPoissonProblem phi_1_problem;
  phi_1_problem.set_bulk_mesh(fake_args.phi_1_mesh_pt);
  phi_1_problem.set_flux_mesh_factory(fake_args.phi_1_flux_mesh_factory_fct_pt);
  problem.set_phi_1_problem_pt(&phi_1_problem);

  // Create and set phi sub problem
  GenericPoissonProblem phi_problem;
  phi_problem.set_bulk_mesh(fake_args.phi_mesh_pt);
  problem.set_phi_problem_pt(&phi_problem);

  // Create and set the BEM handler
  BoundaryElementHandler bem_handler;
  problem.bem_handler_pt() = &bem_handler;
  problem.bem_handler_pt()->Bem_element_factory = fake_args.bem_element_factory_fct_pt;

  // Set up the magnetic parameters
  problem.set_mag_parameters_pt(Factories::magnetic_parameters_factory("simple-llg"));

  // Customise details of the output
  problem.Doc_info.set_directory(fake_args.outdir);
  problem.Doc_info.output_jacobian = fake_args.output_jacobian;

  // Finished customising the problem, now we can build it.
  problem.build();


  // Initialise problem and output initial conditions
  // ============================================================
  problem.initialise_dt(1e-6);
  problem.set_initial_condition(Inputs::exact_M);

  problem.initial_doc();


  // Solve it
  // ============================================================
  problem.magnetostatics_solve();

  // Output results
  problem.doc_solution();

  // Compute field
  Vector<double> average_field = problem.average_magnetostatic_field();
  double rel_err = std::abs((average_field[0] - (-1.0/3.0))/average_field[0]);
  double abs_err = std::abs(average_field[0] - (-1.0/3.0));

  std::cout << average_field << std::endl;

  // Check close/convergence to analytical value
  std::cout
    << "Field values are [" << average_field[0] << ","
    << average_field[1] << "," << average_field[2] << "]\n\n"

    << "Error compared to an exact sphere is:\n"
    << "abs error = " << abs_err << "\n"
    << "rel error = " << rel_err << std::endl;

  // Return error status. We don't expect it to be very accurate because we
  // use a pretty bad approximation to a sphere for speed reasons.
  if(rel_err > 0.25)
    {
      std::cout << std::endl
                << "Relative error is too large!" << std::endl;
      return 1;
    }

  // Also check that y,z components are *all* close to zero (not average)



  // Shut down oomph-lib's MPI
  MPI_Helpers::finalize();

  // Also check that std-dev of x is not too high?
  return 0;
}
