
// Floating point error checks
#include <fenv.h>

#include "generic.h"
#include "meshes/simple_rectangular_tri_mesh.h"

#include "../../implicit_llg_problem.h"


/*

  Idea: run llg on 2D mesh with only exchange and small applied field +
  very rough initial m. All m vectors should align to the small applied
  field.

*/

using namespace oomph;
using namespace MathematicalConstants;

namespace Inputs
{

  void initial_m(const double& t, const Vector<double> &x, Vector<double> &m)
  {
    m.assign(3,0.0);

    // m[0] = sin(x[0]*2*Pi);
    // m[1] = cos(x[1]*2*Pi);
    // m[2] = 1.0 - m[0] - m[1];

    m[0] = x[0]/4 + x[1]/4 + 1.0/4.0;
    m[1] = (1 - m[0]);
    m[2] = 0.0;

    VectorOps::normalise(m);
  }

  //??ds applied field?
  void applied_field(const double& t, const Vector<double> &x,
                     Vector<double> &h_app)
  {
    h_app.assign(3,0.0);
    h_app[0] = 2.0;
  }

}


int main()
{
  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  // Build problem
  ImplicitLLGProblem<TMicromagElement<2,2> >
    problem(10, 10, 1.0, 1.0, Inputs::applied_field, true);

  // Fiddle with exchange strength.
  problem.mag_parameters_pt()->exchange_constant() = mag_parameters::mu0/2;
  // Complicated because of normalisation. To get |Hex| = 1 we have to get
  // the exchange constant to mu0/2.

  // Make it nice and quick to solve (since this is a test).
  problem.mag_parameters_pt()->gilbert_damping() = 0.8;

  //??ds better solver in the end?

  // Initialise
  problem.set_initial_condition(Inputs::initial_m);
  double eps = 1e-4;
  double tmax = 7;
  double dt = 0.03; // initial suggestion for dt

  // Set up output
  DocInfo doc_info;
  doc_info.set_directory("results");

  // problem.doc_solution(doc_info);

  // Do one timestep with known dt to test intermediate results
  {
    std::cout << "step number = " << doc_info.number()
              << ", time = " << problem.time_pt()->time()
              << std::endl;

    problem.unsteady_newton_solve(dt);

    std::cout << " dt = " << dt << std::endl;

    // Output
    std::ofstream first_solution("validation/generated_solution_at_t0.03.dat");
    problem.mesh_pt()->output(first_solution,2);
    first_solution.close();
  }

  // Timestep to end
  while(problem.time_pt()->time() < tmax)
    {
      std::cout << "step number = " << doc_info.number()
                << ", time = " << problem.time_pt()->time()
                << std::endl;

      // Step (use recommended next time step from previous step).
      dt = problem.adaptive_unsteady_newton_solve(dt,eps);

      std::cout << " dt = " << dt << std::endl;

      // // Output
      // problem.doc_solution(doc_info);
    }


  // Check solution--should have converged to M = [1,0,0] by now.
  Vector<double> m;
  problem.mean_magnetisation(m); m[0] -= 1;

  double tol = 1e-3;
  for(unsigned j=0; j<3; j++)
    {
      if( std::abs(m[j]) > tol) return 2;
    }

  return 0;
}
