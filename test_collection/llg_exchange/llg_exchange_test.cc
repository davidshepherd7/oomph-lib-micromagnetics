
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
  double dt = 0.1;
  double tmax = 0.1;

  unsigned nx = 10, ny = nx;
  double lx = 1.0, ly = lx;

  void initial_m(const double& t, const Vector<double> &x, Vector<double> &m)
  {
    m.assign(3,0.0);

    // m[0] = sin(x[0]*2*Pi);
    // m[1] = cos(x[1]*2*Pi);
    // m[2] = 1.0 - m[0] - m[1];

    m[0] = x[0]/2 + x[1]/2;
    m[1] = (1 - m[0]);
    m[2] = 0.0;

    VectorOps::normalise(m);
  }

  //??ds applied field?
  void applied_field(const double& t, const Vector<double> &x,
                     Vector<double> &h_app)
  {
    h_app.assign(3,0.0);
    h_app[0] = 0.01;
  }

}


int main()
{
  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  double dt = Inputs::dt;
  double tmax = Inputs::tmax;

  // Build problem
  ImplicitLLGProblem<TMicromagElement<2,2> >
    problem(Inputs::nx, Inputs::ny, Inputs::lx, Inputs::ly,
            Inputs::applied_field);


  // Fiddle with exchange strength.
  problem.mag_parameters_pt()->exchange_constant() = mag_parameters::mu0/5;
  // Complicated because of normalisation. To get |Hex| = 1 we have to get
  // the exchange constant to mu0. Here we want |Hex| = 0.2 so divide by
  // 5...
  problem.mag_parameters_pt()->gilbert_damping() = 0.1;

  //??ds better solver in the end?

  // Initialise
  problem.initialise_dt(dt);
  problem.set_initial_condition(Inputs::initial_m);

  // Set up output
  DocInfo doc_info;
  doc_info.set_directory("results");

  problem.doc_solution(doc_info);

  // Timestep
  unsigned nstep = int(tmax/dt);
  for(unsigned t=0; t < nstep; t++)
    {
      std::cout << "Timestep " << t << std::endl;

      // Step
      problem.unsteady_newton_solve(dt);

      // Output
      problem.doc_solution(doc_info);
      doc_info.number()++;
    }

  // // ??ds
  // // Check solution
  // double avg_m = problem.average_magnetisation();

  // error checks not done yet so pretend to fail

  return 244;
}
