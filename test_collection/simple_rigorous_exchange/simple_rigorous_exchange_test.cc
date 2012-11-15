/*

  If only exchange field is turned on this initial magnetisation should not
  change over time. Checks that exchange has been implemented sort of right
  (at least for 1D variations).


  2d
  ================================================================

  There is a jump in error after the first Newton step (still only
  ~1e-6). This could indicate slightly wrong residuals or it could be the
  natural limitations due to rounding errors.

  Try tightening newton tolerence to see? Nope, no effect :(

  Strongly reducting eps seems to reduce this error! Note we need eps much
  smaller than we want the real error to be indicating that our error
  estimator is rubbish! Quite possibly this is due to not using mid-point
  method.


  3d
  ================================================================

  Tests with 3d are pretty hard: we can't currently run fast at a decent
  mesh refinement so the test will likely not pass. However we can check if
  the Newton method converges and the error at least doesn't go entirely
  crazy, which is sort of useful...

  3d passes for a = 0.1, nx = 10, x,y and z variations.

*/

#include <fenv.h>

#include "generic.h"
#include "../../implicit_llg_problem.h"

#include "meshes/simple_cubic_tet_mesh.h"

using namespace oomph;


namespace Inputs
{

  double a = 1.0;

  void steady_state_initial_m_x_variations
  (const double& t, const Vector<double> &x, Vector<double> &m)
  {
    m.assign(3,0.0);

    m[0] = std::sin(a*x[0]);
    m[1] = std::cos(a*x[0]);
    m[2] = 0.0;

    // should be normalised already
  }


  void steady_state_initial_m_y_variations
  (const double& t, const Vector<double> &x, Vector<double> &m)
  {
    m.assign(3,0.0);

    m[0] = std::sin(a*x[1]);
    m[1] = std::cos(a*x[1]);
    m[2] = 0.0;

    // should be normalised already
  }

  void steady_state_initial_m_z_variations
  (const double& t, const Vector<double> &x, Vector<double> &m)
  {
    m.assign(3,0.0);

    m[0] = std::sin(a*x[2]);
    m[1] = std::cos(a*x[2]);
    m[2] = 0.0;

    // should be normalised already
  }


  // Turn off field
  void no_applied_field(const double& t, const Vector<double> &x,
                        Vector<double> &h_app)
  {
    h_app.assign(3,0.0);
  }

  // Function pointer for initial magnetisation.
  typedef void (*InitialMFctPt)(const double& t, const Vector<double> &x,
                                Vector<double> &m);

}

using namespace Inputs;


int main(int argc, char *argv[])
{
  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  // General parameters
  // ============================================================

  // Try with a few values of a
  Vector<double> as;
  // as.push_back(2.1);
  as.push_back(0.1);
  as.push_back(1.0);
  as.push_back(MathematicalConstants::Pi);
  as.push_back(std::sqrt(2));

  // Try with a few intial conditions
  Vector<Inputs::InitialMFctPt> initial_conds;
  // initial_conds.push_back(&Inputs::steady_state_initial_m_y_variations);
  // initial_conds.push_back(&Inputs::steady_state_initial_m_x_variations);

  double eps = 1e-10; // Need a really small eps to stay in steady state,
                      // possibly because we need to use mid-point
                      // really...
  double tmax = 10;
  double  dt = 0.03;

  // Set up output
  DocInfo doc_info;
  doc_info.set_directory("results");


  // Decide which tests to run
  bool run_2d = true;
  bool run_3d = true;


  // A 2D problem
  // ============================================================

  if(run_2d)
    {
      // Build problem
      ImplicitLLGProblem<TMicromagElement<2,2> >
        problem(30, 30, 1.0, 1.0, Inputs::no_applied_field, true);

      problem.mag_parameters_pt()->set_simple_llg_parameters();
      problem.renormalise_each_time_step() = true;

      // Mangle the mesh a little (so that we have dmdn non-zero on some
      // boundaries).
      for(unsigned nd=0, nnode=problem.bulk_mesh_pt()->nnode(); nd<nnode; nd++)
        {
          Node* nd_pt = problem.bulk_mesh_pt()->node_pt(nd);
          nd_pt->x(0) += nd_pt->x(1);
        }

      // //??ds check with a tighter newton tol:
      // problem.newton_solver_tolerance() = 1e-10;

      for(unsigned i_ic=0; i_ic<initial_conds.size(); i_ic++)
        {
          for(unsigned i=0; i<as.size(); i++)
            {

              // Change to a new a parameter
              Inputs::a = as[i];

              std::cout << std::endl << std::endl;
              std::cout <<  "====================" << std::endl;
              std::cout <<  "New problem, a = " << Inputs::a << std::endl;
              std::cout <<  "====================" << std::endl;

              // Reset problem
              problem.set_initial_condition(initial_conds[i_ic]);
              problem.doc_solution(doc_info);
              problem.time() = 0.0;
              dt = 0.03;

              // Allow double the initial discretisation error in final error.
              double errtol = problem.compare_m_with_function(initial_conds[i_ic])
                * 2;

              std::cout << std::endl
                        << "error = " << problem.compare_m_with_function(initial_conds[i_ic])
                        << std::endl;

              // Timestep to end
              while(problem.time_pt()->time() < tmax)
                {
                  std::cout << "step number = " << doc_info.number()
                            << ", time = " << problem.time_pt()->time()
                            << std::endl;

                  // Step (use recommended next time step from previous step).
                  dt = problem.adaptive_unsteady_newton_solve(dt,eps);

                  std::cout << " dt = " << dt << std::endl;

                  // Output
                  problem.doc_solution(doc_info);

                  std::cout << "error = " << problem.compare_m_with_function(initial_conds[i_ic])
                            << std::endl;
                }

              // Check that m hasn't changed, error if it has changed more than tol
              if((problem.compare_m_with_function(initial_conds[i_ic]))
                 > errtol)
                {
                  return i_ic*as.size() + i+1;
                }
              else
                {
                  std::cout << "Ok for a = " << Inputs::a << std::endl;
                }

            }
        }

    }

  // Try a 3d version
  // ============================================================
  if(run_3d)
    {

      // Build mesh
      BDF<2> bdf2(true);
      unsigned nx = 10, ny = nx, nz = nx;
      double lx = 1.0, ly = lx, lz = lx;
      SimpleCubicTetMesh<TMicromagElement<3,2> > threedmesh(nx,ny,nz,lx,ly,lz,&bdf2);
      threedmesh.setup_boundary_element_info();


      // Build problem
      ImplicitLLGProblem<TMicromagElement<3,2> > threedproblem;
      threedproblem.bulk_mesh_pt() = &threedmesh;
      threedproblem.add_time_stepper_pt(&bdf2);
      threedproblem.applied_field_fct_pt() = &Inputs::no_applied_field;
      threedproblem.mag_parameters_pt()->set_simple_llg_parameters();
      threedproblem.renormalise_each_time_step() = true;

      threedproblem.build();

      // Add a z varying inital condition to the list
      initial_conds.push_back(Inputs::steady_state_initial_m_z_variations);

      // Make a new, nicer list of a values because we can't refine well in z
      // direction.
      Vector<double> threedas;
      threedas.push_back(0.1);

      for(unsigned i_ic=0; i_ic<initial_conds.size(); i_ic++)
        {
          for(unsigned i=0; i<threedas.size(); i++)
            {

              // Change to a new a parameter
              Inputs::a = threedas[i];

              std::cout << std::endl << std::endl;
              std::cout <<  "====================" << std::endl;
              std::cout <<  "New threedproblem, a = " << Inputs::a << std::endl;
              std::cout <<  "====================" << std::endl;

              // Reset threedproblem
              threedproblem.set_initial_condition(initial_conds[i_ic]);
              threedproblem.doc_solution(doc_info);
              threedproblem.time() = 0.0;
              dt = 0.03;

              // Allow double the initial discretisation error in final error.
              double errtol = threedproblem.compare_m_with_function(initial_conds[i_ic])
                * 4;

              std::cout << std::endl
                        << "error = " << threedproblem.compare_m_with_function(initial_conds[i_ic])
                        << std::endl;

              // Timestep to end
              while(threedproblem.time_pt()->time() < tmax)
                {
                  std::cout << "step number = " << doc_info.number()
                            << ", time = " << threedproblem.time_pt()->time()
                            << std::endl;

                  // Step (use recommended next time step from previous step).
                  dt = threedproblem.adaptive_unsteady_newton_solve(dt,eps);

                  std::cout << " dt = " << dt << std::endl;

                  // Output
                  threedproblem.doc_solution(doc_info);

                  std::cout << "error = " << threedproblem.compare_m_with_function(initial_conds[i_ic])
                            << std::endl;
                }

              // Check that m hasn't changed, error if it has changed more than tol
              if((threedproblem.compare_m_with_function(initial_conds[i_ic]))
                 > errtol)
                {
                  return 300 + i_ic*threedas.size() + i+1;
                }
              else
                {
                  std::cout << "Ok for a = " << Inputs::a << std::endl;
                }

            }
        }

    }

  return 0;
}
