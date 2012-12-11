

#include "generic.h"
#include "../../semi_implicit_problem.h"

// Mesh
#include "meshes/simple_rectangular_tri_mesh.h"
#include "meshes/tetgen_mesh.h"
#include "meshes/simple_cubic_tet_mesh.h"
#include "meshes/simple_cubic_mesh.h"

// Floating point error checks
#include <fenv.h>



/*

  Idea: calculate magnetostatic field then run implicit llg with the
  magnetostatic field included in the applied field.

*/


using namespace oomph;
using namespace MathematicalConstants;
using namespace StringConversion;

namespace Inputs
{

  void initial_m(const double& t, const Vector<double> &x, Vector<double> &m)
  {
    m.assign(3,0.0);

    m[0] = 1.0;
    m[1] = 0.0;
    m[2] = 1.0;
    
    // #warning wrong initial m
    // m[0] = sin(x[0]*2*Pi);
    // m[1] = cos(x[1]*2*Pi);
    // m[2] = 1.0 - m[0] - m[1];

    // m[0] = x[0]/2 + x[1]/2;
    // m[1] = (1 - m[0]);
    // m[2] = 0.0;

    VectorOps::normalise(m);
  }

  // Turn off field
  void no_applied_field(const double& t, const Vector<double> &x,
                        Vector<double> &h_app)
  {
    h_app.assign(3,0.0);
  }


  // void domain_wall_applied_field(const double& t, const Vector<double> &x,
  //                       Vector<double> &h_app)
  // {
  //   h_app.assign(3,0.0);
  //   h_app[0] = std::tanh(x[2] - 50) * 60000;
  // }

}


int main(int argc, char** argv)
{
  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);


  // Default values
  double lx = 30, ly = lx, lz = 100;
  unsigned nx = 7;
  std::string outdir("results");

  // Get command line args
  CommandLineArgs::setup(argc,argv);
  CommandLineArgs::specify_command_line_flag("-nx", &nx);
  CommandLineArgs::specify_command_line_flag("-outdir", &outdir);
  CommandLineArgs::parse_and_assign();
  CommandLineArgs::output();

  const unsigned dim = 3, nnode1d = 2;
  unsigned ny = nx, nz = (int(lz/lx) + 1) * nx; // calculate relative to nx

  BDF<2> ts(true);

  Mesh *phi1_mesh_pt(0), *phi_mesh_pt(0), *llg_mesh_pt(0);

  // Structured meshes
  phi1_mesh_pt = new SimpleCubicMesh<QMagnetostaticFieldElement<dim,nnode1d> >
    (nx,ny,nz, lx,ly,lz);
  phi_mesh_pt = new SimpleCubicMesh<QMagnetostaticFieldElement<dim,nnode1d> >
    (nx,ny,nz, lx,ly,lz);
  llg_mesh_pt = new SimpleCubicMesh<QSemiImplicitMicromagElement<dim,nnode1d> >
    (nx,ny,nz, lx,ly,lz, &ts);
      
  // For some reason we have to do this manually...
  llg_mesh_pt->setup_boundary_element_info();
  phi_mesh_pt->setup_boundary_element_info();
  phi1_mesh_pt->setup_boundary_element_info();
    

  // Build problem
  SemiImplicitHybridMicromagneticsProblem
    <QMagnetostaticFieldElement<dim,nnode1d>,
     QSemiImplicitMicromagElement<dim,nnode1d> > problem
  (phi1_mesh_pt, phi_mesh_pt, llg_mesh_pt, &Inputs::no_applied_field);


  // Set up the magnetic parameters
  problem.mag_parameters_pt()->set_nmag_rectangle();

  // problem.llg_sub_problem_pt()->linear_solver_pt() = new GMRES<CRDoubleMatrix>;

  // Set up time stepping parameters
  problem.set_initial_condition(Inputs::initial_m);
  double dt = 0.03; // initial suggestion for dt
  const double tmax = 100;
  const double eps = 1e-5;

  problem.llg_sub_problem_pt()->max_newton_iterations() = 20;
  problem.llg_sub_problem_pt()->max_residuals() = 40;

  // Since we are using bdf2 we probably need to keep renormalising...
  problem.llg_sub_problem_pt()->renormalise_each_time_step() = true;

  // Set up output
  DocInfo doc_info;
  doc_info.set_directory(outdir);
  problem.doc_solution(doc_info);

  // ConvergenceData conv_dat;
  // problem.llg_sub_problem_pt()->convergence_data_pt() = &conv_dat;
  // problem.llg_sub_problem_pt()->convergence_data_pt()->write_headers("convergence_data");

  // Timestep to end
  while(problem.llg_sub_problem_pt()->time() < tmax)
    {
      std::cout << "step number = " << doc_info.number()
                << ", time = " << problem.llg_sub_problem_pt()->time()
                << ", dt = " << dt << std::endl;
      // Step
      dt = problem.semi_implicit_step(dt,eps);

      // Output
      problem.doc_solution(doc_info);
    }

  // conv_dat.output(std::cout);

}
