

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

  CommandLineArgs::setup(argc,argv);

  const unsigned dim = 3, nnode1d = 2;

  BDF<2> ts(true);

  // nmag style meshes
  double lx = 30, ly = lx, lz = 100;
  unsigned nx = 7, ny = nx, nz = (int(lz/lx) + 1) * nx;

  // SimpleCubicTetMesh<TMagnetostaticFieldElement<dim,nnode1d> >
  //   phi_1_mesh(nx,ny,nz, lx,ly,lz);
  // SimpleCubicTetMesh<TMagnetostaticFieldElement<dim,nnode1d> >
  //   phi_mesh(nx,ny,nz, lx,ly,lz);
  // SimpleCubicTetMesh<TSemiImplicitMicromagElement<dim,nnode1d> >
  //   llg_mesh(nx,ny,nz, lx,ly,lz, &ts);

  // // // Unstructured meshes
  // // TetgenMesh<TMagnetostaticFieldElement<dim,nnode1d> >
  // //   phi_1_mesh("cubeoid.1.node", "cubeoid.1.ele", "cubeoid.1.face");
  // // TetgenMesh<TMagnetostaticFieldElement<dim,nnode1d> >
  // //   phi_mesh("cubeoid.1.node", "cubeoid.1.ele", "cubeoid.1.face");
  // // TetgenMesh<TSemiImplicitMicromagElement<dim,nnode1d> >
  // //   llg_mesh("cubeoid.1.node", "cubeoid.1.ele", "cubeoid.1.face", &ts);

  // // For some reason we have to do this manually...
  // llg_mesh.setup_boundary_element_info();
  // phi_mesh.setup_boundary_element_info();
  // phi_1_mesh.setup_boundary_element_info();

  // SemiImplicitHybridMicromagneticsProblem<TMagnetostaticFieldElement<dim,nnode1d>,
  //                                         TSemiImplicitMicromagElement<dim,nnode1d>
  //                                         >
  // problem(&phi_1_mesh,
  //         &phi_mesh,
  //         &llg_mesh,
  //         //&Inputs::domain_wall_applied_field);
  //          &Inputs::no_applied_field);



  SimpleCubicMesh<QMagnetostaticFieldElement<dim,nnode1d> >
    phi_1_mesh(nx,ny,nz, lx,ly,lz);
  SimpleCubicMesh<QMagnetostaticFieldElement<dim,nnode1d> >
    phi_mesh(nx,ny,nz, lx,ly,lz);
  SimpleCubicMesh<QSemiImplicitMicromagElement<dim,nnode1d> >
    llg_mesh(nx,ny,nz, lx,ly,lz, &ts);

  // For some reason we have to do this manually...
  llg_mesh.setup_boundary_element_info();
  phi_mesh.setup_boundary_element_info();
  phi_1_mesh.setup_boundary_element_info();

  SemiImplicitHybridMicromagneticsProblem<QMagnetostaticFieldElement<dim,nnode1d>,
                                          QSemiImplicitMicromagElement<dim,nnode1d>
                                          >
  problem(&phi_1_mesh,
          &phi_mesh,
          &llg_mesh,
          //&Inputs::domain_wall_applied_field);
          &Inputs::no_applied_field);


  // Set up the magnetic parameters
  problem.mag_parameters_pt()->set_nmag_rectangle();

  // Debugging stuff:
  // ============================================================ 

// #warning using -1*alpha
//   problem.mag_parameters_pt()->gilbert_damping() *= -1;

  #warning turned off boundary exchange effects
  problem.mag_parameters_pt()->boundary_exchange_debug_coeff() = 0.0;

//   # warning flipped sign of hms
//   problem.mag_parameters_pt()->magnetostatic_debug_coeff() *= -1;

  // # warning turned off hms and added a field
  // problem.mag_parameters_pt()->magnetostatic_debug_coeff() = 0.0;

  // # warning reduced hms 
  // problem.mag_parameters_pt()->magnetostatic_debug_coeff() = 0.0;

  // #warning strengthened exchange
  // problem.mag_parameters_pt()->exchange_debug_coeff() = 100.0;

  // #warning disabled exchange
  // problem.mag_parameters_pt()->exchange_debug_coeff() = 0.0;

  // ============================================================ 



  // problem.llg_sub_problem_pt()->linear_solver_pt() = new GMRES<CRDoubleMatrix>;

  // Set up time stepping
  problem.set_initial_condition(Inputs::initial_m);
  double dt = 0.03; // initial suggestion for dt
  const double tmax = 100;
  const double eps = 1e-5;
  bool adaptive = true;

  // //??ds limit dt to check if problems are due to non-conservation of m on
  // // long timesteps....
  // problem.llg_sub_problem_pt()->maximum_dt() = 0.1;

  problem.llg_sub_problem_pt()->max_newton_iterations() = 20;
  problem.llg_sub_problem_pt()->max_residuals() = 40;

  // Since we are using bdf2 we probably need to keep renormalising...
  problem.llg_sub_problem_pt()->renormalise_each_time_step() = true;

  //??ds
  // DoubleVector dummy;
  // CRDoubleMatrix J;
  // problem.llg_sub_problem_pt()->get_jacobian(dummy,J);
  // J.sparse_indexed_output("jacobian");

  // // Check mesh is where I think it is
  // std::ofstream mesh_file("flux_mesh");
  // Mesh* flux_pt = problem.llg_sub_problem_pt()->surface_exchange_mesh_pt();
  // std::cout <<  flux_pt->nelement() << std::endl;
  // for(unsigned e=0, ne=flux_pt->nelement(); e < ne; e++)
  //   {
  //     FiniteElement* ele_pt = flux_pt->finite_element_pt(e);      
  //     for(unsigned nd=0, nnode=ele_pt->nnode(); nd<nnode; nd++)
  //       {          
  //         Node* nd_pt = ele_pt->node_pt(nd);
  //         mesh_file << nd_pt->x(0) << " " << nd_pt->x(1) << " " << nd_pt->x(2) << std::endl;
  //       }
  //   }
  // mesh_file.close();

  // dump BEM
  problem.bem_matrix_pt()->output("bem_matrix");

  // Set up output
  DocInfo doc_info;
  doc_info.set_directory("results");
  // problem.doc_solution(doc_info);

  // ConvergenceData conv_dat;
  // problem.llg_sub_problem_pt()->convergence_data_pt() = &conv_dat;
  // problem.llg_sub_problem_pt()->convergence_data_pt()->write_headers("convergence_data");

  // Timestep to end
  if(adaptive)
    {
      while(problem.llg_sub_problem_pt()->time() < tmax)
        {
          std::cout << "step number = " << doc_info.number()
                    << ", time = " << problem.llg_sub_problem_pt()->time()
                    << ", dt = " << dt << std::endl;
          // Step
          dt = problem.semi_implicit_step(dt,eps);

          problem.doc_solution(doc_info);
        }
    }
  // else
  //   {
  //     while(problem.llg_sub_problem_pt()->time() < tmax)
  //       {
  //         std::cout << "step number = " << doc_info.number()
  //                   << ", time = " << problem.llg_sub_problem_pt()->time()
  //                   << " dt = " << dt << std::endl;
  //         // Step
  //         problem.semi_implicit_step(dt);

  //         problem.doc_solution(doc_info);
  //       }
  //   }

  // conv_dat.output(std::cout);

}
