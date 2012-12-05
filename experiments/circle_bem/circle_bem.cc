

#include "generic.h"
#include "../../semi_implicit_problem.h"

#include "meshes/triangle_mesh.h"

// Floating point error checks
#include <fenv.h>

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
    m[2] = 0.0;

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


  void domain_wall_applied_field(const double& t, const Vector<double> &x,
                                 Vector<double> &h_app)
  {
    h_app.assign(3,0.0);
    h_app[0] = std::tanh(x[2] - 50) * 60000;
  }

}


int main(int argc, char** argv)
{
  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  std::string mesh_name("circle.fig.1"), corners_file_name("");

  CommandLineArgs::setup(argc,argv);
  CommandLineArgs::specify_command_line_flag("-mesh", &mesh_name); 
  CommandLineArgs::specify_command_line_flag("-corners", &corners_file_name); 
  CommandLineArgs::parse_and_assign();
  CommandLineArgs::output();

  const unsigned dim = 2, nnode1d = 2;

  BDF<2> ts(true); 

  TriangleMesh<TMagnetostaticFieldElement<dim,nnode1d> >
    phi_1_mesh(mesh_name+".node", mesh_name+".ele", mesh_name+".poly");
  TriangleMesh<TMagnetostaticFieldElement<dim,nnode1d> >
    phi_mesh(mesh_name+".node", mesh_name+".ele", mesh_name+".poly");
  TriangleMesh<TSemiImplicitMicromagElement<dim,nnode1d> >
    llg_mesh(mesh_name+".node", mesh_name+".ele", mesh_name+".poly", &ts);


//   // Flip meshes
// #warning flipping meshes
//   FiniteElement::Accept_negative_jacobian=true;
//   for(unsigned nd=0, nnode=phi_1_mesh.nnode(); nd<nnode; nd++)
//     {          
//       Node* nd_pt = phi_1_mesh.node_pt(nd); 
//       nd_pt->x(0) *= -1;
//     }

//   for(unsigned nd=0, nnode=phi_mesh.nnode(); nd<nnode; nd++)
//     {          
//       Node* nd_pt = phi_mesh.node_pt(nd); 
//       nd_pt->x(0) *= -1;
//     }
//   for(unsigned nd=0, nnode=llg_mesh.nnode(); nd<nnode; nd++)
//     {          
//       Node* nd_pt = llg_mesh.node_pt(nd); 
//       nd_pt->x(0) *= -1;
//     }

  // For some reason we have to do this manually...
  llg_mesh.setup_boundary_element_info();
  phi_mesh.setup_boundary_element_info();
  phi_1_mesh.setup_boundary_element_info();

  // Corner (fractional) angles for regular polygon with n sides
  double n = 24;
  double angle_deg = ((n - 2) * 180 / n);
  double fractional_angle = angle_deg/360;

  // Read in corner location data (if file specified)
  std::ifstream corners_file(corners_file_name.c_str());
  Vector< std::pair<Vector<double>, double> > corner_node_position_angle;
  Vector< std::pair<Vector<double>, double> >* 
    corner_node_position_angle_pt = &corner_node_position_angle;

  if(corners_file.is_open())
    {
      while ( corners_file.good() )
        {
          Vector<double> x(2,0.0);
          corners_file >> x[0] >> x[1];

          std::cout << x << std::endl;

          // // Flip mesh
          // #warning flipping angles too
          // x[0] *= -1;
      
          std::pair<Vector<double>, double> node(x, fractional_angle);
          corner_node_position_angle.push_back(node);
        }
    }
  else
    {
      std::cout << "Failed to find corners file." << std::endl;
      corner_node_position_angle_pt = 0;
    }

  SemiImplicitHybridMicromagneticsProblem<TMagnetostaticFieldElement<dim,nnode1d>,
                                          TSemiImplicitMicromagElement<dim,nnode1d>
                                          >
  problem(&phi_1_mesh,
          &phi_mesh,
          &llg_mesh,
          //&Inputs::domain_wall_applied_field);
          &Inputs::no_applied_field,
          corner_node_position_angle_pt);


  // Set up time stepping
  problem.set_initial_condition(Inputs::initial_m);
  double dt = 0.03; // initial suggestion for dt
  const double tmax = 100;
  const double eps = 1e-3;
  bool adaptive = true;

  problem.llg_sub_problem_pt()->max_newton_iterations() = 20;
  problem.llg_sub_problem_pt()->max_residuals() = 40;

  // Since we are using bdf2 we probably need to keep renormalising...
  problem.llg_sub_problem_pt()->renormalise_each_time_step() = true;

  // dump BEM
  problem.bem_matrix_pt()->output("bem_matrix");

  // Set up output
  DocInfo doc_info;
  doc_info.set_directory("results");

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

          // print average field
          Vector<double> avg_field;
          problem.average_magnetostatic_field(avg_field);

          std::cout << avg_field << std::endl;
#warning exit after one solve
          return 0;
        }
    }

}


