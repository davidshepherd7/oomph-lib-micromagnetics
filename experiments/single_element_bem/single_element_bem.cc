/*
  description of file goes here
*/

#include "generic.h"
#include "../../semi_implicit_problem.h"

// Floating point error checks
#include <fenv.h>

// Mesh
#include "./single_element_mesh.h"

using namespace oomph;
using namespace MathematicalConstants;
using namespace StringConversion;


namespace Inputs
{

  void initial_m(const double& t, const Vector<double> &x, Vector<double> &m)
  {
    m.assign(3,0.0);

    m[0] = 0.0;
    m[1] = 0.0;
    m[2] = 1.0;

    VectorOps::normalise(m);
  }

  // Turn off field
  void no_applied_field(const double& t, const Vector<double> &x,
                        Vector<double> &h_app)
  {
    h_app.assign(3,0.0);
  }

}


int main(int argc, char *argv[])
{
  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  // Store command line args
  CommandLineArgs::setup(argc,argv);

  // Make single element meshes
  SingleTetElementMesh<TMagnetostaticFieldElement<3,2> >
    phi_mesh(1.0, 1.0, 1.0), phi1_mesh(1.0, 1.0, 1.0);

  SingleTetElementMesh<TSemiImplicitMicromagElement<3,2> > 
    llg_mesh(1.0, 1.0, 1.0) ;

  //??ds - must be a more elegant way to have two identical meshes... or
  // avoid having two altogether... how can we be sure that we don't break
  // everything by applying changes (e.g. refinement) to only one mesh?

  phi1_mesh.output("mesh.dat",2);


  SemiImplicitHybridMicromagneticsProblem<TMagnetostaticFieldElement<3,2>,
                                          TSemiImplicitMicromagElement<3,2>
                                          >
    problem(&phi1_mesh, &phi_mesh, &llg_mesh, &Inputs::no_applied_field,0,false);


  for(unsigned nd=0, nnode=problem.bem_mesh_pt()->nnode(); nd<nnode; nd++)
    {          
      Node* nd_pt = problem.bem_mesh_pt()->node_pt(nd); 
      Vector<double> x(3,0.0);
      nd_pt->position(x);
      std::cout << nd << " " << x << std::endl;
    }

  for(unsigned b=0, nb=phi1_mesh.nboundary(); b<nb; b++)
    {
      std::cout <<  "" << std::endl;
      std::cout << b << std::endl;
      for(unsigned nd=0, nnd=phi1_mesh.nboundary_node(b); nd<nnd; nd++)
        {
          Node* nd_pt = phi1_mesh.boundary_node_pt(b,nd);
          Vector<double> x(3,0.0);
          nd_pt->position(x);
          std::cout <<  nd << " " << x << std::endl;
        }
    }

// Set up the magnetic parameters
problem.mag_parameters_pt()->set_nmag_rectangle();

// Set up time stepping
problem.set_initial_condition(Inputs::initial_m);

// Set up output
DocInfo doc_info;
doc_info.set_directory("results");


// Solve it for a tiny step (don't care about m)
problem.semi_implicit_step(1e-9);

// // Output BEM
// problem.bem_handler_pt()->boundary_matrix_pt()->output("bem_matrix");

// Output results
problem.doc_solution(doc_info);

return 0;
}
