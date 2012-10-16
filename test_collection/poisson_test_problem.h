#ifndef OOMPH_POISSON_TEST_PROBLEM_H
#define OOMPH_POISSON_TEST_PROBLEM_H

/*
  description of file goes here
*/

#include "generic.h"

#include "../generic_poisson_problem.h"

// Mesh for running Poisson tests
#include "meshes/simple_rectangular_quadmesh.h"

// =================================================================
/// Functions for Poisson tests
// =================================================================
namespace TanhSolnForPoisson
{

 /// Parameter for steepness of "step"
 double Alpha=1.0;

 /// Parameter for angle Phi of "step"
 double TanPhi=0.0;

 /// Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0] = tanh(1.0-Alpha*(TanPhi*x[0]-x[1]));
 }

 /// Source function required to make the solution above an exact solution
 void source_function(const Vector<double>& x, double& source)
 {
  source = 2.0*tanh(-1.0+Alpha*(TanPhi*x[0]-x[1]))*
   (1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*
   Alpha*Alpha*TanPhi*TanPhi+2.0*tanh(-1.0+Alpha*(TanPhi*x[0]-x[1]))*
   (1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*Alpha*Alpha;
 }

 /// Flux required by the exact solution on a boundary on which x is fixed
 void prescribed_flux_on_fixed_x_boundary(const Vector<double>& x,
                                          double& flux)
 {
  //The outer unit normal to the boundary is (1,0)
  double N[2] = {1.0, 0.0};
  //The flux in terms of the normal is
  flux =
   -(1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*Alpha*TanPhi*N[0]+(
                                                                          1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*Alpha*N[1];
 }

} // end of namespace

using namespace oomph;

namespace oomph
{
 // =================================================================
 /// Basic Poisson problem for tests
 // =================================================================
 class GenericPoissonForTests :
 public GenericPoissonProblem<QPoissonElement<2,3> >
  {
  public:
   /// Constructor - build a Poisson problem with some of each type of
   /// boundary and preset source function.
   GenericPoissonForTests()
    {
     // Make a mesh
     Mesh* Bulk_mesh_pt =
      new SimpleRectangularQuadMesh<QPoissonElement<2,3> >(4,4,1.0,2.0);

     // Assign bulk mesh (flux mesh is automatically dealt with)
     this->set_bulk_mesh(Bulk_mesh_pt);

     // Assign b.c.s
     for(unsigned b=0; b < this->bulk_mesh_pt()->nboundary(); b++)
      {
       if(b != 1) this->set_dirichlet_boundary(b, &TanhSolnForPoisson::get_exact_u);
      }
     this->set_neumann_boundary
      (1, &TanhSolnForPoisson::prescribed_flux_on_fixed_x_boundary);

     // Assign function pointers
     this->set_source_fct_pt(&TanhSolnForPoisson::source_function);
     this->exact_solution_fct_pt() = &TanhSolnForPoisson::get_exact_u;

     // Finish building the problem
     this->build();
    }


   int run_test()
   {
    // Self test
    if (self_test() != 0)
     {
      std::cerr << "Self test failed" << std::endl;
      return 5;
     }

    // Set the orientation of the "step" to 45 degrees
    TanhSolnForPoisson::TanPhi=1.0;

    // Initial value for the steepness of the "step"
    TanhSolnForPoisson::Alpha=1.0;

    // List of roughly the error norms allowed for different alpha values
    // (based on previous runs).
    Vector<double> Allowed_error;
    Allowed_error.push_back(0.01);
    Allowed_error.push_back(0.05);
    Allowed_error.push_back(0.1);
    Allowed_error.push_back(0.5);

    // Do a couple of solutions for different forcing functions
    unsigned nstep=4;
    for (unsigned istep=0;istep<nstep;istep++)
     {
      // Increase the steepness of the step:
      TanhSolnForPoisson::Alpha+=2.0;

      // Solve the problem
      newton_solve();

      if(get_error_norm() > Allowed_error[istep])
       {
        std::cerr << "Generic poisson test FAILED: error norm is too large"
                  << std::endl;
        return 1;
       }
     }

    return 0;
   }
  };

}

#endif
