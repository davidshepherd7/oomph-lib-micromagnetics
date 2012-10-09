
/*
description of file goes here
*/

#include "generic.h"

#include "../generic_poisson_problem.h"

#include "meshes/simple_rectangular_quadmesh.h"

using namespace oomph;

using namespace std;

//===== start_of_namespace=============================================
/// Namespace for exact solution for Poisson equation with "sharp step" 
//=====================================================================
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


int main()
{

 // Make a mesh
 Mesh* Bulk_mesh_pt = 
  new SimpleRectangularQuadMesh<QPoissonElement<2,3> >(4,4,1.0,2.0);

 // Make a problem
 GenericPoissonProblem<QPoissonElement<2,3> > problem;

 // Assign bulk mesh (flux mesh is automatically dealt with)
 problem.set_bulk_mesh(Bulk_mesh_pt);

 // Assign b.c.s
 for(unsigned b=0; b < problem.bulk_mesh_pt()->nboundary(); b++)
  {
   if(b != 1) problem.set_dirichlet_boundary(b, &TanhSolnForPoisson::get_exact_u);
  }
 problem.set_neumann_boundary
  (1, &TanhSolnForPoisson::prescribed_flux_on_fixed_x_boundary);

 // Assign function pointers
 problem.set_source_fct_pt(&TanhSolnForPoisson::source_function);
 problem.exact_solution_fct_pt() = &TanhSolnForPoisson::get_exact_u;
 
 // Finish building the problem
 problem.build();


 // Go!


 // Create label for output
 //------------------------
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("results");

 // Step number
 doc_info.number()=0;


 // Check if we're ready to go:
 //----------------------------
 cout << "\n\n\nProblem self-test ";
 if (problem.self_test()==0) 
  {
   cout << "passed: Problem can be solved." << std::endl;
  }
 else 
  {
   throw OomphLibError("Self test failed",
                       "main()",
                       OOMPH_EXCEPTION_LOCATION);
  }

 
 // Set the orientation of the "step" to 45 degrees
 TanhSolnForPoisson::TanPhi=1.0;
 
 // Initial value for the steepness of the "step"
 TanhSolnForPoisson::Alpha=1.0; 

 // Do a couple of solutions for different forcing functions
 //---------------------------------------------------------
 unsigned nstep=4;
 for (unsigned istep=0;istep<nstep;istep++)
  {
   // Increase the steepness of the step:
   TanhSolnForPoisson::Alpha+=2.0;

   cout << "\n\nSolving for TanhSolnForPoisson::Alpha="
        << TanhSolnForPoisson::Alpha << std::endl << std::endl;

   // Solve the problem
   problem.newton_solve();

   //Output solution
   problem.doc_solution(doc_info);
 
   //Increment counter for solutions 
   doc_info.number()++; 
  }


}
