
#include "generic.h"
#include "micromag.h"
#include "../../src/meshes/tetgen_mesh.h"

using namespace TimingHelpers;
using namespace Factories;



int main(int argc, char *argv[])
{

  // First make sure that we can create/destroy empty h matrices correctly
  // ============================================================
  {
    HMatrix hmat;
  }


  // Make the problem to get bem mesh etc.
  // ============================================================

  // Need a dummy problem so that we can get equation numbers... :(
  MMArgs args;
  args.parse(argc, argv);
  LLGProblem problem;
  problem.Residual_calculator_pt
    = LLGFactories::residual_calculator_factory("llg");
  problem.add_time_stepper_pt(args.time_stepper_pt);
  args.assign_specific_parameters(&problem);
  problem.build(args.mesh_pts);

  // No corner data
  CornerDataInput input_corner_data;

  // Add all boundaries of all meshes to bem boundary list
  BemBoundaryData bem_boundaries;
  for(unsigned msh=0, nmsh=args.mesh_pts.size(); msh<nmsh; msh++)
    {
      Mesh* mesh_pt = args.mesh_pts[msh];
      for(unsigned b=0, nb=mesh_pt->nboundary(); b<nb; b++)
        {
          bem_boundaries.push_back(std::make_pair(b, mesh_pt));
        }
    }

  // Get the phi/phi1 indicies
  MicromagEquations* e_pt = checked_dynamic_cast<MicromagEquations*>
    (args.mesh_pts[0]->element_pt(0));
  unsigned bem_phi_index = e_pt->phi_index_micromag();
  unsigned bem_phi_1_index = e_pt->phi_1_index_micromag();

  // Build the handler. Disable corner angles in bem matrix, to make
  // matrices easier to compare.
  BoundaryElementHandler Bem_handler;
  Factories::bem_handler_factory(Bem_handler, bem_boundaries,
                                 bem_phi_index,
                                 bem_phi_1_index,
                                 input_corner_data,
                                 false, true, -1);


  // Some useful pointers
  DoubleMatrixBase* old_bem_matrix_pt = Bem_handler.bem_matrix_pt();
  const Mesh* surface_mesh_pt = Bem_handler.bem_mesh_pt();


  // Test the differences in the matrices by multiplying them each with a
  // random vector
  // ============================================================


  // Vectors for multiplication test
  const unsigned nbemnodes = surface_mesh_pt->nnode();
  LinearAlgebraDistribution dist(0, nbemnodes, false);
  DoubleVector x(dist), H_soln(dist), dense_H_soln(dist), oomph_soln(dist);

  // Fill x vector with random entries
  Vector<double> v = VectorOps::random_vector(nbemnodes);
  x.initialise(v);

  // Multiply with pure-oomph bem matrix
  old_bem_matrix_pt->multiply(x, oomph_soln);


  // Use the bem mesh to make a h-matrix
  HMatrix hmat;
  hmat.build(*surface_mesh_pt, &Bem_handler);

  // Dump rank info
  outputrank_supermatrix(hmat.supermatrix_pt(), "rank.ps");

  // Do the multiply
  hmat.multiply(x, H_soln);


  // Copy H matrix into an oomph-lib matrix
  DenseDoubleMatrix double_matrix;
  hmat.to_dense(double_matrix);

  // Multiply
  double_matrix.multiply(x, dense_H_soln);



  // Dump matrices to files
  double_matrix.output("Validation/new_bem_matrix");
  checked_dynamic_cast<DenseDoubleMatrix*>(old_bem_matrix_pt)
    ->output("Validation/old_bem_matrix");



  // Some analysis of how good the approximation is
  // ============================================================

  double H_soln_vs_dense = rel_two_norm_diff(H_soln, dense_H_soln);
  double H_soln_vs_oomph = rel_two_norm_diff(H_soln, oomph_soln);
  double dense_vs_oomph = rel_two_norm_diff(dense_H_soln, oomph_soln);

  std::cout << std::endl;
  std::cout << "muliplication 2 norm of H-matrix vs dense version: "
            << H_soln_vs_dense << std::endl;
  std::cout << "muliplication 2 norm of H-matrix vs oomph bem: "
            << H_soln_vs_oomph << std::endl;
  std::cout << "muliplication 2 norm of dense version of H-matrix vs oomph bem: "
            << dense_vs_oomph << std::endl;

  if(H_soln_vs_dense > 1e-5)
    {
      std::cout << "Test failed due to error in dense H." << std::endl;
      return 1;
    }

  if(H_soln_vs_oomph > 1e-5)
    {
      std::cout << "Test failed due to error in oomph vs H." << std::endl;
      return 2;
    }


  return 0;
}
