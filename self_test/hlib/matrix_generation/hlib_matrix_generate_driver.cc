
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
    HMatrix hmat(50);
  }

  // Make the problem to get bem mesh etc.
  // ============================================================

  // Need a dummy problem so that we can get equation numbers :(
  MMArgs args;
  args.parse(argc, argv);
  LLGProblem problem;
  problem.Residual_calculator_pt
    = LLGFactories::residual_calculator_factory("llg");
  problem.add_time_stepper_pt(args.time_stepper_pt);
  args.assign_specific_parameters(&problem);
  // disable corner angles in bem matrx, to make matrices easier to compare
  problem.Disable_bem_corners = true;
  problem.build(args.mesh_pts);

  // Some useful pointers
  DenseDoubleMatrix* old_bem_matrix_pt = problem.Bem_handler_pt->bem_matrix_pt();
  const Mesh* surface_mesh_pt = problem.Bem_handler_pt->bem_mesh_pt();

  // Vectors for multiply test
  const unsigned nbemnodes = surface_mesh_pt->nnode();
  LinearAlgebraDistribution dist(0, nbemnodes, false);
  DoubleVector x(dist), H_soln(dist), dense_H_soln(dist), oomph_soln(dist);

  // Random entries for x vector
  Vector<double> v = VectorOps::random_vector(nbemnodes);
  x.initialise(v);


  // Multiply with pure-oomph bem matrix
  old_bem_matrix_pt->multiply(x, oomph_soln);


  // Make a H matrix and use it in a multiply
  // ============================================================

  // Use the bem mesh to make a h-matrix

  // unsigned nmin = 500; // ??ds works
  unsigned nmin = 30;
  HMatrix hmat(nmin);
  hmat.build(*surface_mesh_pt, problem.Bem_handler_pt);

  // Dump rank info (needs to be in a variable, not a string constant
  // beacuse the C code could technical write to the string [not const]).
  char rank_file_name[] = "rank.ps";
  outputrank_supermatrix(hmat.supermatrix_pt(), rank_file_name);

  // Do the multiply
  hmat.multiply(x, H_soln);

  // Flip the sign ??ds not sure why we need to do this
  H_soln *= -1;


  // Make a dense matrix from the H matrix and try multiplying with that
  // ============================================================

  // Copy H matrix into an oomph-lib matrix
  DenseDoubleMatrix double_matrix;
  hmat.todense(double_matrix);

  // Flip sign ??ds not sure why we need this to make it match
  for(unsigned i=0; i<double_matrix.nrow(); i++)
    {
      for(unsigned j=0; j<double_matrix.ncol(); j++)
        {
          double_matrix(i,j) *= -1;
        }
    }

  // Multiply
  double_matrix.multiply(x, dense_H_soln);



  // // Output the timing and error results
  // // ============================================================
  // std::cout << "H-matrix build time " << hbuildstop - hbuildstart << std::endl;
  // std::cout << "dense build time " << dbuildstop - dbuildstart << std::endl;
  // std::cout << "oldbem build time " << oldbembuildstop - oldbembuildstart << std::endl;
  // std::cout << std::endl;
  // std::cout << "H-matrix multiply time " <<  hstop - hstart << std::endl;
  // std::cout << "dense multiply time " << dstop - dstart << std::endl;
  // std::cout << "old bem dense multiply time " << oldbemstop - oldbemstart
  //           << std::endl;


  // Dump matrices to files
  // ============================================================
  double_matrix.output("Validation/new_bem_matrix");
  old_bem_matrix_pt->output("Validation/old_bem_matrix");


  double H_soln_vs_dense = rel_two_norm_diff(H_soln, dense_H_soln);
  double H_soln_vs_oomph = rel_two_norm_diff(H_soln, oomph_soln);
  double dense_vs_oomph = rel_two_norm_diff(dense_H_soln, oomph_soln);


  // Some analysis
  // ============================================================
  std::cout << std::endl;
  std::cout << "muliplication 2 norm of H-matrix vs dense version: "
            << H_soln_vs_dense << std::endl;
  std::cout << "muliplication 2 norm of H-matrix vs oomph bem: "
            << H_soln_vs_oomph << std::endl;
  std::cout << "muliplication 2 norm of dense version of H-matrix vs oomph bem: "
            << dense_vs_oomph << std::endl;

  // ??ds need to fix ordering of H-matrix multiply output (and input?)
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
