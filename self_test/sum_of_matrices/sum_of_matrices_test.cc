

#include "generic.h"
#include "../../vector_helpers.h"

// Accumulate command
#include <numeric>

// Floating point debugging
#include <fenv.h>

using namespace oomph;

namespace oomph
{

  void identity_map(std::map<long unsigned, long unsigned> &a, const unsigned long &n)
  {
    for(unsigned long k=0; k< n; k++)
      {
        std::pair<long unsigned, long unsigned> my_pair(k,k);
        a.insert(my_pair);
      }
  }


} // End of oomph namespace


int main()
{
  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  const unsigned mat_size = 50;
  const unsigned n_sparse_values = mat_size;
  const unsigned dense_block_size = mat_size/10;
  const unsigned max_val = 100;
  const unsigned n_tests = 10;

  Vector<double> norm_error(0);

  for(unsigned i_rand=0; i_rand<n_tests; i_rand++)
    {
      srand(i_rand);

      std::cout << "Matrix size is " << mat_size << std::endl;

      // A linear algebra "distribution": specifies no inter-processor communicator
      // and mat_size rows.
      OomphCommunicator communicator;
      LinearAlgebraDistribution dist(&communicator,mat_size,false);

      // Create empty sum of matrices matrix
      SumOfMatrices sum_matrix;


      // A random compressed row form sparse matrix, to avoid overlaps we
      // have exactly one element per row.
      // ============================================================
      CRDoubleMatrix main_matrix(&dist);
      VectorOps::random_single_element_per_row_cr_matrix(main_matrix, mat_size,
                                                         n_sparse_values, max_val);
      sum_matrix.main_matrix_pt() = &main_matrix;

      // Create a diagonal sparse matrix and add to sum.
      // ============================================================
      CRDoubleMatrix id_matrix(&dist);
      VectorOps::diag_cr_matrix(id_matrix, mat_size, 2*max_val);
      NodeGlobalNumbersLookup row_map, col_map;
      row_map.build_identity_map(mat_size);
      col_map.build_identity_map(mat_size);
      sum_matrix.add_matrix(&id_matrix,&row_map,&col_map);

      // Create a random dense matrix and add to sum in top right hand
      // corner.
      // ============================================================
      DenseDoubleMatrix dense_block(dense_block_size,dense_block_size,0.0);
      for(unsigned i=0; i<dense_block_size; i++)
        {
          for(unsigned j=0; j<dense_block_size; j++)
            {
              dense_block(i,j) = rand() % max_val;
            }
        }
      NodeGlobalNumbersLookup dense_row_map, dense_col_map;
      dense_row_map.build_identity_map(dense_block_size);
      dense_col_map.build_identity_map(dense_block_size);

      sum_matrix.add_matrix(&dense_block, &dense_row_map, &dense_col_map);

      // sum_matrix.sparse_indexed_output(std::cout);

      // A random rhs vector
      // ============================================================
      DoubleVector rhs;
      rhs.build(dist,1.0);
      for(unsigned i=0; i<mat_size; i++) rhs[i] = rand() % max_val;



      // Make a total matrix to compare
      // ============================================================
      Vector<int> sum_rows, sum_cols, sum_row_start;
      Vector<double> sum_values;
      VectorOps::get_as_indicies(sum_matrix, sum_values, sum_cols, sum_rows);

      VectorOps::rowindex2rowstart(sum_rows, sum_row_start);
      CRDoubleMatrix total_matrix;
      total_matrix.build(&dist, sum_matrix.ncol(), sum_values, sum_cols, sum_row_start);

      // total_matrix.sparse_indexed_output(std::cout);

      DoubleVector direct_solution(dist);
      total_matrix.solve(rhs,direct_solution);
      // direct_solution.output(std::cout);


      // Try solving with GMRES:
      // ============================================================
      GMRES<CRDoubleMatrix> gmres_solver;
      gmres_solver.tolerance() = 1e-8;
      // gmres_solver.Setup_preconditioner_before_solve
      DoubleVector solve_solution(dist,1.0);
      gmres_solver.solve(&sum_matrix, rhs, solve_solution);

      // solve_solution.output(std::cout);

      // Compare results
      // ============================================================
      DoubleVector diff(solve_solution);
      diff -= direct_solution;
      double a = diff.max();
      std::cout << "\nmax error is: " <<  a << std::endl;

      norm_error.push_back(a);
    }

  double mean_norm_error = std::accumulate(norm_error.begin(), norm_error.end(), 0.0)
    / norm_error.size();
  std::cout << std::endl << std::endl
            << "mean(max errors) is: " << mean_norm_error << std::endl;

  if(mean_norm_error > 1e-6)
    {
      return 1;
    }
  return 0;
}
