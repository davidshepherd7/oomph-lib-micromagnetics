/*
  description of file goes here
*/

#include "generic.h"
#include "../my_general_header.h"
#include <numeric>

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
  const unsigned mat_size = 100;
  const unsigned n_sparse_values = mat_size;
  const unsigned dense_block_size = mat_size/10;
  const unsigned max_val = 100;
  const unsigned n_tests = 20;

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

      // A random compressed row form sparse matrix
       // to avoid overlaps we have exactly one element per row
      Vector<int> main_row_index(n_sparse_values), main_col_index(n_sparse_values);
      Vector<double> main_values(n_sparse_values);

      for(unsigned i=0; i<n_sparse_values; i++)
	{
	  main_row_index[i] = i;
	  main_col_index[i] = (rand() % mat_size);
	  main_values[i] = (rand() % max_val);
	}

      std::sort(main_row_index.begin(),main_row_index.end());

      // std::cout << main_row_index << std::endl;
      // std::cout << main_col_index << std::endl;
      // std::cout << main_values << std::endl;

      CRDoubleMatrix main_matrix(main_row_index,main_col_index,main_values,mat_size,mat_size);
      sum_matrix.main_matrix_pt() = &main_matrix;

      // Add a diagonal sparse matrix
      Vector<int> id_row_index(mat_size), id_col_index(mat_size);
      Vector<double> id_values(mat_size,2*max_val);

      for(unsigned i=0; i<mat_size; i++)
	{
	  id_row_index[i] = i;
	  id_col_index[i] = i;
	}

      CRDoubleMatrix id_matrix(id_row_index,id_col_index,id_values,mat_size,mat_size);
      std::map<long unsigned, long unsigned> row_map, col_map;
      identity_map(row_map, mat_size); identity_map(col_map, mat_size);
      sum_matrix.add_matrix(&id_matrix,&row_map,&col_map,0);

      // A random dense matrix
      DenseDoubleMatrix dense_block(dense_block_size,dense_block_size,0.0);
      for(unsigned i=0; i<dense_block_size; i++)
      	for(unsigned j=0; j<dense_block_size; j++)
      	  dense_block(i,j) = rand() % max_val;

      std::map<long unsigned, long unsigned> dense_row_map, dense_col_map;
      for(unsigned j=0; j<dense_block_size; j++)
	{
	  // Get random row and col which are not already used
	  unsigned new_rand_row, new_rand_col;
	  do
	    {
	      new_rand_row = rand() % mat_size;
	    }
	  while(dense_row_map.find(new_rand_row) != dense_row_map.end());

	  do
	    {
	      new_rand_col = rand() % mat_size;
	    }
	  while(dense_col_map.find(new_rand_col) != dense_col_map.end());

	  std::pair<long unsigned, long unsigned>
	    row_map_pair(new_rand_row,j), col_map_pair(new_rand_col,j);
	  dense_row_map.insert(row_map_pair);
	  dense_col_map.insert(col_map_pair);
	}

      sum_matrix.add_matrix(&dense_block, &dense_row_map, &dense_col_map);

      // A random rhs vector
      DoubleVector rhs;
      rhs.build(dist,1.0);
      for(unsigned i=0; i<mat_size; i++) rhs[i] = rand() % max_val;



      // Make a total matrix to compare
      ////////////////////////////////////////////////////////////
      Vector<int> sum_rows, sum_cols;
      Vector<double> sum_values;
      sum_matrix.get_as_indicies(sum_rows,sum_cols,sum_values);
      CRDoubleMatrix total_matrix;
      total_matrix.build(&dist,sum_rows,sum_cols,sum_values,
			 sum_matrix.nrow(),sum_matrix.ncol());

      DoubleVector direct_solution(dist);
      total_matrix.solve(rhs,direct_solution);
      // direct_solution.output(std::cout);


      // Try solving with GMRES:
      //////////////////////////////////////////////////////

      GMRES<CRDoubleMatrix> gmres_solver;
      gmres_solver.tolerance() = 1e-8;
      // gmres_solver.Setup_preconditioner_before_solve
      DoubleVector solve_solution(dist,1.0);
      gmres_solver.solve(&sum_matrix, rhs, solve_solution);

      // solve_solution.output(std::cout);


      DoubleVector diff(solve_solution);
      diff -= direct_solution;
      double a = diff.max();
      std::cout << "\nmax error is: " <<  a << std::endl;

      norm_error.push_back(a);
    }

  double mean_norm_error = std::accumulate(norm_error.begin(), norm_error.end(), 0.0)
    / norm_error.size();
  std::cout << std::endl << std::endl << "mean(max error) is:" << mean_norm_error << std::endl;

  return 0;
}

