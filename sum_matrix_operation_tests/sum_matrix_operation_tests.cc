/*
  description of file goes here
*/

#include "generic.h"
#include "../my_general_header.h"

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
  const unsigned mat_size = 5;
  const unsigned n_sparse_values = 10;
  const unsigned max_val = 10000;
  srand(22);

  // A linear algebra "distribution": specifies no inter-processor communicator
  // and mat_size rows.
  OomphCommunicator communicator;
  LinearAlgebraDistribution dist(&communicator,mat_size,false);

  // Create empty sum of matrices matrix
  SumOfMatrices sum_matrix;

  // A random dense matrix
  DenseDoubleMatrix main_matrix(mat_size,mat_size,0.0);
  for(unsigned i=0; i<mat_size; i++)
    for(unsigned j=0; j<mat_size; j++)
      main_matrix(i,j) = rand() % max_val;

  // Add some diagonal terms to make it well conditioned
  for(unsigned i=0; i<mat_size; i++)
    main_matrix(i,i) += 2*max_val;

  sum_matrix.main_matrix_pt() = &main_matrix;


  // A random compressed row form sparse matrix ??ds note that this will probably
  // give overlapping values, slightly dodgy but works for now.
  Vector<int> cr_row_index(n_sparse_values), cr_col_index(n_sparse_values);
  Vector<double> cr_values(n_sparse_values);

  for(unsigned i=0; i<n_sparse_values; i++)
    {
      cr_row_index[i] = (rand() % mat_size);
      cr_col_index[i] = (rand() % mat_size);
      cr_values[i] = (rand() % max_val);
    }

  std::sort(cr_row_index.begin(),cr_row_index.end());

  std::cout << cr_row_index << std::endl;
  std::cout << cr_col_index << std::endl;
  std::cout << cr_values << std::endl;

  CRDoubleMatrix a_cr_matrix(cr_values,cr_row_index,cr_col_index,mat_size,mat_size);
  std::map<long unsigned, long unsigned> row_map, col_map;
  identity_map(row_map, mat_size); identity_map(col_map, mat_size);

  std::cout << col_map << std::endl;


  sum_matrix.add_matrix(&a_cr_matrix, &row_map, &col_map);
  // a_cr_matrix.sparse_indexed_output(std::cout);


  // A random rhs vector
  DoubleVector rhs;
  rhs.build(dist,1.0);
  for(unsigned i=0; i<mat_size; i++) rhs[i] = rand() % max_val;



  // Make a total matrix to compare
  ////////////////////////////////////////////////////////////
  DenseDoubleMatrix total_matrix(mat_size,mat_size,0.0);
  for(unsigned i=0; i<mat_size; i++)
    for(unsigned j=0; j<mat_size; j++)
      total_matrix(i,j) = main_matrix(i,j); // manual copy..

  // Add on sparse matrix values
  for(unsigned i=0; i<n_sparse_values; i++)
    {
      unsigned row = cr_row_index[i];
      unsigned col = cr_col_index[i];
      total_matrix(row,col) += cr_values[i];
    }

  main_matrix.output(std::cout);
  a_cr_matrix.sparse_indexed_output(std::cout);
  total_matrix.output(std::cout);

  DoubleVector direct_solution(dist);
  total_matrix.solve(rhs,direct_solution);
  direct_solution.output(std::cout);


  // Try solving with GMRES:
  //////////////////////////////////////////////////////

  GMRES<CRDoubleMatrix> gmres_solver;
  // gmres_solver.Setup_preconditioner_before_solve
  DoubleVector solve_solution(dist,1.0);
  gmres_solver.solve(&sum_matrix, rhs, solve_solution);

  solve_solution.output(std::cout);


  DoubleVector diff(solve_solution);
  diff -= direct_solution;
  std::cout << "\nError is:" << std::endl;
  diff.output(std::cout);


  return 0;
}

