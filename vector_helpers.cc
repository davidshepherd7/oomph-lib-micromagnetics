#include "vector_helpers.h"
#include "../../src/generic/linear_solver.h"

using namespace oomph;

namespace VectorOps
{

  void multiple_rhs_solve_hack(CRDoubleMatrix& A,
                               DoubleMatrixBase& B,
                               DenseDoubleMatrix& X)
  {
    // Create factorisation of A
    SuperLUSolver S;
    S.factorise(&A);

    // For each col in B do a backsub and put results into X
    for(unsigned j=0; j<B.ncol(); j++)
      {
        DoubleVector b(A.distribution_pt()), x(A.distribution_pt());

        // Copy to vector
        for(unsigned i=0; i<B.nrow(); i++) {b[i] = B(i, j);}

        S.backsub(b, x);

        // Copy from vector
        for(unsigned i=0; i<B.nrow(); i++) {X(i, j) = x[i];}
      }
  }

}
