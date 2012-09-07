#ifndef OOMPH_MICROMAGNETICS_PRECONDITIONERS_H
#define OOMPH_MICROMAGNETICS_PRECONDITIONERS_H

/*
  description of file goes here
*/

#include "generic.h"

using namespace oomph;

namespace oomph
{

  // ============================================================
  /// A very simple preconditioner... probably no use for most cases.
  // ============================================================
  class DiagonalPreconditioner : public Preconditioner
  {
  public:

    DiagonalPreconditioner() {}

    ~DiagonalPreconditioner() {}

    /// \short Apply the preconditioner. This method should apply the
    /// preconditioner operator to the vector r and return the vector z.
    void preconditioner_solve(const DoubleVector &r, DoubleVector &z)
    {
      unsigned nrow = inverse_diagonal_values.nrow();

#ifdef PARANOID
      if(r.nrow() != nrow)
	{
	  std::ostringstream error_msg;
	  error_msg << "Different numbers of rows in diagonal matrix and input vector ("
		    << nrow << " vs "
		    << r.nrow() << ").";
	  throw OomphLibError(error_msg.str(),
			      "DiagonalPreconditioner::preconditioner_solve",
			      OOMPH_EXCEPTION_LOCATION);
	}
#endif

      // (inverted) diagonal matrix * vector multiplication is very simple:
      for(unsigned j=0; j<nrow; j++)
	z[j] = r[j] * inverse_diagonal_values[j];
    }

    /// \short Setup the preconditioner (emtpy).
    void setup(Problem* problem_pt, DoubleMatrixBase* matrix_pt)
    {
      // Get distribution info from the Jacobian
      unsigned nrow = matrix_pt->nrow();

      // Set up the inverse_diagonal_values vector
      //??dsparallel - probably not the right dist pt
      inverse_diagonal_values.build(problem_pt->dof_distribution_pt(),0.0);

      // Fill in the inverse_diagonal_values vector
      for(unsigned j=0; j<nrow; j++)
	inverse_diagonal_values[j] = 1/(matrix_pt->operator()(j,j));

      inverse_diagonal_values.output("preconditioner_values");
    }

    /// \short Clean up memory (empty).
    void clean_up_memory(){};

  private:

    /// Store D^{-1} where D is just the diagonal entries of the Jacobian.
    DoubleVector inverse_diagonal_values;

  };

  // ============================================================
  ///
  // ============================================================
  template<typename MAIN_MATRIX_TYPE>
  class MicromagneticsBlockDiagonalPreconditioner :
    public BlockDiagonalPreconditioner<MAIN_MATRIX_TYPE>
  {

    /// \short Setup the preconditioner
    void setup(Problem* problem_pt,
	       DoubleMatrixBase* matrix_pt)
    {
      std::cout << "Calling BlockDiagonalPreconditioner but ignoring the dense sub-block."
		<< std::endl;
      SumOfMatrices* sum_matrix_pt =
	dynamic_cast<SumOfMatrices*>(matrix_pt);

      // Set up everything ignoring the added matrices
      //??ds no need to do explicit down-cast right?
      BlockDiagonalPreconditioner<MAIN_MATRIX_TYPE>::
	setup(problem_pt, sum_matrix_pt->main_matrix_pt());
    }

  };

  // // ============================================================
  // ///
  // // ============================================================
  // class SomePreconditioner : public Preconditioner
  // {
  // public:

  //   SomePreconditioner() {}

  //   ~SomePreconditioner() {}

  //   /// \short Apply the preconditioner. This method should apply the
  //   /// preconditioner operator to the vector r and return the vector z.
  //   void preconditioner_solve(const DoubleVector &r, DoubleVector &z)
  //   {

  //   }

  //   /// \short Setup the preconditioner (emtpy).
  //   void setup(Problem* problem_pt, DoubleMatrixBase* matrix_pt)
  //   {

  //   }

  //   /// \short Clean up memory (empty).
  //   void clean_up_memory(){};

  // };


} // End of oomph namespace

#endif
