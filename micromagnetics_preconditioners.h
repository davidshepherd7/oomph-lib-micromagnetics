#ifndef OOMPH_MICROMAGNETICS_PRECONDITIONERS_H
#define OOMPH_MICROMAGNETICS_PRECONDITIONERS_H

/*
  description of file goes here
*/

#include "generic.h"

using namespace oomph;

namespace oomph
{

  //??ds these classes should maybe just be class specialisations to SumOfMatrices?

  // then again maybe not - not all SUmOfMatrices will have ignoreable extra matrices

  // ============================================================
  ///
  // ============================================================
  template<typename MAIN_MATRIX_TYPE>
  class MMBlockDiagonalPreconditioner :
    public BlockDiagonalPreconditioner<MAIN_MATRIX_TYPE>
  {
  public:

    //??ds add broken constructors etc.

    /// \short Setup the preconditioner
    void setup(Problem* problem_pt,
               DoubleMatrixBase* matrix_pt)
    {
      std::cout << "Calling BlockDiagonalPreconditioner but ignoring the dense sub-block."
                << std::endl;
      SumOfMatrices* sum_matrix_pt =
        dynamic_cast<SumOfMatrices*>(matrix_pt);

      // Set up everything ignoring the added matrices
      BlockDiagonalPreconditioner<MAIN_MATRIX_TYPE>::
        setup(problem_pt, sum_matrix_pt->main_matrix_pt());
    }

  };

  // ============================================================
  ///
  // ============================================================
  template<typename MAIN_MATRIX_TYPE>
  class MMBlockTriangularPreconditioner :
    public BlockTriangularPreconditioner<MAIN_MATRIX_TYPE>
  {
  public:

    //??ds add broken constructors etc.

    /// \short Setup the preconditioner
    void setup(Problem* problem_pt,
               DoubleMatrixBase* matrix_pt)
    {
      std::cout << "Calling BlockTriangularPreconditioner but ignoring the dense sub-block."
                << std::endl;
      SumOfMatrices* sum_matrix_pt =
        dynamic_cast<SumOfMatrices*>(matrix_pt);

      // Set up everything ignoring the added matrices
      BlockTriangularPreconditioner<MAIN_MATRIX_TYPE>::
        setup(problem_pt, sum_matrix_pt->main_matrix_pt());
    }

    // void block_jacobian_output(Problem* problem_pt, DoubleMatrixBase* matrix_pt,
    //                           const std::string& dirname)
    // {
    //   SumOfMatrices* sum_matrix_pt =
    //    dynamic_cast<SumOfMatrices*>(matrix_pt);

    //   // Set up everything ignoring the added matrices
    //   BlockTriangularPreconditioner<MAIN_MATRIX_TYPE>::
    //    block_jacobian_output(problem_pt, sum_matrix_pt->main_matrix_pt(), dirname);
    // }

  };

  // ============================================================
  ///
  // ============================================================
  template<typename MAIN_MATRIX_TYPE>
  class MMExactBlockPreconditioner :
    public ExactBlockPreconditioner<MAIN_MATRIX_TYPE>
  {
  public:

    //??ds add broken constructors etc.

    /// \short Setup the preconditioner
    void setup(Problem* problem_pt,
               DoubleMatrixBase* matrix_pt)
    {
      std::cout << "Calling ExactBlockPreconditioner but ignoring the dense sub-block."
                << std::endl;
      SumOfMatrices* sum_matrix_pt =
        dynamic_cast<SumOfMatrices*>(matrix_pt);

      // Set up everything ignoring the added matrices
      ExactBlockPreconditioner<MAIN_MATRIX_TYPE>::
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
