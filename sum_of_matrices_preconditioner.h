#ifndef OOMPH_SUM_OF_MATRICES_PRECONDITIONER_H
#define OOMPH_SUM_OF_MATRICES_PRECONDITIONER_H

#include "../../src/generic/preconditioner.h"

using namespace oomph;

namespace oomph
{

  /// Preconditioner which wraps around an underlying preconditioner but
  /// only applies it to the main matrix of a sum of matrices.
  class MainMatrixOnlyPreconditioner : public Preconditioner
  {
  public:

    MainMatrixOnlyPreconditioner(Preconditioner* underlying_preconditioner_pt)
    {
      Underlying_prec_pt = underlying_preconditioner_pt;
    }

    /// Clean up the real preconditioner
    virtual ~MainMatrixOnlyPreconditioner()
    {
      delete Underlying_prec_pt;
      Underlying_prec_pt = 0;
    }

    void preconditioner_solve(const DoubleVector &r, DoubleVector &z)
      {
        Underlying_prec_pt->preconditioner_solve(r, z);
      }

    void setup(SumOfMatrices* sum_matrix_pt,
               const OomphCommunicator* comm_pt = 0)
    {
      Underlying_prec_pt->setup(sum_matrix_pt->main_matrix_pt(),
                                          comm_pt);
    }
    void setup() {Underlying_prec_pt->setup();}

    void clean_up_memory() {Underlying_prec_pt->clean_up_memory();}

    DoubleMatrixBase* matrix_pt() const
    {return Underlying_prec_pt->matrix_pt();}

    /// If we get a sum of matrices pointer then pass down the main matrix
    /// pt only, otherwise just pass down the matrix pt.
    void set_matrix_pt(DoubleMatrixBase* matrix_pt)
    {
      SumOfMatrices* s_matrix_pt = dynamic_cast<SumOfMatrices*>(matrix_pt);
      if(s_matrix_pt)
        {
          Underlying_prec_pt->set_matrix_pt(s_matrix_pt->main_matrix_pt());
        }
      else
        {
          Underlying_prec_pt->set_matrix_pt(matrix_pt);
        }
    }
    void set_matrix_pt(SumOfMatrices* sum_matrix_pt)
    {Underlying_prec_pt->set_matrix_pt(sum_matrix_pt->main_matrix_pt());}

    const OomphCommunicator* comm_pt() const
    {return Underlying_prec_pt->comm_pt();}
    void set_comm_pt(const OomphCommunicator* const comm_pt)
    {Underlying_prec_pt->set_comm_pt(comm_pt);}

  private:

    Preconditioner* Underlying_prec_pt;

  };

} // End of oomph namespace

#endif
