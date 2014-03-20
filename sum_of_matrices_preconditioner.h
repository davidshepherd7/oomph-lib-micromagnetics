#ifndef OOMPH_SUM_OF_MATRICES_PRECONDITIONER_H
#define OOMPH_SUM_OF_MATRICES_PRECONDITIONER_H

#include "../../src/generic/preconditioner.h"
#include "../../src/generic/sum_of_matrices.h"

#include "vector_helpers.h"


using namespace oomph;

namespace oomph
{

  ///
  class SoMPreconditioner : public Preconditioner
  {
  public:
    /// Constructor
    SoMPreconditioner()
    {
      Underlying_prec_pt = 0;
    }

    /// Virtual destructor, clean up the real preconditioner.
    virtual ~SoMPreconditioner()
    {
      delete Underlying_prec_pt;
      Underlying_prec_pt = 0;
    }

    DoubleMatrixBase* matrix_pt() const
    {return Underlying_prec_pt->matrix_pt();}

    /// If we get a sum of matrices pointer then pass down the main matrix
    /// pt only, otherwise just pass down the matrix pt.
    void set_matrix_pt(DoubleMatrixBase* matrix_pt)
    {
      SumOfMatrices* s_matrix_pt = dynamic_cast<SumOfMatrices*>(matrix_pt);
      if(s_matrix_pt)
        {
          set_preconditioner_matrix_pt(s_matrix_pt);
        }
      else
        {
          Underlying_prec_pt->set_matrix_pt(matrix_pt);
        }
    }
    void set_matrix_pt(SumOfMatrices* sum_matrix_pt)
    {set_preconditioner_matrix_pt(sum_matrix_pt);}

    const OomphCommunicator* comm_pt() const
    {return Underlying_prec_pt->comm_pt();}
    void set_comm_pt(const OomphCommunicator* const comm_pt)
    {Underlying_prec_pt->set_comm_pt(comm_pt);}

    Preconditioner* underlying_prec_pt() const
    {
#ifdef PARANOID
      if(Underlying_prec_pt == 0)
        {
          std::string err = "Underlying_prec_pt is null!";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif
      return Underlying_prec_pt;
    }

    void set_underlying_prec_pt(Preconditioner* underlying_prec_pt)
    {Underlying_prec_pt = underlying_prec_pt;}

    void preconditioner_solve(const DoubleVector &r, DoubleVector &z)
    {
      underlying_prec_pt()->preconditioner_solve(r, z);
    }

    void setup(SumOfMatrices* sum_matrix_pt,
               const OomphCommunicator* comm_pt = 0)
    {
      underlying_prec_pt()->setup(sum_matrix_pt->main_matrix_pt(),
                                  comm_pt);
    }
    void setup() {underlying_prec_pt()->setup();}

    virtual void clean_up_memory() {underlying_prec_pt()->clean_up_memory();}


  protected:
    /// The preconditioner to use on the main matrix
    Preconditioner* Underlying_prec_pt;

    /// Set the appropriate matrix, constructed from the sum matrix, in the
    /// underlying preconditioner.
    virtual void set_preconditioner_matrix_pt(SumOfMatrices* s_matrix_pt) = 0;

  private:

    /// Broken copy constructor
    SoMPreconditioner(const SoMPreconditioner& dummy)
    {BrokenCopy::broken_copy("SoMPreconditioner");}

    /// Broken assignment operator
    void operator=(const SoMPreconditioner& dummy)
    {BrokenCopy::broken_assign("SoMPreconditioner");}

  };


  /// Preconditioner which wraps around an underlying preconditioner but
  /// only applies it to the main matrix of a sum of matrices.
  class MainMatrixOnlyPreconditioner : public SoMPreconditioner
  {
  public:

    MainMatrixOnlyPreconditioner() {}

    virtual ~MainMatrixOnlyPreconditioner() {}

    /// Just use the main matrix
    void set_preconditioner_matrix_pt(SumOfMatrices* s_matrix_pt)
      {
        Underlying_prec_pt->set_matrix_pt(s_matrix_pt->main_matrix_pt());
      }

  private:

    /// Broken copy constructor
    MainMatrixOnlyPreconditioner(const MainMatrixOnlyPreconditioner& dummy)
    {BrokenCopy::broken_copy("MainMatrixOnlyPreconditioner");}

    /// Broken assignment operator
    void operator=(const MainMatrixOnlyPreconditioner& dummy)
    {BrokenCopy::broken_assign("MainMatrixOnlyPreconditioner");}
  };


} // End of oomph namespace

#endif
