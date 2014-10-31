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
      underlying_prec_pt()->setup(sum_matrix_pt->main_matrix_pt());
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


  /// Preconditioner which merges the diagonals of the added matrices into
  /// a new main matrix for preconditioning.
  class MainMatrixAndDiagsPreconditioner : public SoMPreconditioner
  {
  public:

    MainMatrixAndDiagsPreconditioner()
    {
      Preconditioner_matrix_pt = 0;
    }

    virtual ~MainMatrixAndDiagsPreconditioner()
    {
      delete Preconditioner_matrix_pt; Preconditioner_matrix_pt = 0;
    }

  protected:

    /// ??ds
    void set_preconditioner_matrix_pt(SumOfMatrices* s_matrix_pt)
    {
      CRDoubleMatrix* main_pt = checked_dynamic_cast<CRDoubleMatrix*>
        (s_matrix_pt->main_matrix_pt());

      // Copy main matrix into the new one
      Preconditioner_matrix_pt = new CRDoubleMatrix(*main_pt);

      for(unsigned k=0; k<s_matrix_pt->n_added_matrix(); k++)
        {
          DoubleMatrixBase* added_pt = s_matrix_pt->added_matrix_pt(k);

          // Create diagonal CR version of added matrix 1
          std::list<VectorOps::RowColVal> rcvs;
          for(unsigned j=0; j<added_pt->nrow(); j++)
            {
              VectorOps::RowColVal rcv;
              rcv.row = s_matrix_pt->row_map_pt(k)->added_to_main(j);
              rcv.col = s_matrix_pt->col_map_pt(k)->added_to_main(j);
              rcv.val = (*added_pt)(j, j);

              rcvs.push_back(rcv);
            }

          // Build the cr matrix of the diag of this added matrix
          CRDoubleMatrix diag_added_matrix;
          VectorOps::rowcolvals_to_crmatrix(rcvs, main_pt->distribution_pt(),
                                            main_pt->ncol(),
                                            diag_added_matrix);

          // Add to the total matrix
          VectorOps::cr_matrix_add(*Preconditioner_matrix_pt, diag_added_matrix,
                                   *Preconditioner_matrix_pt);
        }

      // Finally assign it in the preconditioner
      Underlying_prec_pt->set_matrix_pt(Preconditioner_matrix_pt);
    }

  private:

    CRDoubleMatrix* Preconditioner_matrix_pt;

    /// Broken copy constructor
    MainMatrixAndDiagsPreconditioner(const MainMatrixAndDiagsPreconditioner& dummy)
    {BrokenCopy::broken_copy("MainMatrixAndDiagsPreconditioner");}

    /// Broken assignment operator
    void operator=(const MainMatrixAndDiagsPreconditioner& dummy)
    {BrokenCopy::broken_assign("MainMatrixAndDiagsPreconditioner");}
  };


  /// Given a preconditioner:
  /// 1) if it's a the right type of preconditioner return it
  /// 2) otherwise if its a SoM preconditioner containing the right type
  ///    of preconditioner then return a pointer to the underlying
  ///    preconditioner.
  /// 3) otherwise return null
  template<class T>
  T smart_cast_preconditioner(Preconditioner* prec_pt)
  {
    T bp_pt = dynamic_cast<T> (prec_pt);
    if(bp_pt != 0)
      {
        return bp_pt;
      }
    else if(dynamic_cast<SoMPreconditioner*>(prec_pt) != 0)
      {
        SoMPreconditioner* som_main_prec_pt
          = dynamic_cast<SoMPreconditioner*>(prec_pt);

        return dynamic_cast<T>
          (som_main_prec_pt->underlying_prec_pt());
      }
    else
      {
        return 0;
      }
  }


} // End of oomph namespace

#endif
