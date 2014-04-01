#ifndef OOMPH_LLG_PRECONDITIONERS_H
#define OOMPH_LLG_PRECONDITIONERS_H


#include "../../src/generic/preconditioner.h"
#include "../../src/generic/linear_solver.h"
#include "../../src/generic/general_purpose_preconditioners.h"
#include "../../src/generic/general_purpose_block_preconditioners.h"

#include "vector_helpers.h"
#include "oomph_factories.h"
#include "llg_factories.h"

namespace oomph
{
  using namespace VectorOps;

  class MagnetostaticsPreconditioner :
    public BlockTriangularPreconditioner<CRDoubleMatrix>
  {
public:
    MagnetostaticsPreconditioner()
    {
      Llg_preconditioner_pt = 0;
    }

    /// Virtual destructor. Everything is deleted (recursively) by the
    /// destructors in general purpose block preconditioners.
    virtual ~MagnetostaticsPreconditioner() {}

    Preconditioner* Llg_preconditioner_pt;


    void build(bool exact_phi1=false, bool exact_phi=false)
    {
      //??ds use exact by default for now, remove later?
      if(Llg_preconditioner_pt == 0)
        {
          Llg_preconditioner_pt = new SuperLUPreconditioner;
        }

      {
        // this->upper_triangular();
        this->lower_triangular(); //??ds why lower

        Vector<unsigned> dof_to_block(7);

        // phi
        dof_to_block[0] = 0; // phi bulk
        dof_to_block[5] = 0; // phi bound

        // m
        dof_to_block[2] = 0;
        dof_to_block[3] = 0;
        dof_to_block[4] = 0;

        // phi1
        dof_to_block[1] = 1; // phi1 bulk
        dof_to_block[6] = 1; // phi1 bound

        set_dof_to_block_map(dof_to_block);
      }


      // phi1 block
      // ============================================================
      if(exact_phi1)
        {
          set_subsidiary_preconditioner_pt(new SuperLUPreconditioner, 1);
        }
      else
        {
          set_subsidiary_preconditioner_pt
            (Factories::preconditioner_factory("poisson-amg"), 1);
        }


      // llg + phi block
      // ============================================================
      BlockTriangularPreconditioner<CRDoubleMatrix>* llg_phi_prec_pt = 0;
      {
        Vector<unsigned> master_to_subs_block_map(5);
        master_to_subs_block_map[0] = 2; // m first
        master_to_subs_block_map[1] = 3;
        master_to_subs_block_map[2] = 4;
        master_to_subs_block_map[3] = 0; // Then phi
        master_to_subs_block_map[4] = 5;


        Vector<unsigned> dof_to_block(5);
        // m + phi bulk in block 0
        dof_to_block[0] = 0;
        dof_to_block[1] = 0;
        dof_to_block[2] = 0;
        dof_to_block[3] = 0;

        // phi boundary values in block 1
        dof_to_block[4] = 1;


        llg_phi_prec_pt = new BlockTriangularPreconditioner<CRDoubleMatrix>;
        llg_phi_prec_pt->set_dof_to_block_map(dof_to_block);
        // llg_phi_prec_pt->lower_triangular(); //??ds should be lower, ordering weird...

        llg_phi_prec_pt->turn_into_subsidiary_block_preconditioner
          (this, master_to_subs_block_map);
      }
      set_subsidiary_preconditioner_pt(llg_phi_prec_pt, 0);


      // phi boundary block (just -I)
      // ----------------------------------------------------------------
      llg_phi_prec_pt->set_subsidiary_preconditioner_pt
        (new SuperLUPreconditioner, 1);


      // llg + phi bulk block
      // ----------------------------------------------------------------
      BlockTriangularPreconditioner<CRDoubleMatrix>* llg_phi_bulk_prec_pt = 0;
      {
        Vector<unsigned> master_to_subs_block_map(4);
        master_to_subs_block_map[0] = 0; // m first
        master_to_subs_block_map[1] = 1;
        master_to_subs_block_map[2] = 2;
        master_to_subs_block_map[3] = 3; // Then phi bulk

        Vector<unsigned> dof_to_block(4);
        dof_to_block[0] = 0;   // m in block 0
        dof_to_block[1] = 0;
        dof_to_block[2] = 0;
        dof_to_block[3] = 1;    // phi in block 1


        llg_phi_bulk_prec_pt = new BlockTriangularPreconditioner<CRDoubleMatrix>;
        llg_phi_bulk_prec_pt->set_dof_to_block_map(dof_to_block);
        // llg_phi_bulk_prec_pt->lower_triangular(); //??ds should be upper, ordering weird...

        llg_phi_bulk_prec_pt->turn_into_subsidiary_block_preconditioner
          (llg_phi_prec_pt, master_to_subs_block_map);
      }
      llg_phi_prec_pt->set_subsidiary_preconditioner_pt
        (llg_phi_bulk_prec_pt, 0);

      // llg block
      llg_phi_bulk_prec_pt->set_subsidiary_preconditioner_pt(Llg_preconditioner_pt, 0);

      // phi bulk block
      if(exact_phi)
        {
          llg_phi_bulk_prec_pt->set_subsidiary_preconditioner_pt
            (new SuperLUPreconditioner, 1);
        }
      else
        {
          llg_phi_bulk_prec_pt->set_subsidiary_preconditioner_pt
            (Factories::preconditioner_factory("poisson-amg"), 1);
        }

    }

    void setup()
    {
      BlockTriangularPreconditioner<CRDoubleMatrix>::setup();
    }
  };



  /// ??ds
  class LLGBlockPreconditioner :
    public BlockTriangularPreconditioner<CRDoubleMatrix>
  {
  public:
    /// Constructor
    LLGBlockPreconditioner()
    {
      Include_jcc = false;
      Use_schur_complement = false;
    }

    /// Virtual destructor
    virtual ~LLGBlockPreconditioner() {}

    void build()
    {
      this->upper_triangular();

      // Group first two directions of m into one block.
      Vector<unsigned> dof_to_block(7);

      // ??ds modify to accept only 3 blocks (the 3 m)
      // dof_to_block[0] = 0;
      // dof_to_block[1] = 0;
      // dof_to_block[2] = 1;

      // m
      dof_to_block[2] = 0;
      dof_to_block[3] = 0;
      dof_to_block[4] = 1;

      // phis
      dof_to_block[0] = 2; // phi bulk
      dof_to_block[5] = 2; // phi bound
      dof_to_block[1] = 2; // phi1 bulk
      dof_to_block[6] = 2; // phi1 bound

      set_dof_to_block_map(dof_to_block);

      // m_aabb block exact
      set_subsidiary_preconditioner_pt(new SuperLUPreconditioner, 0);

      // phi blocks: nothing here I hope!
      set_subsidiary_preconditioner_pt(new IdentityPreconditioner, 2);

      // m_cc block: Use an exact solve, but the block is replaced by its
      // Schur complement during setup.
      set_subsidiary_preconditioner_pt(new SuperLUPreconditioner, 1);
    }

    // This is inefficient all over the place, don't use it for real!
    void setup()
    {
      if(Use_schur_complement)
        {
          // Construct Schur complement:
          CRDoubleMatrix A = get_block(0, 0);
          CRDoubleMatrix B = get_block(0, 1);
          DenseDoubleMatrix AinvB; multiple_rhs_solve_hack(A, B, AinvB);
          CRDoubleMatrix C = get_block(1, 0);

          // Convert to a cr matrix because multiply needs it...
          CRDoubleMatrix AinvB_cr;
          std::list<RowColVal> rcv = get_as_indicies(AinvB);
          rowcolvals_to_crmatrix(rcv,
                                 block_distribution_pt(1),
                                 block_distribution_pt(1)->nrow(),
                                 AinvB_cr);

          CRDoubleMatrix CAinvB;
          C.multiply(AinvB_cr, CAinvB);

          // multiply by -1
          for(unsigned j=0; j<CAinvB.nnz(); j++) {CAinvB.value()[j] *= -1;}

          // Add J_cc block if requested
          if(Include_jcc)
            {
              CRDoubleMatrix D = get_block(1,1);
              cr_matrix_add(D, CAinvB, CAinvB);
            }

          // Replace the block with the Schur complement
          // this->set_replace_block(1,1, &CAinvB);
          throw OomphLibError("Not implemented (yet?).", OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      // Do the rest of the setup
      BlockTriangularPreconditioner<CRDoubleMatrix>::setup();
    }

    bool Include_jcc;
    bool Use_schur_complement;

  private:
    /// Broken copy constructor
    LLGBlockPreconditioner(const LLGBlockPreconditioner& dummy)
    {BrokenCopy::broken_copy("LLGBlockPreconditioner");}

    /// Broken assignment operator
    void operator=(const LLGBlockPreconditioner& dummy)
    {BrokenCopy::broken_assign("LLGBlockPreconditioner");}

  };



} // End of oomph namespace

#endif
