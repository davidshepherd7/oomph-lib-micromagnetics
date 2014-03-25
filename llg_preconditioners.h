#ifndef OOMPH_LLG_PRECONDITIONERS_H
#define OOMPH_LLG_PRECONDITIONERS_H


#include "../../src/generic/preconditioner.h"
#include "../../src/generic/linear_solver.h"
#include "../../src/generic/general_purpose_preconditioners.h"
#include "../../src/generic/general_purpose_block_preconditioners.h"

#include "micromag_factories.h"

namespace oomph
{

  class LLGBlockPreconditioner : public BlockTriangularPreconditioner<CRDoubleMatrix>
  {
public:
    LLGBlockPreconditioner() {}


    void build(bool exact_phi1=false, bool exact_phi=false,
               bool exact_llg=true, bool schur_complement_llg=false)
    {
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
      if(schur_complement_llg)
        {
          throw OomphLibError("Not yet implemented",
                              OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
        }
      if(exact_llg)
        {
          llg_phi_bulk_prec_pt->set_subsidiary_preconditioner_pt
            (new SuperLUPreconditioner, 0);
        }
      else
        {
          llg_phi_bulk_prec_pt->set_subsidiary_preconditioner_pt
            (Factories::preconditioner_factory("ilu-2"), 0);
        }

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

    ~LLGBlockPreconditioner() {}

    void setup()
    {
      BlockTriangularPreconditioner<CRDoubleMatrix>::setup();
    }
  };


} // End of oomph namespace

#endif
