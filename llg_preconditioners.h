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

  class MagnetostaticsBlockPreconditioner :
    public BlockTriangularPreconditioner<CRDoubleMatrix>
  {
public:
    MagnetostaticsBlockPreconditioner()
    {
      Llg_preconditioner_pt = 0;
      Drop_P = false;
    }

    /// Virtual destructor. Everything is deleted (recursively) by the
    /// destructors in general purpose block preconditioners.
    virtual ~MagnetostaticsBlockPreconditioner() {}

    Preconditioner* Llg_preconditioner_pt;


    void build(bool exact_phi1=false, bool exact_phi=false)
    {
      //??ds use exact by default for now, remove later?
      if(Llg_preconditioner_pt == 0)
        {
          Llg_preconditioner_pt = new SuperLUPreconditioner;
        }

      {
        this->lower_triangular();

        Vector<unsigned> dof_to_block(7);

        // phi
        dof_to_block[0] = 1; // phi bulk
        dof_to_block[5] = 0; // phi bound

        // m
        dof_to_block[2] = 1;
        dof_to_block[3] = 1;
        dof_to_block[4] = 1;

        // phi1
        dof_to_block[1] = 2; // phi1 bulk
        dof_to_block[6] = 2; // phi1 bound

        set_dof_to_block_map(dof_to_block);
      }


      // phi boundary block (just -I)
      // ============================================================
      set_subsidiary_preconditioner_pt(new SuperLUPreconditioner, 0);


      // llg + phi bulk block
      // ============================================================
      BlockTriangularPreconditioner<CRDoubleMatrix>* llg_phi_bulk_prec_pt = 0;
      {
        Vector<unsigned> master_to_subs_block_map(4);
        master_to_subs_block_map[0] = 2; // m first
        master_to_subs_block_map[1] = 3;
        master_to_subs_block_map[2] = 4;
        master_to_subs_block_map[3] = 0; // Then phi bulk

        Vector<unsigned> dof_to_block(4);
        dof_to_block[0] = 0;   // m in block 0
        dof_to_block[1] = 0;
        dof_to_block[2] = 0;
        dof_to_block[3] = 1;   // phi in block 1


        llg_phi_bulk_prec_pt = new BlockTriangularPreconditioner<CRDoubleMatrix>;
        llg_phi_bulk_prec_pt->set_dof_to_block_map(dof_to_block);

        if(Drop_P)
          {
            llg_phi_bulk_prec_pt->lower_triangular();
          }
        else
          {
            llg_phi_bulk_prec_pt->upper_triangular();
          }

        llg_phi_bulk_prec_pt->turn_into_subsidiary_block_preconditioner
          (this, master_to_subs_block_map);


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

      // pass it to the main preconditioner
      set_subsidiary_preconditioner_pt(llg_phi_bulk_prec_pt, 1);



      // phi1 block
      // ============================================================
      if(exact_phi1)
        {
          set_subsidiary_preconditioner_pt(new SuperLUPreconditioner, 2);
        }
      else
        {
          set_subsidiary_preconditioner_pt
            (Factories::preconditioner_factory("poisson-amg"), 2);
        }

    }

    void setup()
    {
      BlockTriangularPreconditioner<CRDoubleMatrix>::setup();
    }

    bool Drop_P;
  };


  /// Wrapper class to allow llg block preconditioners to ignore the
  /// magnetostatics dofs without needing to modify the dof to block map.
  class DummyPinnedMsPreconditioner :
    public BlockPreconditioner<CRDoubleMatrix>
  {
  public:
    /// Constructor
    DummyPinnedMsPreconditioner()
    {
      Real_preconditioner_pt = 0;
    }

    /// Virtual destructor
    virtual ~DummyPinnedMsPreconditioner()
    {
      delete Real_preconditioner_pt; Real_preconditioner_pt = 0;
    }

    virtual void clean_up_memory()
    {
      if(Real_preconditioner_pt != 0)
        {
          Real_preconditioner_pt->clean_up_memory();
        }
    }

    void setup()
    {
      /// Call block setup because this is a master block preconditioner of
      /// "Real_preconditioner_pt".
      block_setup();

      // Check that we don't have any phi dofs
#ifdef PARANOID
      if(get_block(0, 0).nrow() != 0)
        {
          std::string err = "Non-empty phi block!";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      // Call setup for the real preconditioner
      Real_preconditioner_pt->setup(matrix_pt());
    }

    void preconditioner_solve(const DoubleVector& r, DoubleVector& z)
    {
      Real_preconditioner_pt->preconditioner_solve(r, z);
    }

    Preconditioner* Real_preconditioner_pt;

  private:
    /// Broken copy constructor
    DummyPinnedMsPreconditioner(const DummyPinnedMsPreconditioner& dummy)
    {BrokenCopy::broken_copy("DummyPinnedMsPreconditioner");}

    /// Broken assignment operator
    void operator=(const DummyPinnedMsPreconditioner& dummy)
    {BrokenCopy::broken_assign("DummyPinnedMsPreconditioner");}

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
      J_aabb_prec_pt = 0;
    }

    /// Virtual destructor
    virtual ~LLGBlockPreconditioner() {}

    void build()
    {
#ifdef PARANOID
      if(J_aabb_prec_pt == 0)
        {
          std::string err = "No preconditioner set for aabb block.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      this->upper_triangular();

      // Group first two directions of m into one block.
      Vector<unsigned> dof_to_block(3);
      dof_to_block[0] = 0;
      dof_to_block[1] = 0;
      dof_to_block[2] = 1;
      set_dof_to_block_map(dof_to_block);


      // J_aabb block
      // ============================================================
      set_subsidiary_preconditioner_pt(J_aabb_prec_pt, 0);


      // m_cc block
      // ============================================================

      // Use an exact solve, but the block is replaced by its Schur
      // complement during setup.
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

    Preconditioner* J_aabb_prec_pt;

  private:
    /// Broken copy constructor
    LLGBlockPreconditioner(const LLGBlockPreconditioner& dummy)
    {BrokenCopy::broken_copy("LLGBlockPreconditioner");}

    /// Broken assignment operator
    void operator=(const LLGBlockPreconditioner& dummy)
    {BrokenCopy::broken_assign("LLGBlockPreconditioner");}

  };



  /// Class to handle ??ds
  class LLGSubBlockPreconditioner
    : public BlockPreconditioner<CRDoubleMatrix>
  {
  public:
    /// Constructor
    LLGSubBlockPreconditioner()
    {
      Jaa_pt = 0;
      Prec1_pt = 0;
      Prec2_pt = 0;
    }

    /// Virtual destructor
    virtual ~LLGSubBlockPreconditioner()
    {
      clean_up_memory();
      delete Prec1_pt; Prec1_pt = 0;
      delete Prec2_pt; Prec2_pt = 0;
    }

    void clean_up_memory()
    {
      if(Prec1_pt != 0) {Prec1_pt->clean_up_memory();}
      if(Prec2_pt != 0) {Prec2_pt->clean_up_memory();}
      delete Jaa_pt; Jaa_pt = 0;
    }

    void build()
    {
      // J_aa block: exact solve
      Prec1_pt = new SuperLUPreconditioner;

      // J_bb block: Use an exact solve, but the block is replaced an
      // approximation to the Schur complement during setup.
      Prec2_pt = new SuperLUPreconditioner;
    }

    void setup()
    {
      // clean the memory
      this->clean_up_memory();

      // Initialise block stuff, identity map between blocks
      block_setup();

      // Get blocks
      Jaa_pt = new CRDoubleMatrix;
      get_block(0, 0, *Jaa_pt);
      CRDoubleMatrix Jab = get_block(0, 1); // Jba == -Jab


      // Setup first block's preconditioner
      // ============================================================
      Prec1_pt->setup(Jaa_pt);


      // Setup second block's preconditioner
      // ============================================================

      // Add matrices, result goes into Jab
      cr_matrix_add(*Jaa_pt, Jab, Jab);

      // Construct preconditioner for application of inverse
      Prec2_pt->setup(&Jab);
    }


    void preconditioner_solve(const DoubleVector &r, DoubleVector &z)
    {
      // Get the right hand side vector in block form.
      Vector<DoubleVector> block_r;
      this->get_block_vectors(r, block_r);

      // Make sure output vector is built
      if (!z.built()) { z.build(this->distribution_pt(),0.0); }

      // Vector of vectors for storage of solution in block form.
      Vector<DoubleVector> block_z(2);

      // Apply first prec
      Prec1_pt->preconditioner_solve(block_r[0], block_z[0]);

      // Apply second prec
      DoubleVector temp1, temp2;
      Prec2_pt->preconditioner_solve(block_r[1], temp1);
      Jaa_pt->multiply(temp1, temp2);
      Prec2_pt->preconditioner_solve(temp2, block_z[1]);

      // copy solution in block vectors block_z back to z
      this->return_block_vectors(block_z,z);
    }

    Preconditioner* Prec1_pt;
    Preconditioner* Prec2_pt;

    CRDoubleMatrix* Jaa_pt;

  private:
    /// Broken copy constructor
    LLGSubBlockPreconditioner(const LLGSubBlockPreconditioner& dummy)
    {BrokenCopy::broken_copy("LLGSubBlockPreconditioner");}

    /// Broken assignment operator
    void operator=(const LLGSubBlockPreconditioner& dummy)
    {BrokenCopy::broken_assign("LLGSubBlockPreconditioner");}

  };




} // End of oomph namespace

#endif
