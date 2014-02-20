

#include "my_generic_problem.h"
#include "micromagnetics_element.h"

#include "../../src/generic/problem.h"
#include "../../src/generic/sum_of_matrices.h"

// just for prec cast, move?
#include "my_general_header.h"


namespace oomph
{

  /// Integrate a function given by func_pt over every element in a mesh
  /// and return the total. This should probably be in the mesh class but
  /// that's core oomph-lib so I'll leave it here.
  double MyProblem::integrate_over_mesh(const ElementalFunction* func_pt,
                                        const Mesh* const mesh_pt) const
  {
    double result = 0;
    for(unsigned e=0, ne=mesh_pt->nelement(); e < ne; e++)
      {
        MicromagEquations* ele_pt
          = checked_dynamic_cast<MicromagEquations*>
          (mesh_pt->element_pt(e));
        result += ele_pt->integrate_over_element(func_pt);
      }
    return result;
  }

  /// \short Integrate a function given by func_pt over every element
  /// in every bulk mesh in this problem.
  double MyProblem::integrate_over_problem(const ElementalFunction* func_pt) const
  {
    double result = 0;
    for(unsigned j=0; j<this->nsub_mesh(); j++)
      {
        if(mesh_pt(j)->finite_element_pt(0)->dim() == this->dim())
          {
            result += integrate_over_mesh(func_pt, mesh_pt(j));
          }
      }
    return result;
  }


  void MyProblem::dump_current_mm_or_jacobian_residuals(const std::string& label)
  {
    // We actually want the mass matrix if we are doing explicit steps
    if(explicit_flag())
      {
        CRDoubleMatrix M;
        DoubleVector residuals;

        // This is how you get mass matrices, ugly...
        AssemblyHandler* old_assembly_handler_pt = this->assembly_handler_pt();
        ExplicitTimeStepHandler ehandler;
        this->assembly_handler_pt() = &ehandler;
        this->get_jacobian(residuals, M);
        this->assembly_handler_pt() = old_assembly_handler_pt;

        M.sparse_indexed_output(Doc_info.directory() + "/massmatrix_" + label,
                                Output_precision, true);
        residuals.output(Doc_info.directory() + "/explicit_residual_" + label,
                         Output_precision);
      }
    else
      {
        CRDoubleMatrix J, *J_pt;
        SumOfMatrices J2;
        DoubleVector residuals;

        // Check if we are using som-gmres, if so we don't want to get
        // the whole matrix because it will be very slow! Just do outputs
        // with the main matrix instead (I'm assuming it's a CR matrix).
        if(dynamic_cast<GMRES<SumOfMatrices>* >(linear_solver_pt()) != 0)
          {
            this->get_jacobian(residuals, J2);
            J_pt = checked_dynamic_cast<CRDoubleMatrix*>(J2.main_matrix_pt());
          }
        else
          {
            this->get_jacobian(residuals, J);
            J_pt = &J;
          }

        J_pt->sparse_indexed_output(Doc_info.directory() + "/jacobian_" + label,
                                    Output_precision, true);
        residuals.output(Doc_info.directory() + "/residual_" + label,
                         Output_precision);

        // Also dump blocks if we have a block preconditioner
        IterativeLinearSolver* its_pt = iterative_linear_solver_pt();
        if(its_pt != 0)
          {
            // Try to get a block preconditioner from the preconditioner
            BlockPreconditioner<CRDoubleMatrix>* bp_pt
              = smart_cast_preconditioner<BlockPreconditioner<CRDoubleMatrix>*>
              (its_pt->preconditioner_pt());

            if(bp_pt != 0)
              {
                // Set up blocks
                bp_pt->set_matrix_pt(J_pt);
                bp_pt->set_comm_pt(J_pt->distribution_pt()->communicator_pt());
                bp_pt->block_setup();

                // Dump the blocks
                std::string basefname = Doc_info.directory() + "/J_" + label;
                bp_pt->output_blocks_to_files(basefname, Output_precision);
              }
          }
      }
  }

  /// \short Perform set up of problem.
  void MyProblem::build(Vector<Mesh*>& bulk_mesh_pts)
  {
    // Copy the first mesh's first timestepper to the problem

    FiniteElement* fele_pt = dynamic_cast<FiniteElement*>
      (bulk_mesh_pts[0]->element_pt(0));

    // Finite element mesh: grab ts from node
    if(fele_pt != 0)
      {
        TimeStepper* ts_pt = bulk_mesh_pts[0]->node_pt(0)->time_stepper_pt();
        this->add_time_stepper_pt(ts_pt);

        // ??ds assumed any timesteppers hiding somewhere else are added elsewhere

#ifdef PARANOID
        for(unsigned j=0; j<bulk_mesh_pts.size(); j++)
          {
            if(bulk_mesh_pts[j]->node_pt(0)->time_stepper_pt()
               != ts_pt)
              {
                std::string err = "Multiple timesteppers, you need to do somedhing more fancy here";
                throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                    OOMPH_CURRENT_FUNCTION);
              }
          }
#endif
      }

    // Non finite element mesh: grab ts from internal data
    else
      {
        TimeStepper* ts_pt = bulk_mesh_pts[0]->element_pt(0)->
          internal_data_pt(0)->time_stepper_pt();
        this->add_time_stepper_pt(ts_pt);

        // ??ds again assumed any timesteppers hiding somewhere else are added elsewhere

#ifdef PARANOID
        for(unsigned j=0; j<bulk_mesh_pts.size(); j++)
          {
            if(bulk_mesh_pts[j]->element_pt(0)->
               internal_data_pt(0)->time_stepper_pt() != ts_pt)
              {
                std::string err = "Multiple timesteppers? you need to do somedhing more fancy here";
                throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                    OOMPH_CURRENT_FUNCTION);
              }
          }
#endif
      }


    // Push all the meshes into the problem's sub mesh list
    for(unsigned j=0; j<bulk_mesh_pts.size(); j++)
      {
        add_sub_mesh(bulk_mesh_pts[j]);
      }

    // If we have an iterative solver with a block preconditioner then
    // add all the meshes to the block preconditioner as well.
    IterativeLinearSolver* its_pt = iterative_linear_solver_pt();
    if(its_pt != 0)
      {
        // Try to get a block preconditioner from the preconditioner
        BlockPreconditioner<CRDoubleMatrix>* bp_pt
          = smart_cast_preconditioner<BlockPreconditioner<CRDoubleMatrix>*>
          (its_pt->preconditioner_pt());

        if(bp_pt != 0)
          {
            // Set up meshes
            bp_pt->set_nmesh(nsub_mesh());
            for(unsigned i=0; i< nsub_mesh(); i++)
              {
                bp_pt->set_mesh(i, mesh_pt(i));
              }
          }
      }

    // Get the problem dimension
    if(fele_pt != 0)
      {
        Dim = fele_pt->nodal_dimension();
      }
    else
      {
        // Presumably if a "bulk" mesh contains non-finite elements
        // then this is not a pde, so no dimension as such.
        Dim = 0;
      }

    if(!Disable_explicit_solver_optimisations)
      {
        // Set the solver for explicit timesteps (mass matrix) to CG with a
        // diagonal predconditioner.
        IterativeLinearSolver* expl_solver_pt = new CG<CRDoubleMatrix>;
        expl_solver_pt->preconditioner_pt() =
          new MatrixBasedLumpedPreconditioner<CRDoubleMatrix>;

        // If it takes more than 100 iterations then something has almost
        // certainly gone wrong!
        expl_solver_pt->max_iter() = 100;
        expl_solver_pt->enable_error_after_max_iter();
        explicit_solver_pt() = expl_solver_pt;

        // expl_solver_pt->enable_doc_convergence_history();

        // Store + re-use the mass matrix used in explicit steps (since we
        // are almost certainly not going to do spatially adaptivity
        // anytime soon this is safe).
        this->enable_mass_matrix_reuse();
      }
  }

}
