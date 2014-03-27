

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

    if(!Disable_mass_matrix_solver_optimisations)
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
        mass_matrix_solver_pt() = expl_solver_pt;

        // expl_solver_pt->enable_doc_convergence_history();

        // Store + re-use the mass matrix used in explicit steps (since we
        // are almost certainly not going to do spatially adaptivity
        // anytime soon this is safe).
        this->enable_mass_matrix_reuse();
      }
  }

  /// Hook to be overloaded with any calculations needed after setting of
  /// initial conditions.
  void MyProblem::actions_after_set_initial_condition() {}

  void MyProblem::initial_doc()
  {
#ifdef PARANOID
    if(*(Doc_info.directory().end()-1) == '/')
      {
        std::string error_msg = "Don't put a / on the end of results dir";
        throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

    // Output Jacobian if requested
    if((to_lower(Doc_info.output_jacobian) == "at_start")
       || (to_lower(Doc_info.output_jacobian) == "always"))
      {
        dump_current_mm_or_jacobian_residuals("at_start");
      }

    // pvd file
    // ============================================================
    // Write start of .pvd XML file
    std::ofstream pvd_file((Doc_info.directory() + "/" + "soln.pvd").c_str(),
                           std::ios::out);
    pvd_file << "<?xml version=\"1.0\"?>" << std::endl
             << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">"
             << std::endl
             << "<Collection>" << std::endl;
    pvd_file.close();


    // Trace file
    // ============================================================

    // Clear (by overwriting) and write headers
    std::ofstream trace_file((Doc_info.directory() + "/" + Trace_filename).c_str());
    trace_file
      << "DocInfo_numbers"
      << Trace_seperator << "times"
      << Trace_seperator << "dts"
      << Trace_seperator << "error_norms"

      << Trace_seperator << "n_newton_iters"
      << Trace_seperator << "n_solver_iters"

      << Trace_seperator << "solver_times"
      << Trace_seperator << "jacobian_setup_times"
      << Trace_seperator << "preconditioner_setup_times"

      << Trace_seperator << "LTE_norms"
      << Trace_seperator << "trace_values"

      << Trace_seperator << "unix_timestamp"
      << Trace_seperator << "dummy"
      << Trace_seperator << "solution_norms"
      << Trace_seperator << "total_step_time"


      // Reserved slots in case I think of more things to add later
      << Trace_seperator << "dummy"
      << Trace_seperator << "dummy"
      << Trace_seperator << "dummy"
      << Trace_seperator << "dummy"
      << Trace_seperator << "dummy"
      << Trace_seperator << "dummy";

    // Other headers depending on the derived class
    write_additional_trace_headers(trace_file);

    // Finish the line and close
    trace_file << std::endl;
    trace_file.close();


    // Info file
    // ============================================================
    std::ofstream info_file((Doc_info.directory() + "/" + Info_filename).c_str());
    info_file
      << "real_time " << real_date_time() << std::endl
      << "unix_time " << std::time(0) << std::endl
      << "driver_name " << problem_name() << std::endl
      << "initial_nnode " << mesh_pt()->nnode() << std::endl
      << "initial_nelement " << mesh_pt()->nelement() << std::endl
      << "initial_nsub_mesh " << nsub_mesh() << std::endl;

    info_file << Doc_info.args_str;
    info_file.close();


    // Write initial solution and anything else problem specific
    // (e.g. more trace file headers)
    this->doc_solution();
    initial_doc_additional();
  }

  void MyProblem::doc_solution()
  {
    bool doc_this_step = true;
    if(!is_steady())
      {
        doc_this_step = should_doc_this_step(time_pt()->dt(), time());
      }

    const std::string dir = Doc_info.directory();
    const std::string num = to_string(Doc_info.number());

    if(Always_write_trace || doc_this_step)
      {
        // Always output trace file data
        write_trace();
      }

    // Output full set of data if requested for this timestep
    if(doc_this_step)
      {

        // Solution itself
        std::ofstream soln_file((dir + "/" + "soln" + num + ".dat").c_str(),
                                std::ios::out);
        soln_file.precision(Output_precision);
        doc_solution_additional(soln_file);
        soln_file.close();

        if(!is_steady())
          {
            // Write the simulation time and filename to the pvd file
            std::ofstream pvd_file((dir + "/" + "soln.pvd").c_str(),
                                   std::ios::app);
            pvd_file.precision(Output_precision);

            pvd_file << "<DataSet timestep=\"" << time()
                     << "\" group=\"\" part=\"0\" file=\"" << "soln"
                     << num << ".vtu"
                     << "\"/>" << std::endl;

            pvd_file.close();
          }


        // Maybe dump the restart data
        if(Dump)
          {
            std::ofstream dump_file((dir + "/" + "dump" + num + ".dat").c_str(),
                                    std::ios::out);
            this->dump(dump_file);
          }


        Doc_info.number()++;
      }
  }

  double MyProblem::min_element_size()
  {
    // Check that we have finite elements ??ds this will still go wrong
    // if there are only some none finite elements in the mesh...

    //??ds what happens with face elements?
    FiniteElement* test_pt = dynamic_cast<FiniteElement*>
      (mesh_pt(0)->element_pt(0));

    if(test_pt != 0)
      {
        double min_size = mesh_pt(0)->finite_element_pt(0)->size();

        // Loop over all meshes in problem
        for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
          {
            Mesh* mesh_pt = this->mesh_pt(msh);
            for(unsigned ele=0, nele=mesh_pt->nelement(); ele<nele; ele++)
              {
                FiniteElement* ele_pt = mesh_pt->finite_element_pt(ele);
                double new_size = ele_pt->size();
                if(new_size < min_size)
                  {
                    min_size = new_size;
                  }
              }
          }

        return min_size;
      }
    // If it's not a finite element then we can't get a size so return
    // a dummy value.
    else
      {
        return Dummy_doc_data;
      }
  }

  void MyProblem::write_trace()
  {
    std::ofstream trace_file((Doc_info.directory() + "/" + Trace_filename).c_str(),
                             std::ios::app);
    trace_file.precision(Output_precision);

    double time = Dummy_doc_data, dt = Dummy_doc_data,
      lte_norm = Dummy_doc_data;
    if(!is_steady())
      {
        time = this->time();
        dt = this->time_pt()->dt();
        lte_norm = this->lte_norm();
      }

    // Write out data that can be done for every problem
    trace_file
      << Doc_info.number()
      << Trace_seperator << time
      << Trace_seperator << dt
      << Trace_seperator << get_error_norm()

      << Trace_seperator << Nnewton_iter_taken
      << Trace_seperator << Solver_iterations

      << Trace_seperator << Solver_times
      << Trace_seperator << Jacobian_setup_times
      << Trace_seperator << Preconditioner_setup_times

      << Trace_seperator << lte_norm
      << Trace_seperator << trace_values()

      << Trace_seperator << std::time(0)
      << Trace_seperator << Dummy_doc_data
      << Trace_seperator << get_solution_norm()
      << Trace_seperator << Total_step_time

      // Reserved slots in case I think of more things to add later
      << Trace_seperator << Dummy_doc_data
      << Trace_seperator << Dummy_doc_data
      << Trace_seperator << Dummy_doc_data
      << Trace_seperator << Dummy_doc_data
      << Trace_seperator << Dummy_doc_data
      << Trace_seperator << Dummy_doc_data;

    //??ds residuals?

    // Add problem specific data
    write_additional_trace_data(trace_file);

    // Finish off this line
    trace_file << std::endl;
    trace_file.close();
  }

  void MyProblem::set_up_impulsive_initial_condition()
  {

#ifdef PARANOID
    if(nglobal_data() != 0)
      {
        std::string err = "Problem has global data which cannot be set from function pt.";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif
    unsigned nprev_steps=this->time_stepper_pt()->nprev_values();
    for(unsigned t=0; t< nprev_steps; t++)
      {
        // Loop over all nodes in all meshes in problem and set values.
        for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
          {
            Mesh* mesh_pt = this->mesh_pt(msh);

            for(unsigned nd=0, nnd=mesh_pt->nnode(); nd<nnd; nd++)
              {
                Node* nd_pt = mesh_pt->node_pt(nd);
                for(unsigned j=0, nj=nd_pt->nvalue(); j<nj; j++)
                  {
                    nd_pt->set_value(t, j, nd_pt->value(0, j));
                  }
              }

#ifdef PARANOID
            for(unsigned ele=0, nele=mesh_pt->nelement(); ele<nele; ele++)
              {
                FiniteElement* ele_pt = mesh_pt->finite_element_pt(ele);
                if(ele_pt->ninternal_data() != 0)
                  {
                    std::string err =
                      "Element with non-nodal data, cannot set via function...";
                    throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                        OOMPH_CURRENT_FUNCTION);
                  }
              }
#endif

          }
      }

    actions_after_set_initial_condition();
  }

  void MyProblem::set_initial_condition(InitialConditionFct& ic_fpt)
  {
#ifdef PARANOID
    if(nglobal_data() != 0)
      {
        std::string err = "Problem has global data which cannot be set from function pt.";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif

    // Loop over current & previous timesteps
    int nprev_steps=this->time_stepper_pt()->nprev_values();
    for(int t=nprev_steps; t>=0; t--)
      {
        double time = time_pt()->time(t);

        // Loop over all nodes in all meshes in problem and set values.
        for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
          {
            Mesh* mesh_pt = this->mesh_pt(msh);

            for(unsigned nd=0, nnd=mesh_pt->nnode(); nd<nnd; nd++)
              {
                Node* nd_pt = mesh_pt->node_pt(nd);

                // Get the position
                const unsigned dim = nd_pt->ndim();
                Vector<double> x(dim);
                nd_pt->position(t, x);

                // Get the values
                Vector<double> values = ic_fpt(time, x);

#ifdef PARANOID
                if(values.size() != nd_pt->nvalue())
                  {
                    std::string err = "Wrong number of values in initial condition.";
                    throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                        OOMPH_CURRENT_FUNCTION);
                  }
#endif
                // Copy into dofs
                for(unsigned j=0, nj=values.size(); j<nj; j++)
                  {
                    nd_pt->set_value(t, j, values[j]);
                  }
              }

#ifdef PARANOID
            for(unsigned ele=0, nele=mesh_pt->nelement(); ele<nele; ele++)
              {
                FiniteElement* ele_pt = mesh_pt->finite_element_pt(ele);
                if(ele_pt->ninternal_data() != 0)
                  {
                    std::string err =
                      "Element with non-nodal data, cannot set via function...";
                    throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                        OOMPH_CURRENT_FUNCTION);
                  }
              }
#endif

            //??ds can't set external/internal data like this though
          }
      }

    actions_after_set_initial_condition();
  }

}
