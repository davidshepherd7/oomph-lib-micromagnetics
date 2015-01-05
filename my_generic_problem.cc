

#include "my_generic_problem.h"
#include "micromagnetics_element.h"

#include "../../src/generic/problem.h"
#include "../../src/generic/sum_of_matrices.h"

#include "sum_of_matrices_preconditioner.h"

#include "magnetics_helpers.h"

#include "../../src/generic/trapezoid_rule.h"


namespace oomph
{

  /// \short Integrate a function given by func_pt over every element
  /// in every bulk mesh in this problem.
  double MyProblem::integrate_over_problem(const ElementalFunction* func_pt,
                                           const Integral* quadrature_pt) const
  {
    using namespace MManipulation;

    double result = 0;
    for(unsigned j=0; j<this->nsub_mesh(); j++)
      {
        if(mesh_pt(j)->finite_element_pt(0)->dim() == this->dim())
          {
            result += integrate_over_mesh(func_pt, mesh_pt(j), quadrature_pt);
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

        // Check if we are using a SumOfMatrices, if so we don't want to
        // get the whole matrix because it will be very slow! Just do
        // outputs with the main matrix instead (I'm assuming it's a CR
        // matrix). ??ds Due to oomph-lib's stupid implementation of
        // templated linear solvers it's impossible to tell which
        // get_jacobian function to call in the general case. So this will
        // not detect other linear solvers that use SumOfMatrices :(
        if(dynamic_cast<GMRES<SumOfMatrices>* >(linear_solver_pt()) != 0)
          {
            this->get_jacobian(residuals, J2);
            J_pt = checked_dynamic_cast<CRDoubleMatrix*>(J2.main_matrix_pt());
          }
        else
          {
            // Need pointer for use in same way as pointer from
            // SumOfMatricse. Store actual data in J so that variable is
            // scoped properly.
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
        mass_matrix_solver_for_explicit_timestepper_pt() = expl_solver_pt;

        // expl_solver_pt->enable_doc_convergence_history();

        // Store + re-use the mass matrix used in explicit steps (since we
        // are almost certainly not going to do spatially adaptivity
        // anytime soon this is safe).
        this->enable_mass_matrix_reuse();
      }


    // If we requested exact solution output then check we have a solution
    if(Want_doc_exact && Exact_solution_pt == 0)
      {
        std::string warning = "Requested doc'ing exact solution, but we don't have an exact solution function pointer.";
        OomphLibWarning(warning, OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);

      }
  }


  /// Hook to be overloaded with any calculations needed after setting of
  /// initial conditions.
  void MyProblem::actions_after_set_initial_condition()
  {
    // If using TR calculate initial derivative with these initial conditions
    TR* tr_pt = dynamic_cast<TR*>(time_stepper_pt());
    if(tr_pt != 0)
      {
        tr_pt->setup_initial_derivative(this);
      }
  }

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

    const std::string& dir = Doc_info.directory();

    // Output Jacobian if requested
    if((to_lower(Doc_info.output_jacobian) == "at_start")
       || (to_lower(Doc_info.output_jacobian) == "always"))
      {
        dump_current_mm_or_jacobian_residuals("at_start");
      }

    // pvd files
    // ============================================================
    // Write start of .pvd XML file
    std::ofstream pvd_file((dir + "/" + "soln.pvd").c_str(),
                           std::ios::out);
    pvd_file << "<?xml version=\"1.0\"?>" << std::endl
             << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">"
             << std::endl
             << "<Collection>" << std::endl;
    pvd_file.close();

    // Write start of exact.pvd XML file
    if(doc_exact())
      {
        std::ofstream pvd_file((dir + "/" + "exact.pvd").c_str(),
                               std::ios::out);
        pvd_file << "<?xml version=\"1.0\"?>" << std::endl
                 << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">"
                 << std::endl
                 << "<Collection>" << std::endl;
        pvd_file.close();
      }

    // Trace file
    // ============================================================

    // Clear (by overwriting) and write headers
    std::ofstream trace_file((dir + "/" + Trace_filename).c_str());
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
      << Trace_seperator << "newton_residuals"
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
    std::ofstream info_file((dir + "/" + Info_filename).c_str());
    info_file
      << "real_time " << real_date_time() << std::endl
      << "unix_time " << std::time(0) << std::endl
      << "driver_name " << problem_name() << std::endl
      << "initial_nnode " << mesh_pt()->nnode() << std::endl
      << "initial_nelement " << mesh_pt()->nelement() << std::endl
      << "initial_nsub_mesh " << nsub_mesh() << std::endl;

    info_file << Doc_info.args_str;
    info_file.close();

    // If requested then output history values before t=0.
    if(Output_initial_condition)
      {
#ifdef PARANOID
        if(ntime_stepper() != 1)
          {
            std::string err = "More/less that 1 time stepper, not sure what to output.";
            throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif

        // Output info for each history value in order.
        const unsigned nval = time_stepper_pt()->nprev_values();
        for(unsigned it=nval-1; it>0; it--)
          {
            doc_solution(it);
          }
      }

    // Output boundary numbers for each boundary node
    // ============================================================
    if(Should_doc_boundaries)
      {
        this->doc_boundaries(dir + "/nodes_on_boundary");
      }


    // Write initial solution and anything else problem specific
    // (e.g. more trace file headers)
    initial_doc_additional();

    // Doc the initial condition
    this->doc_solution();
  }

  void MyProblem::doc_solution(const unsigned& t_hist,
                               const std::string& prefix)
  {
    bool doc_this_step = true;
    if(!is_steady())
      {
        doc_this_step = should_doc_this_step(time_pt()->dt(t_hist),
                                             time_pt()->time(t_hist));
      }

    const std::string dir = Doc_info.directory();
    const std::string num = to_string(Doc_info.number());

    if(Always_write_trace || doc_this_step)
      {
        // Always output trace file data
        write_trace(t_hist);
      }

    // Output full set of data if requested for this timestep
    if(doc_this_step)
      {

        // Solution itself
        std::ofstream soln_file((dir + "/" + prefix + "soln" + num + ".dat").c_str(),
                                std::ios::out);
        soln_file.precision(Output_precision);
        output_solution(t_hist, soln_file);
        soln_file.close();

        // Exact solution if available and requested
        if(doc_exact())
          {
            std::ofstream exact_file((dir + "/" + prefix + "exact" + num + ".dat").c_str(),
                                     std::ios::out);
            exact_file.precision(Output_precision);
            output_exact_solution(t_hist, exact_file);
            exact_file.close();
          }

        // If not a steady state problem then write time information
        if(!is_steady() && prefix == "")
          {
            // Write the simulation time and filename to the solution pvd
            // file
            std::ofstream pvd_file((dir + "/" + "soln.pvd").c_str(),
                                   std::ios::app);
            pvd_file.precision(Output_precision);

            pvd_file << "<DataSet timestep=\"" << time_pt()->time(t_hist)
                     << "\" group=\"\" part=\"0\" file=\"" << "soln"
                     << num << ".vtu"
                     << "\"/>" << std::endl;

            pvd_file.close();


            // Write the simulation time and filename to the exact solution
            // pvd file
            if(doc_exact())
              {
                std::ofstream exact_pvd_file((dir + "/" + "exact.pvd").c_str(),
                                             std::ios::app);
                exact_pvd_file.precision(Output_precision);

                exact_pvd_file << "<DataSet timestep=\"" << time()
                               << "\" group=\"\" part=\"0\" file=\"" <<  "exact"
                         << num << ".vtu"
                         << "\"/>" << std::endl;

                exact_pvd_file.close();
              }
          }


        // Maybe dump the restart data
        if(Dump)
          {
            std::ofstream dump_file((dir + "/" + "dump" + num + ".dat").c_str(),
                                    std::ios::out);
            this->dump(dump_file);
            dump_file.close();
          }

        // Maybe dump ltes
        if(Output_ltes)
          {
            std::ofstream ltefile((dir + "/ltes" + num + ".csv").c_str(),
                                  std::ios::out);
            output_ltes(t_hist, ltefile);
            ltefile.close();
          }

        // Maybe write out predicted values, check we only have one time
        // stepper first though.
#ifdef PARANOID
        if(Output_predictor_values && ntime_stepper() != 1)
          {
            std::string err = "Can only output predictor values for a single time stepper";
            throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);

            // Otherwise we would have multiple "predictor_time" values
            // below.
          }
#endif

        if(t_hist != 0)
          {
            std::string err = "Can't output history of predicted value: they aren't stored.";
            OomphLibWarning(err, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
          }
        else if(Output_predictor_values && time_stepper_pt()->adaptive_flag())
          {
            const unsigned predictor_time =
              time_stepper_pt()->predictor_storage_index();

            std::ofstream pred_file((dir + "/" + "predsoln" + num + ".dat").c_str(),
                                    std::ios::out);
            pred_file.precision(Output_precision);
            output_solution(predictor_time, pred_file, 2);
            pred_file.close();
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
        double min_size = mesh_pt(0)->finite_element_pt(0)->compute_physical_size();

        // Loop over all meshes in problem
        for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
          {
            Mesh* mesh_pt = this->mesh_pt(msh);
            for(unsigned ele=0, nele=mesh_pt->nelement(); ele<nele; ele++)
              {
                FiniteElement* ele_pt = mesh_pt->finite_element_pt(ele);
                double new_size = ele_pt->compute_physical_size();
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

  void MyProblem::write_trace(const unsigned& t_hist)
  {
    std::ofstream trace_file((Doc_info.directory() + "/" + Trace_filename).c_str(),
                             std::ios::app);
    trace_file.precision(Output_precision);

    double time = Dummy_doc_data, dt = Dummy_doc_data;
    if(!is_steady())
      {
        time = this->time_pt()->time(t_hist);
        dt = this->time_pt()->dt(t_hist);
      }

    // Some values are impossible to get for history values, print dummy
    // values instead.
    double nnewton_iter_taken = Dummy_doc_data,
      lte_norm = Dummy_doc_data,
      real_time = Dummy_doc_data,
      total_step_time = Dummy_doc_data;

    Vector<double> solver_iterations(1, Dummy_doc_data),
      solver_times(1, Dummy_doc_data),
      jacobian_setup_times(1, Dummy_doc_data),
      preconditioner_setup_times(1, Dummy_doc_data),
      max_res(1, Dummy_doc_data);

    if(t_hist == 0)
      {
        nnewton_iter_taken = Nnewton_iter_taken;
        solver_iterations = Solver_iterations;
        solver_times = Solver_times;
        jacobian_setup_times = Jacobian_setup_times;
        preconditioner_setup_times = Preconditioner_setup_times;
        real_time = std::time(0);
        total_step_time = Total_step_time;
        max_res = Max_res;

        if(!is_steady())
          {
            lte_norm = this->lte_norm();
          }
      }

    // Write out data that can be done for every problem
    trace_file
      << Doc_info.number()
      << Trace_seperator << time
      << Trace_seperator << dt
      << Trace_seperator << get_error_norm(t_hist)

      << Trace_seperator << nnewton_iter_taken
      << Trace_seperator << solver_iterations

      << Trace_seperator << solver_times
      << Trace_seperator << jacobian_setup_times
      << Trace_seperator << preconditioner_setup_times

      << Trace_seperator << lte_norm
      << Trace_seperator << trace_values(t_hist)

      << Trace_seperator << real_time
      << Trace_seperator << max_res
      << Trace_seperator << get_solution_norm(t_hist)
      << Trace_seperator << total_step_time

      // Reserved slots in case I think of more things to add later
      << Trace_seperator << Dummy_doc_data
      << Trace_seperator << Dummy_doc_data
      << Trace_seperator << Dummy_doc_data
      << Trace_seperator << Dummy_doc_data
      << Trace_seperator << Dummy_doc_data
      << Trace_seperator << Dummy_doc_data;


    // Add problem specific data
    write_additional_trace_data(t_hist, trace_file);

    // Finish off this line
    trace_file << std::endl;
    trace_file.close();
  }

  void MyProblem::doc_boundaries(const std::string& boundary_file_basename) const
  {
    const unsigned n_boundary = mesh_pt()->nboundary();
    for(unsigned b=0; b<n_boundary; b++)
      {
        std::ofstream boundary_file((boundary_file_basename +
                                     to_string(b) + ".csv").c_str(),
                                    std::ios::out);

        // write headers
        boundary_file << "x,y,z,b" << std::endl;

        const unsigned n_boundary_node = mesh_pt()->nboundary_node(b);
        for(unsigned nd=0; nd<n_boundary_node; nd++)
          {
            Node* node_pt = mesh_pt()->boundary_node_pt(b, nd);
            Vector<double> position = node_pt->position();

            // Write out this node
            for(unsigned i=0; i<dim(); i++)
              {
                boundary_file << position[i] << ",";
              }
            for(unsigned i=dim(); i<3; i++)
              {
                boundary_file << 0.0 << ",";
              }
            boundary_file << b << std::endl;
          }

        boundary_file.close();
      }
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

  void MyProblem::set_initial_condition(const InitialConditionFct& ic)
  {
#ifdef PARANOID
    // Can't set global data from a function of space... have to overload this
    // if you have any
    if(nglobal_data() != 0)
      {
        std::string err = "Problem has global data which cannot be set from function pt.";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif

    // Loop over current & previous timesteps (so we need t<nprev_values+1).
    const int nprev_values = this->time_stepper_pt()->nprev_values();
    for(int tindex=0; tindex<nprev_values+1; tindex++)
      {
        double time = time_pt()->time(tindex);

        // Loop over all nodes in all meshes in problem and set values.
        const unsigned nmsh = nsub_mesh();
        for(unsigned msh=0; msh<nmsh; msh++)
          {
            Mesh* mesh_pt = this->mesh_pt(msh);

            for(unsigned nd=0, nnd=mesh_pt->nnode(); nd<nnd; nd++)
              {
                Node* nd_pt = mesh_pt->node_pt(nd);

                // Get the position at present time
                const unsigned dim = nd_pt->ndim();
                Vector<double> x(dim);
                nd_pt->position(x);

                // Set position at tindex time to be the same as at present
                // (impulsive positions).
                for(unsigned j=0; j<dim; j++)
                  {
                    nd_pt->x(tindex, j) = x[j];
                  }

                // Get the values
                Vector<double> values = ic(time, x);

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
                    nd_pt->set_value(tindex, j, values[j]);
                  }
              }

#ifdef PARANOID
            // Can't set internal data like this so check that we have none.
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

}
