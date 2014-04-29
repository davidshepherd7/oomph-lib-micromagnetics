

#include "llg_problem.h"
#include "llg_factories.h"

#include "oomph_factories.h"
#include "micromag_types.h"

#include "micromagnetics_element.h"
#include "magnetostatic_field_flux_element.h"

namespace oomph
{

  /// Function that does the real work of the constructors.
  void LLGProblem::build(Vector<Mesh*>& bulk_mesh_pts)
  {
#ifdef PARANOID
    if(Residual_calculator_pt == 0)
      {
        std::ostringstream error_msg;
        error_msg
          << "Must the following pointers to non-null values before calling build():\n"
          << "Residual_calculator_pt (= " << Residual_calculator_pt << ")"
          << std::endl;

        throw OomphLibError(error_msg.str(), OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

    if(Analytic_ms_fct_pt != 0 && !Disable_ms)
      {
        std::string err = "Other ms must be disabled to use analytical ms.";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif

    // Call the underlying build to deal with adding meshes and time stepper
    MyProblem::build(bulk_mesh_pts);

    // Finish off element build, at this point we should have only micromag
    // elements in the meshes (so we can loop over all meshes) but we don't
    // have a global mesh yet.
    for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
      {
        for(unsigned i=0; i<mesh_pt(msh)->nelement(); i++)
          {
            MicromagEquations* elem_pt =checked_dynamic_cast<MicromagEquations*>
              (mesh_pt(msh)->element_pt(i));

            // Set whether the Jacobian should be finite differenced
            elem_pt->Use_fd_jacobian = Use_fd_jacobian;

            // Set values for magnetic parameters
            elem_pt->magnetic_parameters_pt() = mag_parameters_pt();

            // Set the residual calculation function
            elem_pt->Residual_calculator_pt = Residual_calculator_pt;

            if(Analytic_ms_fct_pt != 0)
              {
                AnalyticalMagnetostatics* ams_pt = new AnalyticalMagnetostatics;
                ams_pt->Magnetostatic_field_fct_pt = Analytic_ms_fct_pt;
                elem_pt->Ms_calc_pt = ams_pt;
              }
          }
      }

    // For debugging we might want to pin m values on the boundary
    if(Pin_boundary_m)
      {
        // Loop over all meshes in problem
        for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
          {
            Mesh* mesh_pt = this->mesh_pt(msh);
            for(unsigned b=0, nb=mesh_pt->nboundary(); b<nb; b++)
              {
                for(unsigned nd=0, nnd=mesh_pt->nboundary_node(b); nd<nnd; nd++)
                  {
                    Node* nd_pt = mesh_pt->boundary_node_pt(b, nd);
                    for(unsigned j=0; j<3; j++)
                      {
                        nd_pt->pin(m_index(j));
                      }
                  }
              }

          }
      }

    // Set up bem stuff if we are doing it
    if(fembem_ms_flag())
      {
        // Figure out how to build the flux meshes that we're going to need
        // for neumann boundaries.
        Flux_mesh_factory_pt = Factories::mm_flux_mesh_factory_factory
          (bulk_mesh_pts[0]->finite_element_pt(0));

        // Loop over all meshes in problem
        for(unsigned msh=0, nmsh=bulk_mesh_pts.size(); msh<nmsh; msh++)
          {
            Vector<unsigned> boundaries;
            for(unsigned b=0; b<bulk_mesh_pts[msh]->nboundary(); b++)
              {boundaries.push_back(b);}

            // Set up neumann condition on phi_1 boundary values (using flux mesh)
            Flux_mesh_pt = flux_mesh_factory(bulk_mesh_pts[msh], boundaries);

            // Add to global mesh
            this->add_sub_mesh(Flux_mesh_pt);

            // Pin a phi_1 value which isn't involved in the boundary element
            // method (we have to pin something to avoid a singular Jacobian,
            // pinning boundary nodes makes BEM much more complex)
            if(pin_a_bulk_phi_1())
              {
#ifdef OOMPH_HAS_MPI
                // In parallel we need to make sure that only one node is
                // pinned in total
                std::string err = "Not implemented!";
                throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                    OOMPH_CURRENT_FUNCTION);
                // Check that processor id is 0, if so then pin as for
                // serial, otherwise do nothing? ??ds Could be problems
                // when nodes duplicated? Not sure how all that works

#else
                Node* pinned_phi_1_node_pt = bulk_mesh_pts[msh]
                 ->get_some_non_boundary_node();
                pinned_phi_1_node_pt->pin(phi_1_index());
                pinned_phi_1_node_pt->set_value(phi_1_index(), 0.0);
#endif
              }

            // Sometimes we don't have any non-boundary nodes (this doesn't
            // work with Hlib yet).
            if(pin_any_phi_1())
              {
#ifdef OOMPH_HAS_MPI
                // In parallel we need to make sure that only one node is
                // pinned in total
                std::string err = "Not implemented!";
                throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                    OOMPH_CURRENT_FUNCTION);
                // Check that processor id is 0, if so then pin as for
                // serial, otherwise do nothing? ??ds Could be problems
                // when nodes duplicated? Not sure how all that works
#else
                // Just grab the first node and pin it
                Node* pinned_phi_1_node_pt = bulk_mesh_pts[msh]->node_pt(0);
                pinned_phi_1_node_pt->pin(phi_1_index());
                pinned_phi_1_node_pt->set_value(phi_1_index(), 0.0);
#endif
              }

          }

      }
    // Otherwise pin all phi and phi_1 dofs to zero
    else if(Disable_ms)
      {
        oomph_info << "Pinning phi values in main problem's meshes." << std::endl;

        // Loop over all meshes in problem
        for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
          {
            for(unsigned nd=0, nnode=mesh_pt(msh)->nnode(); nd<nnode; nd++)
              {
                Node* nd_pt = mesh_pt(msh)->node_pt(nd);

                nd_pt->pin(phi_index());
                nd_pt->pin(phi_1_index());
                nd_pt->set_value(phi_index(),0.0);
                nd_pt->set_value(phi_1_index(),0.0);
              }
          }
      }
    else
      {
        std::string err = "Not sure how to set up ms...";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }



    // Select solver parameters to use for phi solves.
    // ============================================================

    // Start with current ones as defaults.
    get_solver_parameters(Phi_seg_solve_parameters);

    if(!Disable_magnetostatic_solver_optimistations)
      {
        // Optimisations for linear problems
        Phi_seg_solve_parameters.jacobian_reuse_is_enabled = true;
        Phi_seg_solve_parameters.problem_is_nonlinear = false;

        // A good solver
        Phi_seg_solve_parameters.linear_solver_pt
          = Factories::linear_solver_factory("cg");
        checked_dynamic_cast<IterativeLinearSolver*>(Phi_seg_solve_parameters.linear_solver_pt)
          ->preconditioner_pt() = Factories::preconditioner_factory("poisson-amg");
      }
    else
      {
        Phi_seg_solve_parameters.linear_solver_pt
          = Factories::linear_solver_factory("superlu");
      }

    // Similarly for phi1 (keep them separate because stored Jacobians
    // differ).
    get_solver_parameters(Phi_1_seg_solve_parameters);
    if(!Disable_magnetostatic_solver_optimistations)
      {
        Phi_1_seg_solve_parameters.jacobian_reuse_is_enabled = true;
        Phi_1_seg_solve_parameters.problem_is_nonlinear = false;
        Phi_1_seg_solve_parameters.linear_solver_pt
          = Factories::linear_solver_factory("cg");
        checked_dynamic_cast<IterativeLinearSolver*>(Phi_1_seg_solve_parameters.linear_solver_pt)
          ->preconditioner_pt() = Factories::preconditioner_factory("poisson-amg");
      }
    else
      {
        Phi_seg_solve_parameters.linear_solver_pt
          = Factories::linear_solver_factory("superlu");
      }



    // Finish building
    // ============================================================

    // Build the global mesh
    this->build_global_mesh();

    // Number the equations
    this->assign_eqn_numbers();

    // Write out some stuff
    mag_parameters_pt()->output(*oomph_info.stream_pt());
    oomph_info << "LLG Number of equations: " << ndof() << std::endl;
    oomph_info << "Number of sub meshes: " << this->nsub_mesh() << std::endl;

  }

  /// \short Error for adaptive timestepper (rms of nodal error determined by
  /// comparison with explicit timestepper result).
  double LLGProblem::global_temporal_error_norm()
  {
    double global_error = 0.0;

    // Loop over the nodes
    for(unsigned i=0, ni=mesh_pt()->nnode(); i<ni; i++)
      {
        Node* nd_pt = mesh_pt()->node_pt(i);
        for(unsigned j=0; j<3; j++)
          {
            // Get timestepper's error estimate for this direction of m
            // at this point.
            double error = nd_pt->time_stepper_pt()->
              temporal_error_in_value(nd_pt, m_index(j));

            //Add the square of the individual error to the global error
            global_error += error*error;
          }
      }

    // Divide by the number of data points
    global_error /= 3*double(mesh_pt()->nnode());

    return std::sqrt(global_error);
  }


  /// \short Loop over all nodes in bulk mesh and get magnetisations
  Vector<Vector<double> > LLGProblem::get_nodal_magnetisations(unsigned i_time) const
  {
    // Get the space needed
    unsigned nnode = mesh_pt()->nnode();
    Vector< Vector<double> > m_list(nnode, Vector<double>(3, 0.0));

    for(unsigned nd=0, nnd=mesh_pt()->nnode(); nd<nnd; nd++)
      {
        for(unsigned j=0; j<3; j++)
          {
            m_list[nd][j] =
              mesh_pt()->node_pt(nd)->value(i_time, m_index(j));
          }
      }

    return m_list;
  }

  /// \short Solve for the magnetostatic field.
  void LLGProblem::magnetostatics_solve()
  {
    // Do nothing if no solve is needed
    if(!fembem_ms_flag()) return;

    // Check we're not inside a segregated solve already
    check_not_segregated(OOMPH_CURRENT_FUNCTION);

    Inside_segregated_magnetostatics = true;

    // We really need c++11, this array initialisation is ridiculous!
    Vector<unsigned> non_phi_1_indices, non_phi_indices;
    non_phi_1_indices.push_back(phi_index());
    non_phi_1_indices.push_back(m_index(0));
    non_phi_1_indices.push_back(m_index(1));
    non_phi_1_indices.push_back(m_index(2));

    non_phi_indices.push_back(phi_1_index());
    non_phi_indices.push_back(m_index(0));
    non_phi_indices.push_back(m_index(1));
    non_phi_indices.push_back(m_index(2));

    oomph_info << std::endl
               << "Decoupled BEM solve" << std::endl
               << "--------------------------" <<std::endl;


    SolverParameters previous_solver_parameters;
    get_solver_parameters(previous_solver_parameters);


    // solve for phi1
    // ============================================================
    oomph_info << "solving phi1" << std::endl;


    set_solver_parameters(Phi_1_seg_solve_parameters);

    segregated_pin_indices(non_phi_1_indices);
    newton_solve();
    undo_segregated_pinning();

    get_solver_parameters(Phi_1_seg_solve_parameters);



    // pin and set boundary values of phi via bem
    // ============================================================

    oomph_info << "solving BEM" << std::endl;
    double t_start = TimingHelpers::timer();

    // Get bem values. Note that dofs must be in the same equation
    // numbering as when the bem handler was built at this point for this
    // to work (due to how the lookup schemes work). In particular the
    // pinning/segregated pinning MUST be the same. Additionally the phi
    // dofs cannot be pinned (although this can be hacked around by setting
    // up the lookup scheme to use a different index to the real index and
    // setting the pinned values by hand).
    bem_handler_pt()->get_bem_values_and_copy_into_values();

    double t_end = TimingHelpers::timer();
    oomph_info << "BEM time taken: " << t_end - t_start << std::endl;


    // solve for phi
    // ============================================================

    oomph_info << "solving phi" << std::endl;

    // boundary values of phi need to be pinned, use segregated pinning
    // number so that it can be easily undone.
    for(unsigned j=0; j<bem_handler_pt()->Bem_boundaries.size(); j++)
      {
        const Mesh* mesh_pt = bem_handler_pt()->Bem_boundaries[j].second;
        unsigned b = bem_handler_pt()->Bem_boundaries[j].first;

        for(unsigned nd=0, nnd=mesh_pt->nboundary_node(b); nd<nnd; nd++)
          {
            Node* nd_pt = mesh_pt->boundary_node_pt(b, nd);
            if(!nd_pt->is_pinned(phi_index()))
              {
                nd_pt->eqn_number(phi_index())
                  = Data::Is_segregated_solve_pinned;
              }
          }
      }

    set_solver_parameters(Phi_seg_solve_parameters);

    segregated_pin_indices(non_phi_indices);
    newton_solve();
    undo_segregated_pinning();

    get_solver_parameters(Phi_seg_solve_parameters);



    // Done
    // ============================================================

    set_solver_parameters(previous_solver_parameters);

    // oomph_info << "mean field is " << average_magnetostatic_field() << std::endl;

    Inside_segregated_magnetostatics = false;
  }


  /// Linearly extrapolate phi
  void LLGProblem::extrapolate_phi(const double& new_dt, const double& prev_dt)
  {
    // Don't change phi_1 values because we don't need them except for
    // calculating phi.

    double dtn = time_stepper_pt()->time_pt()->dt();
    double dtnm1 = time_stepper_pt()->time_pt()->dt(1);
    const unsigned phi_index = this->phi_index();

    // Loop over all meshes in problem
    for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
      {
        Mesh* mesh_pt = this->mesh_pt(msh);
        for(unsigned nd=0, nnd=mesh_pt->nnode(); nd<nnd; nd++)
          {
            Node* nd_pt = mesh_pt->node_pt(nd);

            double phi_nm1 = nd_pt->value(2, phi_index);
            double phi_n = nd_pt->value(1, phi_index);
            double phi_np1 = ((dtn + dtnm1)/dtnm1)*phi_n - (dtn/dtnm1)*phi_nm1;

            nd_pt->set_value(0, phi_index, phi_np1);
          }
      }

  }


  /// \short Abs of mean difference of actual m and m given by a function
  /// at the middle of each element.
  double LLGProblem::compare_m_with_function(const SolutionFunctor& fct) const
  {
    double diff = 0.0;

    // Compare at middle of element
    Vector<double> s(3,0.0);
    for(unsigned j=0; j<dim(); j++) s[j] = 0.5;

    // Sum the difference over all bulk elements in problem
    unsigned bulk_ele_count = 0;

    // Loop over all meshes in problem
    for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
      {
        // Skip non-bulk meshes
        if((mesh_pt(msh)->nnode() == 0)
           || (mesh_pt(msh)->node_pt(0)->ndim() != Dim)) continue;

        for(unsigned e=0, ne=mesh_pt(msh)->nelement(); e < ne; e++)
          {
            // Get m and x
            MMInterpolator intp(mesh_pt(msh)->finite_element_pt(e), s);
            Vector<double> numerical_m = intp.m();
            Vector<double> x = intp.x();

            Vector<double> exact_m = fct(time(), x);

            for(unsigned j=0; j<3; j++)
              {
                diff += std::abs(numerical_m[j] - exact_m[j]);
              }

            bulk_ele_count++;
          }
      }

    // Divide to get the mean
    diff /= (3.0 * double(bulk_ele_count));

    return diff;
  }


  double LLGProblem::mean_norm_m_error() const
  {
    double temp_error = 0;

    // Loop over all nodes in problem
    for(unsigned nd=0, nnd=mesh_pt()->nnode(); nd<nnd; nd++)
      {
        Node* nd_pt = mesh_pt()->node_pt(nd);
        Vector<double> m_values(3, 0.0);
        for(unsigned j=0; j<3; j++) m_values[j] = nd_pt->value(m_index(j));

        double m_2norm = VectorOps::two_norm(m_values);
        temp_error += std::abs(m_2norm - 1);
      }

    return temp_error/double(mesh_pt()->nnode());
  }

  /// \short Get a vector of the mean of the nodal magnetisations
  Vector<double> LLGProblem::mean_magnetisation() const
  {
    Vector<Vector<double> > ms = get_nodal_magnetisations();
    Vector<double> mean_m(3, 0.0);

    unsigned n_ms = ms.size();
    for(unsigned i=0; i<n_ms; i++)
      {
        mean_m[0] += ms[i][0];
        mean_m[1] += ms[i][1];
        mean_m[2] += ms[i][2];
      }

    mean_m[0] /= n_ms;
    mean_m[1] /= n_ms;
    mean_m[2] /= n_ms;

    return mean_m;
  }


  void LLGProblem::norm_m_error(double &m_error_avg, double &m_error_stddev) const
  {
    // Get mean from other function
    m_error_avg = mean_norm_m_error();

    // Calculate std deviation
    // ============================================================
    double temp_stddev = 0.0;
    for(unsigned nd=0, nnd=mesh_pt()->nnode(); nd<nnd; nd++)
      {
        Node* nd_pt = mesh_pt()->node_pt(nd);
        Vector<double> m_values(3,0.0);
        for(unsigned j=0; j<3; j++) m_values[j] = nd_pt->value(m_index(j));

        double m_2norm = VectorOps::two_norm(m_values);
        temp_stddev += std::pow( m_2norm - 1 - m_error_avg, 2);
      }

    temp_stddev /= double(mesh_pt()->nnode());
    m_error_stddev = std::sqrt(temp_stddev);
  }

}
