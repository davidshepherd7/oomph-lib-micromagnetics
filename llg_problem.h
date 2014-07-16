#ifndef OOMPH_IMPLICIT_LLG_EXCHANGE_PROBLEM_H
#define OOMPH_IMPLICIT_LLG_EXCHANGE_PROBLEM_H

// General includes

// oomph-lib includes

// Micromagnetics includes
#include "micromagnetics_element.h"
#include "vector_helpers.h"
#include "magnetics_helpers.h"
#include "my_generic_problem.h"
#include "mallinson_solution.h"
#include "residual_calculator.h"
#include "boundary_element_handler.h"
#include "micromagnetics_flux_element.h"
#include "generic_poisson_problem.h"
#include "micromag_types.h"

namespace oomph
{

  class ResidualCalculator;

  /// Enumeration for how to handle phi_1's singularity
  namespace phi_1_singularity_handling
  {
    enum phi_1_singularity_handling {pin_any, pin_bulk, pin_boundary,
                                     normalise,
                                     nothing,
                                     jacobian};
  }

  // ============================================================
  ///
  // ============================================================
  class LLGProblem : public MyProblem
  {
  public:

    /// Default constructor - do nothing except nulling pointers.
    LLGProblem() :
      Compare_with_mallinson(false),
      Previous_energies(5, 0.0)
    {
      Boundary_solution_pt = 0;

      // Needed for if we want to switch solvers in runs
      Super_LU_solver_pt = new SuperLUSolver;

      // Storage for residual calculator
      Residual_calculator_pt = 0;

      // Initialise storage for energies etc.
      Exchange_energy = MyProblem::Dummy_doc_data;
      Zeeman_energy = MyProblem::Dummy_doc_data;
      Crystalline_anisotropy_energy = MyProblem::Dummy_doc_data;
      Magnetostatic_energy = MyProblem::Dummy_doc_data;

      Abs_damping_error = MyProblem::Dummy_doc_data;
      Rel_damping_error = MyProblem::Dummy_doc_data;


      // Bem stuff
      Bem_handler_pt = 0;
      Flux_mesh_pt = 0;
      Flux_mesh_factory_pt = 0;
      Phi_1_singularity_method = phi_1_singularity_handling::pin_any;

      Relax_magnetisation = false;

      Decoupled_ms = false;
      Extrapolate_decoupled_ms = false;
      Disable_ms = false;
      Analytic_ms_fct_pt = 0;
      Inside_explicit_timestep = false;
      Inside_segregated_magnetostatics = false;
#ifdef PARANOID
      Check_angles = true;
#else
      Check_angles = false;
#endif
      Disable_magnetostatic_solver_optimistations = false;

      Force_gaussian_quadrature_in_energy = false;


      // Debugging switches
      Pin_boundary_m = false;
      Use_fd_jacobian = false;
      Renormalise_each_time_step = false;
    }

    /// Virtual destructor. Policy decision: my problem classes won't call
    /// delete on anything. That's up to the driver code or
    /// whatever. Ideally we should just start using c++11 smart pointers
    /// and not have to worry so much about memory management!
    virtual ~LLGProblem() {}

    std::string problem_name() const override {return "LLG";}

    /// Get the jacobian as a SumOfMatrices. This is probably the best way
    /// to deal with Jacobians involving a dense block (i.e. in fully
    /// implicit bem). If we aren't using that then this is basically the
    /// same as using the cr matrix form but with the Jacobian wrapped
    /// inside a SumOfMatrices class.
    void get_jacobian(DoubleVector &residuals, SumOfMatrices &jacobian) override
    {
      // Get the fem Jacobian and the residuals
      CRDoubleMatrix* fem_jacobian_pt =
        new CRDoubleMatrix(residuals.distribution_pt());
      Problem::get_jacobian(residuals, *fem_jacobian_pt);

      // maybe modify a phi_1 jacobian entry to remove singularity
      maybe_twiddle_phi_1_jacobian_entry(*fem_jacobian_pt);

      // Assign the fem jacobian to be the "main" sumofmatrices matrix. Do
      // delete it when done.
      jacobian.main_matrix_pt() = fem_jacobian_pt;
      jacobian.set_delete_main_matrix();

      // If we're doing bem here then add on the bem matrix and the
      // d(phibound)/d(phibound) block identity matrix.
      if(implicit_ms_flag())
        {
          // Add the bem matrix to the jacobian in the right places. Don't
          // delete it when done.
          jacobian.add_matrix(bem_handler_pt()->bem_matrix_pt(),
                              bem_handler_pt()->row_lookup_pt(),
                              bem_handler_pt()->col_lookup_pt(),
                              false);

          // Create identity CRDoubleMatrix
          unsigned bem_nnode = bem_handler_pt()->bem_mesh_pt()->nnode();
          unsigned nrow = fem_jacobian_pt->nrow();
          LinearAlgebraDistribution* dist_pt = fem_jacobian_pt->distribution_pt();
          Vector<double> values(bem_nnode, -1.0);
          Vector<int> row_index(bem_nnode), row_start;
          for(unsigned nd=0; nd<bem_nnode; nd++)
            {
              unsigned i = bem_handler_pt()->output_lookup_pt()
                ->added_to_main(nd);
              row_index[nd] = i;
            }

          // Sort the row index list then create row starts vector. DO NOT
          // COPY THIS CODE FOR CREATION OF OTHER MATRICES: sorting is only
          // safe because all the values are the same.
          std::sort(row_index.begin(), row_index.end());
          VectorOps::rowindex2rowstart(row_index, nrow, row_start);

          // Create the matrix (col_index == row index since on diagonal)
          CRDoubleMatrix bem_block_identity(dist_pt, nrow, values,
                                            row_index, row_start);


          // Add it on
          VectorOps::cr_matrix_add(*fem_jacobian_pt, bem_block_identity,
                                   *fem_jacobian_pt);
        }

    }

    void maybe_twiddle_phi_1_jacobian_entry(CRDoubleMatrix &jacobian) const
    {
      if(twiddle_phi_1_jacobian())
        {
          // Find a phi 1 equation
          int phi_1_dof = mesh_pt()->node_pt(0)->eqn_number(phi_1_index());

          // If not pinnned (phi1 is either all pinned or not pinned at all,
          // which is why we have this problem in the first place).
          if(phi_1_dof >= 0)
            {
              // Then add some small amount to the Jacobian to make it
              // non-singular, writing to CR matrices is hard...

              bool success = false;
              for(long k=jacobian.row_start()[phi_1_dof];
                  k<jacobian.row_start()[phi_1_dof+1];
                  k++)
                {
                  if(jacobian.column_index()[k] == phi_1_dof)
                    {
                      jacobian.value()[k] += 1e-3;
                      success = true;
                      break;
                    }
                }

              if(!success)
                {
                  std::string err = "Failed to modify phi1 jacobian entry.";
                  throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                                      OOMPH_EXCEPTION_LOCATION);
                }
            }
        }
    }

    void output_solution(const unsigned& t, std::ostream& outstream,
                         const unsigned& npoints=2) const override
    {
      MyProblem::output_solution(t, outstream, npoints);

      // Output error in |m| if requested.
      // ??ds refactor to merge with other addtional doc stuff somehow?
      if(Doc_m_length_error)
        {
          const std::string dir = Doc_info.directory();
          const std::string num = to_string(Doc_info.number());

          std::ofstream ml_file((dir + "/ml_error" + num + ".csv").c_str(),
                                std::ios::out);
          doc_m_length_error(t, ml_file);
          ml_file.close();
        }
    }

    /// Get the Jacobian as a CRDoubleMatrix (the normal matrix format). If
    /// we are using fully implicit bem then this is not a good idea for
    /// "real" problems but useful for tests. Otherwise this is exactly the
    /// same as Problem::get_jacobian().
    void get_jacobian(DoubleVector &residuals, CRDoubleMatrix &jacobian) override
    {
      // If we're calculating ms here then include the bem matrix. Warning:
      // this is not going to be fast! Use the sumofmatrices version
      // instead if possible.
      if(implicit_ms_flag())
        {
          Vector<double> sum_values;
          Vector<int> sum_rows, sum_cols, sum_row_start;
          unsigned ncol = 0;
          LinearAlgebraDistribution dist;

          // These braces make sure that the fem jacobian is destroyed asap
          // to reduce memory usage.
          {
            // Get as a sum of matrices jacobian
            SumOfMatrices sum_jacobian;
            get_jacobian(residuals, sum_jacobian);

            // Copy out the data we need
            ncol = sum_jacobian.ncol();
            dist = *(checked_dynamic_cast<CRDoubleMatrix*>
                     (sum_jacobian.main_matrix_pt())->distribution_pt());

            // // Convert to indicies ??ds SLOW (N^2)
            // sum_jacobian.get_as_indices(sum_cols, sum_rows, sum_values);
            VectorOps::get_as_indicies(sum_jacobian, sum_values, sum_cols,
                                       sum_rows);

          } // sum_jacobian destroyed -> fem_jacobian destroyed, but
            // information we need is still in the vectors.

          // Convert to rowstart form and make a cr matrix out of it
          VectorOps::rowindex2rowstart(sum_rows, ncol, sum_row_start);
          jacobian.build(&dist, ncol, sum_values, sum_cols, sum_row_start);
        }

      // Otherwise just do it as normal
      else
        {
          Problem::get_jacobian(residuals, jacobian);
        }

      // maybe modify a phi_1 jacobian entry to remove singularity
      maybe_twiddle_phi_1_jacobian_entry(jacobian);
    }

    /// Overload to set flag to avoid BEM interfering with get mass matrix
    /// and to swap to ll form of residual if needed.
    void get_dvaluesdt(DoubleVector& f)
      {
        Inside_explicit_timestep = true;

        bool needs_reset = false;
        if(Residual_calculator_pt->use_gilbert_form())
          {
            Residual_calculator_pt->set_use_ll_form();
            needs_reset = true;
          }

        Problem::get_dvaluesdt(f);

        if(needs_reset)
          {
            Residual_calculator_pt->set_use_gilbert_form();
            Inside_explicit_timestep = false;
          }

        Inside_explicit_timestep = false;
      }

    /// Function that does the real work of the constructors.
    void build(Vector<Mesh*>& bulk_mesh_pts) override;

    /// Renormalise magnetisation to 1 (needed with BDF2)
    void renormalise_magnetisation()
    {
      oomph_info << "Renormalising nodal magnetisations." << std::endl;

      // Loop over meshes and renormalise m at each node
      for(unsigned nd=0; nd<mesh_pt()->nnode(); nd++)
        {
          Node* nd_pt = mesh_pt()->node_pt(nd);

          // Get m vector
          Vector<double> m_values(3,0.0);
          for(unsigned j=0; j<3; j++) m_values[j] = nd_pt->value(m_index(j));

          // Normalise
          VectorOps::normalise(m_values);

          // Write m vector
          for(unsigned j=0; j<3; j++) nd_pt->set_value(m_index(j),m_values[j]);
        }
    }

    /// Relax magnetisation using strongly damped llg until torque is small
    void relax(HAppFctPt h_app_relax=HApp::zero)
    {
      double initial_damping = mag_parameters_pt()->Gilbert_damping;
      double initial_dt = time_pt()->dt();
      double initial_time = time_pt()->time();
      HAppFctPt initial_happ = mag_parameters_pt()->Applied_field_fct_pt;

      // Set a high damping and relaxation field
      Magnetic_parameters_pt->Gilbert_damping = 1.0;
      Magnetic_parameters_pt->Applied_field_fct_pt = h_app_relax;


      // Initialise loop variables
      double dt = std::max(initial_dt/100, 1e-4);
      double relax_tol = 1e-4;

      // ??ds big hack here: just relax for 250 time units, enough for
      // mumag4... fix to use torque eventually
      while(time() < 250)
        {
          // Output some basic info
          oomph_info
            << std::endl
            << std::endl
            << "Relaxation step " << N_steps_taken
            << ", time = " << time()
            << ", dt = " << dt << std::endl
            << "=============================================" << std::endl
            << std::endl;

          // Do the newton solve (different ones depending flags set)
          dt = smart_time_step(dt, relax_tol);

          // Output
          doc_solution(); //??ds extra label?
        }


      // Revert time info
      N_steps_taken = 0;
      time_pt()->time() = initial_time;
      initialise_dt(initial_dt);

      // Set initial condition to be the relaxed values all the way back in
      // time.
      set_up_impulsive_initial_condition();

      // Revert the damping and field
      Magnetic_parameters_pt->Gilbert_damping = initial_damping;
      Magnetic_parameters_pt->Applied_field_fct_pt = initial_happ;
    }

    virtual void actions_before_time_integration() override
    {
      MyProblem::actions_before_time_integration();

      if(Relax_magnetisation)
        {
          relax();
        }
    }


    virtual void actions_before_implicit_timestep() override
    {
      MyProblem::actions_before_implicit_timestep();

      // Project phi to correct time if necessary
      if(Decoupled_ms)
        {
          if(Extrapolate_decoupled_ms)
            {
              extrapolate_phi(time_stepper_pt()->time_pt()->dt(),
                              time_stepper_pt()->time_pt()->dt(1));
            }
        }

      // If we are using Dirichlet boundary conditions then update them.
      if(Pin_boundary_m)
        {
          update_dirichlet_conditions();
        }
    }

    void update_dirichlet_conditions()
    {
      const double t = time_pt()->time();

      const unsigned n_msh = nsub_mesh();
      for(unsigned msh=0; msh<n_msh; msh++)
        {
          Mesh* mesh_pt = this->mesh_pt(msh);
          const unsigned n_boundary = mesh_pt->nboundary();
          for(unsigned b=0; b<n_boundary; b++)
            {
              const unsigned n_boundary_node = mesh_pt->nboundary_node(b);
              for(unsigned nd=0; nd<n_boundary_node; nd++)
                {
                  Node* nd_pt = mesh_pt->boundary_node_pt(b, nd);

                  // Get position
                  Vector<double> x(dim(), 0.0); nd_pt->position(x);

                  // Calculate m
                  Vector<double> m = boundary_solution(t, x);

                  // Write to node data
                  for(unsigned j=0; j<3; j++)
                    {
                      unsigned i = m_index(j);
#ifdef PARANOID
                      if(!nd_pt->is_pinned(i))
                        {
                          std::string err
                            = "Boundary m value not pinned but tried to set condition";
                          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                              OOMPH_CURRENT_FUNCTION);
                        }
#endif

                      nd_pt->set_value(i, m[i]);
                    }
                }
            }
        }

    }

    Vector<double> boundary_solution(const double& t, Vector<double>& x) const
    {
#ifdef PARANOID
      if(Boundary_solution_pt == 0)
        {
          std::string err = "Boundary_solution_pt is null!";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif
      return (*Boundary_solution_pt)(t, x);
    }


    virtual void actions_after_implicit_timestep() override
      {
        MyProblem::actions_after_implicit_timestep();

        if(Decoupled_ms)
          {
            // Unpin phis (from m solve)
            undo_segregated_pinning();

            // Solve for the magnetostatic field.
            magnetostatics_solve();
          }

        // check neighbouring magnetisation angles if requested
        maybe_check_angles();

        // Calculate and store the energy (ready to be output)
        calculate_energies();
      }


    virtual void actions_before_newton_step() override
    {
      oomph_info << std::endl
                 << "Newton step " << Nnewton_iter_taken + 1 << std::endl
                 << "---------------------------------------" << std::endl;

      if(!Inside_segregated_magnetostatics)
        {
          // Call base class version
          MyProblem::actions_before_newton_step();
        }
    }

    virtual void actions_after_newton_step() override
    {
      if(!Inside_segregated_magnetostatics)
        {
          // Call base class version
          MyProblem::actions_after_newton_step();
        }
    }


    virtual void actions_before_newton_solve() override
    {
      if(!Inside_segregated_magnetostatics)
        {
          // Call base class version
          MyProblem::actions_before_newton_solve();

          // Update BEM magnetostatics boundary conditions (if we are doing
          // them fully implicitly).
          if(implicit_ms_flag())
            {
              bem_handler_pt()->get_bem_values_and_copy_into_values();
            }

          // For decoupled fem/bem solves pin everything but m
          else if(fembem_ms_flag())
            {
              check_not_segregated(OOMPH_CURRENT_FUNCTION);

              Vector<unsigned> non_m_indices;
              non_m_indices.push_back(phi_1_index());
              non_m_indices.push_back(phi_index());

              segregated_pin_indices(non_m_indices);
            }
        }
    }


    virtual void actions_before_explicit_timestep() override
    {
      // Set this variable to avoid getting BEM in mass matrix (due to a
      // hack in oomph core.. fix that instead?)
      Inside_explicit_timestep = true;

      // Check we're not inside a segregated solve already
      check_not_segregated(OOMPH_CURRENT_FUNCTION);

      MyProblem::actions_before_explicit_timestep();
    }

    virtual void actions_after_explicit_timestep() override
    {
      MyProblem::actions_after_explicit_timestep();

      // We need to keep M normalised...
      if(renormalise_each_time_step())
        {
          renormalise_magnetisation();
        }

      // check neighbouring magnetisation angles if requested
      maybe_check_angles();

      // Explicit timestep is over now
      Inside_explicit_timestep = false;
    }

    virtual void actions_before_explicit_stage() override
    {
      MyProblem::actions_before_explicit_stage();

      // Check we're not inside a segregated solve already
      check_not_segregated(OOMPH_CURRENT_FUNCTION);

      // Segregated-pin phi values ready for m step
      Vector<unsigned> non_m_indices;
      non_m_indices.push_back(phi_1_index());
      non_m_indices.push_back(phi_index());
      segregated_pin_indices(non_m_indices);
    }

    virtual void actions_after_explicit_stage() override
    {
      // Remove m-only segregation
      undo_segregated_pinning();

      // Solve for the new magnetostatic field.
      magnetostatics_solve();

      MyProblem::actions_after_explicit_stage();
    }

    virtual void actions_after_newton_solve() override
    {
      oomph_info << std::endl
                 << "Finalising" << std::endl
                 << "-------------------------------------------" << std::endl;

      if(!Inside_segregated_magnetostatics)
        {
          MyProblem::actions_after_newton_solve();

          // If we're using BDF we need to keep M normalised.
          if(renormalise_each_time_step())
            {
              renormalise_magnetisation();
            }
        }
    }

    void maybe_check_angles()
      {
        if(Check_angles)
          {
            // From the nmag user manual:
            // [Begin quote M Donahue]
            // * if the spin angle is approaching 180 degrees,
            //   then the results are completely bogus.
            // * over 90 degrees the results are highly questionable.
            // * Under 30 degrees the results are probably reliable.
            // [end quote]
            // (the spin angle is the angle between two neighbouring magnetisations).

            // Check this:
            Vector<double> a = elemental_max_m_angle_variations();
            double max_angle_var = *std::max_element(a.begin(), a.end());
            if(max_angle_var > MathematicalConstants::Pi/6)
              {
                std::string error_msg
                  = "Large angle variations of " + to_string(max_angle_var)
                  + " > " + to_string(MathematicalConstants::Pi/6)
                  + " across a single element,\n";
                error_msg += "this often means that your mesh is not sufficiently refined.";
                throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
          }
      }

    virtual void actions_before_newton_convergence_check() override
      {
        // Call base class actions function
        MyProblem::actions_before_newton_convergence_check();

        // Update BEM magnetostatics boundary conditions (if we are doing them
        // fully implicitly).
        if(implicit_ms_flag())
          {
            bem_handler_pt()->get_bem_values_and_copy_into_values();
          }

        // Maybe subtract largest phi value from all phi values to get it
        // O(1) (technically it has an arbitrary added constant).
        if(normalise_phi_1())
          {
            const unsigned n_node = mesh_pt()->nnode();
            const unsigned index = phi_1_index();

            // Find the min phi1 value
            double min_phi_1 = mesh_pt()->node_pt(0)->value(index);
            for(unsigned nd=1; nd<n_node; nd++)
              {
                double this_phi_1 = mesh_pt()->node_pt(nd)
                  ->value(index);
                if(this_phi_1 < min_phi_1)
                  {
                    min_phi_1 = this_phi_1;
                  }
              }

            std::cout << "first min phi " << min_phi_1 << std::endl;

            // subtract it from each
            for(unsigned nd=0; nd<n_node; nd++)
              {
                double this_phi_1 = mesh_pt()->node_pt(nd)->value(index);
                mesh_pt()->node_pt(0)->set_value(index,
                                                 this_phi_1 - min_phi_1);
              }

            // Find the max phi1 value
            double new_max_phi_1 = mesh_pt()->node_pt(0)->value(index);
            for(unsigned nd=1; nd<n_node; nd++)
              {
                double this_phi_1 = mesh_pt()->node_pt(nd)
                  ->value(index);
                if(this_phi_1 > new_max_phi_1)
                  {
                    new_max_phi_1 = this_phi_1;
                  }
              }
            std::cout << new_max_phi_1 << std::endl;
          }
      }


    void write_additional_trace_headers(std::ofstream& trace_file) const override
    {
      trace_file
        << Trace_seperator << "m_length_error_means"
        << Trace_seperator << "m_length_error_std_devs"
        << Trace_seperator << "max_angle_errors"
        << Trace_seperator << "mean_mxs"
        << Trace_seperator << "mean_mys"
        << Trace_seperator << "mean_mzs"
        << Trace_seperator << "exchange_energy"
        << Trace_seperator << "zeeman_energy"
        << Trace_seperator << "crystalline_anisotropy_energy"
        << Trace_seperator << "magnetostatic_energy"
        << Trace_seperator << "total_energy"
        << Trace_seperator << "abs_damping_error"
        << Trace_seperator << "rel_damping_error"
        << Trace_seperator << "h_applied_first_element"
        << Trace_seperator << "m_length_error_maxes"
        ;

    }

    void write_additional_trace_data(const unsigned& t_hist,
                                     std::ofstream& trace_file) const override
    {
      Vector<Vector<double> > ms = MManipulation::nodal_magnetisations(t_hist,
                                                                       *this);
      Vector<double> ml_errors = MManipulation::nodal_m_length_errors(ms);
      Vector<double> angle_variations = elemental_max_m_angle_variations(t_hist);
      Vector<double> mean_m = MManipulation::mean_nodal_magnetisation(ms);
      Vector<double> h_app = first_element_happ(t_hist);


      // Energies only available for present time
      double exchange_energy = Dummy_doc_data, zeeman_energy = Dummy_doc_data,
        crystalline_anisotropy_energy = Dummy_doc_data,
        magnetostatic_energy = Dummy_doc_data,
        micromagnetic_energy = Dummy_doc_data,
        abs_damping_error = Dummy_doc_data, rel_damping_error = Dummy_doc_data;
      if(t_hist == 0)
        {
          exchange_energy = Exchange_energy;
          zeeman_energy = Zeeman_energy;
          crystalline_anisotropy_energy = Crystalline_anisotropy_energy;
          magnetostatic_energy = Magnetostatic_energy;
          micromagnetic_energy = this->micromagnetic_energy();
          abs_damping_error = Abs_damping_error;
          rel_damping_error = Rel_damping_error;
        }

      trace_file
        << Trace_seperator << VectorOps::mean(ml_errors)
        << Trace_seperator << VectorOps::stddev(ml_errors)
        << Trace_seperator
        << *std::max_element(angle_variations.begin(), angle_variations.end())
        << Trace_seperator << mean_m[0]
        << Trace_seperator << mean_m[1]
        << Trace_seperator << mean_m[2]

        << Trace_seperator << exchange_energy
        << Trace_seperator << zeeman_energy
        << Trace_seperator << crystalline_anisotropy_energy
        << Trace_seperator << magnetostatic_energy
        << Trace_seperator << micromagnetic_energy
        << Trace_seperator << abs_damping_error
        << Trace_seperator << rel_damping_error
        << Trace_seperator << h_app
        << Trace_seperator << VectorOps::max(ml_errors)
        ;
    }

    /// Get happ value from first element
    Vector<double> first_element_happ(const unsigned& t_hist=0) const
    {
      double t = time_pt()->time(t_hist);
      Vector<double> dummy_x, dummy_s;
      Vector<double> H_app = ele_pt()->get_applied_field(t, dummy_x);

      return H_app;
    }

    void initial_doc_additional() const override
    {
      // If we have a H matrix then write out its rank data
      if(Bem_handler_pt != 0)
        {
          bem_handler_pt()->maybe_write_h_matrix_data(Doc_info.directory());
        }
    }

    /// \short Return a vector of the maximum angle variation in each
    /// element.
    Vector<double> elemental_max_m_angle_variations(const unsigned& t_hist=0) const
    {
      Vector<double> angles;

      for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
        {
          // Skip non-bulk meshes
          if((mesh_pt(msh)->nnode() == 0)
             || (mesh_pt(msh)->node_pt(0)->ndim() != Dim)) continue;

          for(unsigned e=0, ne=mesh_pt(msh)->nelement(); e < ne; e++)
            {
              MicromagEquations* ele_pt =
                checked_dynamic_cast<MicromagEquations*>
                (mesh_pt(msh)->element_pt(e));

              angles.push_back(ele_pt->max_m_angle_variation(t_hist));
            }
        }
      return angles;
    }

    /// \short Error for adaptive timestepper (rms of nodal error determined by
    /// comparison with explicit timestepper result).
    double global_temporal_error_norm() override;


    virtual void actions_after_set_initial_condition() override
      {
        MyProblem::actions_after_set_initial_condition();

        // Solve for initial field and phi values
        magnetostatics_solve();

        calculate_energies(false);
      }

    virtual void dump(std::ofstream& dump_file) const override
    {
      // Set very high precision to avoid any issues
      dump_file.precision(14);

      // Let base class handle the rest
      MyProblem::dump(dump_file);
    }

    virtual void read(std::ifstream& restart_file) override
    {
      // Let base class handle the rest
      MyProblem::read(restart_file);

      actions_after_set_initial_condition();
    }



    // Lots of magnetisation manipulation functions
    // ============================================================

    /// \short Abs of mean difference of actual m and m given by a function
    /// at the middle of each element.
    double compare_m_with_function(const SolutionFunctorBase& fct) const;


    double get_error_norm(const unsigned& t_hist) const override
    {
      if(Compare_with_mallinson)
        {
          using namespace CompareSolutions;

          double time = time_pt()->time(t_hist);
          Vector<double> m_now = MManipulation::mean_nodal_magnetisation(t_hist, *this);
          double exact_time = switching_time_wrapper(mag_parameters_pt(), m_now);

          return std::abs(exact_time - time);
        }
      else
        {
          return MyProblem::get_error_norm();
        }
    }

    /// ??ds this is a total hack, I should integrate and divide by volume
    /// or something, but instead I've just averaged it element-wise..
    Vector<double> average_magnetostatic_field() const
    {
      // ??ds implement..
      return Vector<double>();
    }

    // Access functions
    // ============================================================

    unsigned m_index(const unsigned &j) const
    {return this->ele_pt()->m_index_micromag(j);}

    unsigned phi_index() const
    {
      return this->ele_pt()->phi_index_micromag();
    }

    unsigned phi_1_index() const
    {
      return this->ele_pt()->phi_1_index_micromag();
    }

    /// \short Set function for Magnetic_parameters_pt.
    void set_mag_parameters_pt(MagneticParameters* _magnetic_parameters_pt)
    {
      Magnetic_parameters_pt = _magnetic_parameters_pt;
    }

    /// \short Const access function for Magnetic_parameters_pt.
    const MagneticParameters* mag_parameters_pt() const
    {
#ifdef PARANOID
      if(Magnetic_parameters_pt == 0)
        {
          std::string error_msg = "Magnetic paramters pointer is null.";
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return Magnetic_parameters_pt;
    }

    /// \short Get a pointer to a MicromagEquations element.
    MicromagEquations* ele_pt() const
    {
      return checked_dynamic_cast<MicromagEquations*>(mesh_pt(0)->element_pt(0));
    }

    BoundaryElementHandlerBase* bem_handler_pt() const
    {
#ifdef PARANOID
      if(Bem_handler_pt == 0)
        {
          std::string err = "Null Bem_handler_pt.";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif
      return Bem_handler_pt;
    }

    Mesh* flux_mesh_factory(Mesh* mesh_pt,
                            Vector<unsigned>& boundaries) const
    {
#ifdef PARANOID
      if(Flux_mesh_factory_pt == 0)
        {
          throw OomphLibError("No flux mesh factory set!",
                              OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif
      return Flux_mesh_factory_pt(mesh_pt, boundaries);
    }

    bool implicit_ms_flag() const
    {
      return (!Decoupled_ms) && (!Disable_ms)
        && (!Inside_explicit_timestep) && (!Inside_segregated_magnetostatics);
    }

    bool fembem_ms_flag() const
    {
      return !Disable_ms;
    }

    bool analytic_ms_flag() const
    {
#ifdef PARANOID
      if(!Disable_ms)
        {
          std::string err = "This is a weird state... prob won't work.";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif
      return Analytic_ms_fct_pt != 0;
    }


    bool renormalise_each_time_step() const
    {
      return Renormalise_each_time_step;
    }

    /// \short Calculate energies and store them for easy reference
    /// (e.g. for output).
    void calculate_energies(bool calculate_effective_damping=true);

    /// \short Calculate the total (micromagnetic) energy for all meshes in
    /// the problem.
    double micromagnetic_energy() const
    {
      return Exchange_energy + Zeeman_energy +
        Crystalline_anisotropy_energy + Magnetostatic_energy;
    }

    /// Output |m| error values of each node along with spatial position
    /// for plotting with paraview. To plot load the csv file, use "Table
    /// To Points", add a 3D view, set the ltes to visible and color by the
    /// lte of interest. Useful to set "Representation" to "Points" and
    /// increase point size so that things are easily visible. Not
    /// implemented for nodes with varying numbers of values (you might
    /// need something fancier than a csv file for this).
    void doc_m_length_error(const unsigned& t_hist, std::ostream& output) const
    {
      // Output position labels
      for(unsigned j=0; j<Dim; j++)
        {
          output << "x" << j << ", ";
        }

      // Output label for |m| error
      output << "ml_error" << ", ";

      output << std::endl;

      // Output actual positions and errors
      for(unsigned i=0, ni=mesh_pt()->nnode(); i<ni; i++)
        {
          Node* nd_pt = mesh_pt()->node_pt(i);

          // Output position of node
          for(unsigned j=0; j<Dim; j++)
            {
              output << nd_pt->x(j) << ", ";
            }

          // and the error
          Vector<double> m(3, 0.0);
          for(unsigned j=0; j<3; j++)
            {
              m[j] = nd_pt->value(t_hist, m_index(j));
            }
          output << std::abs(1 - VectorOps::two_norm(m));

          output << std::endl;
        }
    }

    /// Can we check the solution using Mallinson's exact time + phi
    /// solutions?
    bool Compare_with_mallinson;

    bool Disable_ms;

    MagnetostaticFieldFctPt Analytic_ms_fct_pt;

    bool Check_angles;

    /// Normalise magnetisation problem after each step?
    bool Renormalise_each_time_step;

    bool Pin_boundary_m;
    bool Use_fd_jacobian;

    bool Inside_segregated_magnetostatics;

    LLGResidualCalculator* Residual_calculator_pt;

    /// Should the magnetisation be relaxed before we start time
    /// integration
    bool Relax_magnetisation;

    /// Should we write out spatial values of |m| error
    bool Doc_m_length_error;

  private:

    /// \short Storage for initial linear solver: in case we want to swap
    /// to superlu for large dt.
    LinearSolver* My_linear_solver_pt;
    LinearSolver* Super_LU_solver_pt;

    /// Magnetic parameters storage. ??ds should maybe go in meshes?
    MagneticParameters* Magnetic_parameters_pt;

    /// \short Exchange_energy, computed after previous Newton solve.
    double Exchange_energy;

    /// \short Zeeman energy, computed after previous Newton solve.
    double Zeeman_energy;

    /// \short Crystalline anisotropy energy, computed after previous Newton
    /// solve.
    double Crystalline_anisotropy_energy;

    /// \short Magnetostatic energy, computed after previous Newton solve.
    double Magnetostatic_energy;

    /// \short Mesh for flux elements to impose boundary condition on phi1.
    Mesh* Flux_mesh_pt;

public:

    /// \short Pointer to class for handling BEM
    BoundaryElementHandlerBase* Bem_handler_pt;

    /// \short Pointer to function for creating the flux mesh (can't be
    /// hard coded becuase it depends on the element type, which depends on
    /// dimension etc.)
    FluxMeshFactoryFctPt Flux_mesh_factory_pt;

    /// Pointer to function deteriminig values of Dirichlet boundaries.
    InitialMFct* Boundary_solution_pt;

    /// \short Error in recomputed damping constant for the last time step
    /// (based on actual change in energy).
    double Rel_damping_error;
    double Abs_damping_error;
    std::deque<double> Previous_energies;

    /// Enumeration for how to handle phi_1's singularity
    phi_1_singularity_handling::phi_1_singularity_handling Phi_1_singularity_method;

    /// Should we pin a phi1 value
    bool pin_a_bulk_phi_1() const
      {
        return Phi_1_singularity_method == phi_1_singularity_handling::pin_bulk;
      }

    /// Should we pin a boundary phi1 value (e.g. for testing that this works)
    bool pin_a_boundary_phi_1() const
    {
      return Phi_1_singularity_method == phi_1_singularity_handling::pin_boundary;
    }

    /// Should we pin a phi1 value (with boundary values allowed)
    bool pin_any_phi_1() const
    {
      return Phi_1_singularity_method == phi_1_singularity_handling::pin_any;
    }

    /// Should we normalise phi1 values
    bool normalise_phi_1() const
    {
      return Phi_1_singularity_method == phi_1_singularity_handling::normalise;
    }

    bool twiddle_phi_1_jacobian() const
      {
        return Phi_1_singularity_method == phi_1_singularity_handling::jacobian;
      }


  private:

    /// Inaccessible copy constructor
    LLGProblem(const LLGProblem & dummy)
    {BrokenCopy::broken_copy("LLGProblem");}

    /// Inaccessible assignment operator
    void operator=(const LLGProblem &dummy)
    {BrokenCopy::broken_assign("LLGProblem");}


    // Decoupled ms code
    // ============================================================

  public:

    /// Are we solving for ms "properly" or using a separate solve?
    bool Decoupled_ms;

    /// Should we extrapolate it to the correct time using history?
    bool Extrapolate_decoupled_ms;

    /// Boolean to be set in explicit predictor steps to avoid
    /// get_jacobian(..) getting the BEM wrapped up in its mass
    /// matrix. ??ds should probably modify how explicit time stepping
    /// works instead to avoid this.
    bool Inside_explicit_timestep;

    BEMElementFactoryFctPt Bem_element_factory_pt;

    /// Storage for nodes to unpin after explicit step
    Vector<Node*> unpinned_phi_nodes;
    Vector<Node*> unpinned_phi_1_nodes;

    /// Solver parameter sets to use in segregated phi/phi1 solves.
    SolverParameters Phi_seg_solve_parameters;
    SolverParameters Phi_1_seg_solve_parameters;

    /// Bool to disable optimisations on phi/phi1 segregated solves (linear,
    /// constant Jacobian, CG w/ AMG solver) for debugging.
    bool Disable_magnetostatic_solver_optimistations;

    /// Function to use to create integration schemes for elements.
    NodalQuadratureFactoryFctPt Nodal_quadrature_factory_fpt;

    /// Force energy calculations to ignore the integral_pt() and use
    /// Gaussian quadrature.
    bool Force_gaussian_quadrature_in_energy;


    /// \short Solve for the magnetostatic field.
    void magnetostatics_solve();

    /// Linearly extrapolate phi
    void extrapolate_phi(const double& new_dt, const double& prev_dt);


  };

}

#endif
