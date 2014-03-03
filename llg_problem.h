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
#include "residual_swapping_explicit_timestepper.h"
#include "generic_poisson_problem.h"
#include "micromag_types.h"

namespace oomph
{

  class ResidualCalculator;

  // ============================================================
  ///
  // ============================================================
  class LLGProblem : public MyProblem
  {
  public:

    /// Default constructor - do nothing except nulling pointers.
    LLGProblem() :
      Compare_with_mallinson(false),
      Applied_field_fct_pt(0),
      Previous_energies(5, 0.0)
    {
      // Needed for if we want to switch solvers in runs
      Super_LU_solver_pt = new SuperLUSolver;

      // Storage for residual calculator
      Residual_calculator_pt = 0;

      // Initialise storage for energies etc.
      Exchange_energy = MyProblem::Dummy_doc_data;
      Zeeman_energy = MyProblem::Dummy_doc_data;
      Crystalline_anisotropy_energy = MyProblem::Dummy_doc_data;
      Magnetostatic_energy = MyProblem::Dummy_doc_data;
      Effective_damping_constant = MyProblem::Dummy_doc_data;
      Alt_eff_damp = MyProblem::Dummy_doc_data;

      // Bem stuff
      Bem_handler_pt = 0;
      Flux_mesh_pt = 0;
      Flux_mesh_factory_pt = 0;

      Decoupled_ms = false;
      Extrapolate_decoupled_ms = false;
      Disable_ms = false;
      Analytic_ms_fct_pt = 0;
      Inside_explicit_timestep = false;
#ifdef PARANOID
      Check_angles = true;
#else
      Check_angles = false;
#endif
      Phi_problem_pt = 0;
      Phi_1_problem_pt = 0;

      // Debugging switches
      Pin_boundary_m = false;
      Use_fd_jacobian = false;
      Renormalise_each_time_step = -1;
    }

    std::string problem_name() const {return "LLG";}

    /// Get the jacobian as a SumOfMatrices. This is probably the best way
    /// to deal with Jacobians involving a dense block (i.e. in fully
    /// implicit bem). If we aren't using that then this is basically the
    /// same as using the cr matrix form but with the Jacobian wrapped
    /// inside a SumOfMatrices class.
    void get_jacobian(DoubleVector &residuals, SumOfMatrices &jacobian)
    {
      // Get the fem Jacobian and the residuals
      CRDoubleMatrix* fem_jacobian_pt =
        new CRDoubleMatrix(residuals.distribution_pt());
      Problem::get_jacobian(residuals, *fem_jacobian_pt);

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
          jacobian.add_matrix(Bem_handler_pt->bem_matrix_pt(),
                              Bem_handler_pt->row_lookup_pt(),
                              Bem_handler_pt->col_lookup_pt(),
                              false);

          // Create identity CRDoubleMatrix
          unsigned bem_nnode = Bem_handler_pt->bem_mesh_pt()->nnode();
          unsigned nrow = fem_jacobian_pt->nrow();
          LinearAlgebraDistribution* dist_pt = fem_jacobian_pt->distribution_pt();
          Vector<double> values(bem_nnode, -1.0);
          Vector<int> row_index(bem_nnode), row_start;
          for(unsigned nd=0; nd<bem_nnode; nd++)
            {
              unsigned i = Bem_handler_pt->output_lookup_pt()
                ->added_to_main(nd);
              row_index[nd] = i;
            }

          // Sort the row index list then create row starts vector. DO NOT
          // COPY THIS CODE FOR CREATION OF OTHER MATRICES: sorting is only
          // safe because all the values are the same (1.0).
          std::sort(row_index.begin(), row_index.end());
          VectorOps::rowindex2rowstart(row_index, nrow, row_start);

          // Create the matrix
          CRDoubleMatrix bem_block_identity(dist_pt, nrow, values,
                                            row_index, row_start);


          // Add it on
          VectorOps::cr_matrix_add(*fem_jacobian_pt, bem_block_identity,
                                   *fem_jacobian_pt);
        }

    }

    /// Get the Jacobian as a CRDoubleMatrix (the normal matrix format). If
    /// we are using fully implicit bem then this is not a good idea for
    /// "real" problems but useful for tests. Otherwise this is exactly the
    /// same as Problem::get_jacobian().
    void get_jacobian(DoubleVector &residuals, CRDoubleMatrix &jacobian)
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
    }

    /// Function that does the real work of the constructors.
    void build(Vector<Mesh*>& bulk_mesh_pts);

    /// Destructor
    virtual ~LLGProblem()
    {
      // mesh is cleaned up by problem base class
      // timestepper is cleaned up by problem base class
      delete Magnetic_parameters_pt; Magnetic_parameters_pt = 0;
      delete Residual_calculator_pt; Residual_calculator_pt = 0;
      delete Bem_handler_pt; Bem_handler_pt = 0;
      delete Phi_problem_pt; Phi_problem_pt = 0;
      delete Phi_1_problem_pt; Phi_1_problem_pt = 0;

      // Kill boundary value storage vectors
      for(unsigned j=0; j<Phi_boundary_values_pts.size(); j++)
        {
          delete Phi_boundary_values_pts[j];
        }
    }

    /// Renormalise magnetisation to 1 (needed with BDF2)
    void renormalise_magnetisation()
    {
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


    virtual void actions_before_implicit_timestep()
    {
      MyProblem::actions_before_implicit_timestep();

      // // ??ds adaptive midpoint doesn't restore phi state properly!  so
      // // recalculate ms. This needs to be done more efficiently eventually
      // if(Decoupled_ms)
      //   {
      //     // Solve for the magnetostatic field.
      //     magnetostatics_solve();
      // }

      if(Decoupled_ms && Extrapolate_decoupled_ms)
        {
          // Project to correct time
          extrapolate_phi(time_stepper_pt()->time_pt()->dt(),
                          time_stepper_pt()->time_pt()->dt(1));
        }

    }

    /// Overload shift time values to also shift phi problem if needed.
    void shift_time_values()
    {
      Problem::shift_time_values();

      if(Decoupled_ms)
        {
          // Push old phi values back in time (so that we can use them
          // later to get time derivatives of the field).
          phi_problem_pt()->shift_time_values();
        }
    }

    virtual void actions_after_implicit_timestep()
      {
        MyProblem::actions_after_implicit_timestep();

        if(Decoupled_ms)
          {
            // Solve for the magnetostatic field.
            magnetostatics_solve();
          }
      }


    virtual void actions_before_newton_step()
    {
      // Call base class version
      MyProblem::actions_before_newton_step();

#ifdef PARANOID
      if(Bem_handler_pt == 0 && !Disable_ms)
        {
          std::string err = "No bem handler pointer set but magnetostatics are enabled";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif

      oomph_info << std::endl
                << "Newton step " << Nnewton_iter_taken + 1 << std::endl
                << "---------------------------------------" << std::endl;
    }


    virtual void actions_before_newton_solve()
    {
#ifdef PARANOID
      if(Decoupled_ms && !Decoupled_ms_has_been_calculated)
        {
          std::string err = "Decoupled_ms_has_been_calculated not set";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif

      // Call base class version
      MyProblem::actions_before_newton_solve();

      // Update BEM magnetostatics boundary conditions (if we are doing
      // them fully implicitly).
      maybe_update_bem_boundary_conditions();
    }

    virtual void actions_before_explicit_timestep()
    {
      MyProblem::actions_before_explicit_timestep();

      // If phi values are in the main dofs then we need to pin them, they
      // can't be explicitly timestepped in oomph-lib's framework!
      if(!Decoupled_ms)
        {
          const unsigned phi_index = this->phi_index();
          const unsigned phi_1_index = this->phi_1_index();

          // Preallocate enough storage to make sure we don't do it
          // repeatedly inside the loop.
          const unsigned nnode = mesh_pt()->nnode();
          unpinned_phi_nodes.reserve(nnode);
          unpinned_phi_1_nodes.reserve(nnode);


          // Store list of pointers to nodes with unpinned phi/phi1 and pin
          // their phi values.
          for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
            {
              Mesh* mesh_pt = this->mesh_pt(msh);
              for(unsigned nd=0, nnd=mesh_pt->nnode(); nd<nnd; nd++)
                {
                  Node* nd_pt = mesh_pt->node_pt(nd);

                  if(!nd_pt->is_pinned(phi_index))
                    {
                      unpinned_phi_nodes.push_back(nd_pt);
                      nd_pt->pin(phi_index);
                    }

                  if(!nd_pt->is_pinned(phi_1_index))
                    {
                      unpinned_phi_1_nodes.push_back(nd_pt);
                      nd_pt->pin(phi_1_index);
                    }
                }
            }

          // reassign equation numbers
          std::cout << assign_eqn_numbers() << std::endl;
        }

      // Set this variable to avoid getting BEM in mass matrix (due to a
      // hack in oomph core.. fix that instead?)
      Inside_explicit_timestep = true;
    }

    virtual void actions_after_explicit_timestep()
      {
        MyProblem::actions_after_explicit_timestep();

        // Need to unpin any phi that we pinned earlier
        if(!Decoupled_ms)
          {
            const unsigned phi_index = this->phi_index();
            const unsigned phi_1_index = this->phi_1_index();

            // unpin phi values
            for(unsigned j=0; j<unpinned_phi_nodes.size(); j++)
              {
                unpinned_phi_nodes[j]->unpin(phi_index);
              }

            // unpin phi 1 values
            for(unsigned j=0; j<unpinned_phi_1_nodes.size(); j++)
              {
                unpinned_phi_1_nodes[j]->unpin(phi_1_index);
              }

            // reassign equation numbers
            std::cout << assign_eqn_numbers() << std::endl;

          }

        // We need to keep M normalised...
        oomph_info << "Renormalising nodal magnetisations." << std::endl;
        renormalise_magnetisation();

        // check neighbouring magnetisation angles if requested
        maybe_check_angles();

        // Explicit timestep is over now
        Inside_explicit_timestep = false;
      }

    virtual void actions_before_explicit_stage()
    {
#ifdef PARANOID
      if(Decoupled_ms && !Decoupled_ms_has_been_calculated)
        {
          std::string err = "Decoupled_ms_has_been_calculated not set";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif

      MyProblem::actions_before_explicit_stage();
    }

    virtual void actions_after_explicit_stage()
    {
      // Solve for the new magnetostatic field.
      if(Decoupled_ms)
        {
          magnetostatics_solve();
        }

      MyProblem::actions_after_explicit_stage();
    }

    virtual void actions_after_newton_solve()
    {
      MyProblem::actions_after_newton_solve();

      oomph_info << std::endl
                << "Finalising" << std::endl
                << "-------------------------------------------" << std::endl;

      // If we're using BDF we need to keep M normalised.
      if(renormalise_each_time_step())
        {
          oomph_info << "Renormalising nodal magnetisations." << std::endl;
          renormalise_magnetisation();
        }

      // check neighbouring magnetisation angles if requested
      maybe_check_angles();

      // Calculate and store the energy (ready to be output)
      calculate_energies();
    }

    void maybe_check_angles()
      {
        if(Check_angles)
          {
            std::cout << "angles" << std::endl;
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

    virtual void actions_before_newton_convergence_check()
      {
        // Call base class actions function
        MyProblem::actions_before_newton_convergence_check();

        // Update BEM magnetostatics boundary conditions (if we are doing them
        // fully implicitly).
        maybe_update_bem_boundary_conditions();
      }

    void maybe_update_bem_boundary_conditions()
      {
        if(implicit_ms_flag())
          {
            Bem_handler_pt->get_bem_values_and_copy_into_values();
          }

        // Otherwise don't do anything, if bem is being used
        // semi-implicitly then it's done elsewhere. If it's not being used
        // then nothing to update.
      }


    /// Output solution
    void doc_solution_additional(std::ofstream &some_file) const
    {
      // Number of plot points
      unsigned npts = 2;

      // Output solution with specified number of plot points per element
      mesh_pt()->output(some_file, npts);

      if(Decoupled_ms)
        {
          // Output the magnetostatic field data
          std::ofstream field_file((Doc_info.directory() + "/field"
                                    + Doc_info.number_as_string() + ".dat").c_str());
          phi_problem_pt()->mesh_pt()->output(field_file, npts);
          field_file.close();

          std::ofstream phi1_file((Doc_info.directory() + "/phione"
                                   + Doc_info.number_as_string() + ".dat").c_str());
          phi_1_problem_pt()->mesh_pt()->output(phi1_file, npts);
          phi1_file.close();
        }
    }

    void write_additional_trace_headers(std::ofstream& trace_file) const
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
        << Trace_seperator << "effective_damping"
        << Trace_seperator << "alt_effective_damping"
        << Trace_seperator << "h_applied_first_element"
        ;

    }

    void write_additional_trace_data(std::ofstream& trace_file) const
    {

      // Get average (and standard deviation) of |m| - 1
      double m_error_avg(0), m_error_stddev(0);
      norm_m_error(m_error_avg, m_error_stddev);

      Vector<double> angle_variations = elemental_max_m_angle_variations();
      Vector<double> mean_m = mean_magnetisation();
      Vector<double> h_app = first_element_happ();

      trace_file
        << Trace_seperator << m_error_avg
        << Trace_seperator << m_error_stddev
        << Trace_seperator
        << *std::max_element(angle_variations.begin(), angle_variations.end())
        << Trace_seperator << mean_m[0]
        << Trace_seperator << mean_m[1]
        << Trace_seperator << mean_m[2]

        << Trace_seperator << Exchange_energy
        << Trace_seperator << Zeeman_energy
        << Trace_seperator << Crystalline_anisotropy_energy
        << Trace_seperator << Magnetostatic_energy
        << Trace_seperator << micromagnetic_energy()
        << Trace_seperator << Effective_damping_constant
        << Trace_seperator << Alt_eff_damp
        << Trace_seperator << h_app
        ;
    }

    /// Get happ value from first element
    Vector<double> first_element_happ() const
    {
      double t = time_stepper_pt()->time();
      Vector<double> H_app;
      Vector<double> dummy_x, dummy_s;
      ele_pt()->get_applied_field(t, dummy_x, dummy_s, H_app);

      return H_app;
    }

    void initial_doc_additional() const
    {
      // If we have a H matrix then write out its rank data
      if(Bem_handler_pt != 0)
        {
          Bem_handler_pt->maybe_write_h_matrix_data(Doc_info.directory());
        }
    }

    /// \short Return a vector of the maximum angle variation in each
    /// element.
    Vector<double> elemental_max_m_angle_variations() const
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

              angles.push_back(ele_pt->max_m_angle_variation());
            }
        }
      return angles;
    }

    /// \short Error for adaptive timestepper (rms of nodal error determined by
    /// comparison with explicit timestepper result).
    double global_temporal_error_norm();


    virtual void actions_after_set_initial_condition()
      {
        MyProblem::actions_after_set_initial_condition();

        if(Decoupled_ms)
          {
            // Solve for initial field and phi values
            magnetostatics_solve();

            // Copy phi backwards in time
            phi_problem_pt()->set_up_impulsive_initial_condition();
          }

        calculate_energies(false);
      }

    virtual void dump(std::ofstream& dump_file) const
    {
      // Set very high precision to avoid any issues
      dump_file.precision(14);

      // Dump sub problems if needed
      if(Decoupled_ms)
        {
          phi_problem_pt()->dump(dump_file);
          phi_1_problem_pt()->dump(dump_file);
        }

      // Let base class handle the rest
      MyProblem::dump(dump_file);
    }

    virtual void read(std::ifstream& restart_file)
    {

      // Read sub problems if needed
      if(Decoupled_ms)
        {
          phi_problem_pt()->read(restart_file);
          phi_1_problem_pt()->read(restart_file);
        }

      // Let base class handle the rest
      MyProblem::read(restart_file);


      actions_after_set_initial_condition();
    }



    // Lots of magnetisation manipulation functions
    // ============================================================

    /// \short Loop over all nodes in bulk mesh and get magnetisations
    Vector<Vector<double> > get_nodal_magnetisations(unsigned i_time=0) const;

    /// \short Abs of mean difference of actual m and m given by a function
    /// at the middle of each element.
    double compare_m_with_function(const InitialM::InitialMFctPt fct_pt) const;

    /// \short ??ds
    double mean_norm_m_error() const;

    /// \short ??ds
    Vector<double> mean_magnetisation() const;

    /// Best trace value to use is probably mean of m.
    Vector<double> trace_values() const {return mean_magnetisation();}

    /// \short ??ds
    void get_nodal_two_norms(Vector<double> &output) const
    {
      Vector< Vector<double> > m = get_nodal_magnetisations();
      output.assign(m.size(), 0.0);
      transform(m.begin(), m.end(), output.begin(), VectorOps::two_norm);
    }

    /// \short ??ds
    double mean_nodal_magnetisation_length() const
    {
      Vector<double> ms; get_nodal_two_norms(ms);
      return VectorOps::mean(ms);
    }

    /// \short ??ds
    void norm_m_error(double &m_error_avg, double &m_error_stddev) const;


    double get_error_norm() const
    {
      if(Compare_with_mallinson)
        {
          using namespace CompareSolutions;

          double time = ele_pt()->node_pt(0)->time_stepper_pt()->time();
          Vector<double> m_now = mean_magnetisation();
          double exact_time = switching_time_wrapper(mag_parameters_pt(), m_now);

          return std::abs(exact_time - time);
        }
      else
        {
          return MyProblem::Dummy_doc_data;
        }
    }

    /// ??ds this is a total hack, I should integrate and divide by volume
    /// or something, but instead I've just averaged it element-wise..
    Vector<double> average_magnetostatic_field() const
    {
      if(Decoupled_ms)
        {
          const unsigned nodal_dim = checked_dynamic_cast<MagnetostaticFieldEquations*>
            (phi_problem_pt()->mesh_pt()->element_pt(0))->node_pt(0)->ndim();

          // Pick a point in the middle of the element
          const Vector<double> s(nodal_dim, 0.3);
          Vector<double> total_ms(3, 0.0);

          // // Loop over all elements calculating the value in the middle of the element
          // for(unsigned e=0, ne=phi_problem_pt()->mesh_pt()->nelement(); e < ne; e++)
          //   {
          //     MagnetostaticFieldEquations* ele_pt
          //       = checked_dynamic_cast<MagnetostaticFieldEquations*>
          //       (phi_problem_pt()->mesh_pt()->element_pt(e));

          //     // Get the shape function and eulerian coordinate derivative at
          //     // position s.
          //     unsigned n_node = ele_pt->nnode();
          //     Shape psi(n_node); DShape dpsidx(n_node,nodal_dim);
          //     ele_pt->dshape_eulerian(s,psi,dpsidx);

          //     // Interpolate grad phi
          //     Vector<double> interpolated_dphidx(nodal_dim,0.0);
          //     for(unsigned l=0;l<n_node;l++)
          //       {
          //         double phi_value = ele_pt->raw_nodal_value(l,0);
          //         for(unsigned i=0; i<nodal_dim; i++)
          //           {interpolated_dphidx[i] += phi_value*dpsidx(l,i);}
          //       }

          //     // Add this grad phi to the sum
          //     for(unsigned j=0; j<nodal_dim; j++)
          //       {
          //         total_dphidx[j] += interpolated_dphidx[j];
          //       }
          //   }

          // Loop over all elements calculating the value in the middle of the element
          for(unsigned e=0, ne=mesh_pt()->nelement(); e < ne; e++)
            {
              MicromagEquations* ele_pt = checked_dynamic_cast<MicromagEquations*>
                (mesh_pt()->element_pt(e));

              // Interpolate
              Vector<double> ms;
              ele_pt->get_magnetostatic_field(s, ms);

              // Add this to the sum
              for(unsigned j=0; j<3; j++)
                {
                  total_ms[j] += ms[j];
                }
            }

          // Divide sum by number of elements to get the average. Take the
          // negative to get the field.
          double nele = double(phi_problem_pt()->mesh_pt()->nelement());
          Vector<double> average_magnetostatic_field(3,0.0);
          for(unsigned j=0; j<nodal_dim; j++)
            {
              average_magnetostatic_field[j] = total_ms[j] / nele;
            }

          return average_magnetostatic_field;
        }
      else
        {
          throw OomphLibError("Function not yet implemented for fully coupled ms",
                              OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
        }


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

    /// \short Non-const access function for Applied_field_fct_pt.
    HApp::HAppFctPt& applied_field_fct_pt() {return Applied_field_fct_pt;}

    /// \short Const access function for Applied_field_fct_pt.
    HApp::HAppFctPt applied_field_fct_pt() const {return Applied_field_fct_pt;}

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
        return (!Decoupled_ms) && (!Disable_ms) && (!Inside_explicit_timestep);
      }


    bool renormalise_each_time_step() const
    {
      // If flag not set then only do it for bdf timesteppers
      if(Renormalise_each_time_step == -1)
        {
          return ((dynamic_cast<const BDF<2>*>(time_stepper_pt()) != 0)
                  || (dynamic_cast<const BDF<1>*>(time_stepper_pt()) != 0));
        }
      // Otherwise do what the flag says
      else if((Renormalise_each_time_step == +1) || (Renormalise_each_time_step == 0))
        {
          return bool(Renormalise_each_time_step);
        }
      else
        {
          std::string err = "Undefined state for renormalise ";
          err += to_string(Renormalise_each_time_step);
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
    }

    /// \short Calculate energies and store them for easy reference
    /// (e.g. for output).
    void calculate_energies(bool calculate_effective_damping=true)
    {
      // If you want to turn off energy calculations (e.g. for speed)
      // this is the place to do it. Replace values with
      // MyProblem::Dummy_doc_data.

      // Calculate and store new values
      Exchange_energy = MManipulation::exchange_energy(this);
      Zeeman_energy = MManipulation::zeeman_energy(this);
      Crystalline_anisotropy_energy =
        MManipulation::crystalline_anisotropy_energy(this);
      Magnetostatic_energy = MManipulation::magnetostatic_energy(this);

      // Store energy for damping calculations
      Previous_energies.push_front(micromagnetic_energy());

      // Keep the list of previous energies reasonably small (we only
      // need N for any bdf<N> calculation).
      if(Previous_energies.size() > 5) Previous_energies.pop_back();

      // Calculate and store effective damping if not disabled.
      if(calculate_effective_damping)
        {
          Effective_damping_constant =
            MManipulation::effective_damping_used(this);
          Alt_eff_damp = MManipulation::
            alt_effective_damping_used(this, Previous_energies);
        }
    }


    /// \short Calculate the total (micromagnetic) energy for all meshes in
    /// the problem.
    double micromagnetic_energy() const
    {
      return Exchange_energy + Zeeman_energy +
        Crystalline_anisotropy_energy + Magnetostatic_energy;
    }

    /// Can we check the solution using Mallinson's exact time + phi
    /// solutions?
    bool Compare_with_mallinson;

    bool Disable_ms;

    MagnetostaticFieldFctPt Analytic_ms_fct_pt;

    bool Check_angles;

    /// Normalise magnetisation problem after each step?
    int Renormalise_each_time_step;

    bool Pin_boundary_m;
    bool Use_fd_jacobian;

    LLGResidualCalculator* Residual_calculator_pt;

  private:

    /// \short Storage for initial linear solver: in case we want to swap
    /// to superlu for large dt.
    LinearSolver* My_linear_solver_pt;
    LinearSolver* Super_LU_solver_pt;

    /// Pointer to the applied field.
    HApp::HAppFctPt Applied_field_fct_pt;

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
    BoundaryElementHandler* Bem_handler_pt;

    /// \short Pointer to function for creating the flux mesh (can't be
    /// hard coded becuase it depends on the element type, which depends on
    /// dimension etc.)
    FluxMeshFactoryFctPt Flux_mesh_factory_pt;

    /// \short Recomputed effective damping constant for the last time step
    /// (based on actual change in energy).
    double Effective_damping_constant;
    double Alt_eff_damp;
    std::deque<double> Previous_energies;

  private:

    /// Inaccessible copy constructor
    LLGProblem(const LLGProblem & dummy)
    {BrokenCopy::broken_copy("LLGProblem");}

    /// Inaccessible assignment operator
    void operator=(const LLGProblem &dummy)
    {BrokenCopy::broken_assign("LLGProblem");}


    // Decoupled ms code
    // ============================================================

  private:

    /// Sub problems for the magnetostatics solve
    GenericPoissonProblem* Phi_1_problem_pt;
    GenericPoissonProblem* Phi_problem_pt;

    /// Intermediate storage for results of bem (ideally we would have it
    /// call a function to get the boundary values filled in but c++ member
    /// functions pointers are useless...)
    Vector<DoubleVector*> Phi_boundary_values_pts;


    /// Has decoupled ms been calculated yet? If not then we need to
    /// calculate it before any steps.
    bool Decoupled_ms_has_been_calculated;

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

    GenericPoissonProblem::FluxMeshFactoryFctPt Phi_1_flux_mesh_factory_fct_pt;

    /// Storage for nodes to unpin after explicit step
    Vector<Node*> unpinned_phi_nodes;
    Vector<Node*> unpinned_phi_1_nodes;

    void build_decoupled_ms(Vector<Mesh*>& llg_mesh_pts,
                            Vector<Mesh*>& phi_mesh_pts,
                            Vector<Mesh*>& phi_1_mesh_pts);

    /// \short Solve for the magnetostatic field.
    void magnetostatics_solve();

    GenericPoissonProblem* phi_1_problem_pt() const
    {
#ifdef PARANOID
      if(Phi_1_problem_pt == 0)
        {
          std::string error_msg = "Phi 1 problem pointer is null!";
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return Phi_1_problem_pt;
    }

    /// Linearly extrapolate phi
    void extrapolate_phi(const double& new_dt, const double& prev_dt);

    void set_phi_1_problem_pt(GenericPoissonProblem* p)
    { Phi_1_problem_pt = p;}

    GenericPoissonProblem* phi_problem_pt() const
    {
#ifdef PARANOID
      if(Phi_problem_pt == 0)
        {
          std::string error_msg = "Phi problem pointer is null!";
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return Phi_problem_pt;
    }

    void set_phi_problem_pt(GenericPoissonProblem* p)
    { Phi_problem_pt = p;}


  };


  /// \short Functions for constructing things needed for LLG problems
  namespace LLGFactories
  {
    /// \short Make a mesh as specified by an input argument. Refined
    /// according to the given refinement level (in some way appropriate
    /// for that mesh type). Assumption: this will be passed into a
    /// problem, which will delete the pointer when it's done.
    Mesh* mesh_factory(const std::string& _mesh_name,
                       int refinement_level,
                       TimeStepper* time_stepper_pt,
                       double scaling_factor=1.0,
                       unsigned nnode1d = 2);

    LLGResidualCalculator* residual_calculator_factory(const std::string& residual);

    /// \short Create a variable order quadrature object based on the
    /// dimension and shape of the element. Only works for some element
    /// types.
    Integral* variable_order_integrator_factory(const FiniteElement* const el_pt);

    /// \short Create a function to create bem elements based on the
    /// elements used in the bulk mesh.
    BEMElementFactoryFctPt bem_element_factory_factory
    (const FiniteElement* bulk_ele_pt);

    /// \short very simple function: create a new face element of type
    /// ELEMENT.
    template<class ELEMENT>
    MicromagBEMElementEquations* bem_element_factory(FiniteElement* ele,
                                                     const int& face)
    {
      return new ELEMENT(ele, face);
    }

    FluxMeshFactoryFctPt
    mm_flux_mesh_factory_factory(const FiniteElement* bulk_ele_pt);
  }

  /// \short A namespace full of functions that take some "dynamic"
  /// (i.e. can be calculated at runtime) input and create a new instance
  /// of the appropriate object, using the new command (Factory Method
  /// design pattern).
  ///
  /// Typically these objects are passed straight into other classes and
  /// will be deleted by the destructor of that class. If not it is your
  /// responsibility to make sure the objects are deleted.
  namespace SemiImplicitFactories
  {
    /// \short Make a mesh of Micromag elements as specified by an
    /// input argument. Refined according to the given refinement level (in
    /// some way appropriate for that mesh type).
    Mesh* llg_mesh_factory(const std::string& _mesh_name,
                           int refinement_level,
                           TimeStepper* time_stepper_pt,
                           double scaling_factor = 1.0,
                           unsigned nnode1d = 2);


    /// \short Make a mesh of MagnetostaticField elements as specified by an
    /// input argument. Refined according to the given refinement level (in
    /// some way appropriate for that mesh type).
    Mesh* phi_mesh_factory(const std::string& _mesh_name,
                           int refinement_level,
                           TimeStepper* time_stepper_pt,
                           double scaling_factor = 1.0,
                           unsigned nnode1d = 2);


    /// \short Return a factory function which will create the appropriate
    /// "flux mesh" for the bulk element pointer given.
    GenericPoissonProblem::FluxMeshFactoryFctPt
    phi_1_flux_mesh_factory_factory(const FiniteElement* bulk_phi_1_ele_pt);
  }


  /// Base command line args processing class for pure llg and semi
  /// implicit (and any other magnetism problems).
  class MMArgs : public MyCliArgs
  {
  public:
    /// Constructor: Initialise pointers to null.
    MMArgs() : h_app_fct_pt(0), mag_params_pt(0) {}

    virtual void set_flags()
    {
      MyCliArgs::set_flags();

      specify_command_line_flag("-initial-m", &initial_m_name);
      initial_m_name = "z";

      specify_command_line_flag("-h-app", &h_app_name);
      h_app_name = "minus_z";

      specify_command_line_flag("-mag-params", &mag_params_name);
      mag_params_name = "simple-llg";

      specify_command_line_flag("-renormalise", &renormalise);
      renormalise = -1;

      specify_command_line_flag("-damping", &damping);
      damping = -10;

      // Flags automatically default to false
      specify_command_line_flag("-pin-boundary-m");

      specify_command_line_flag("-hlib-bem", &hlib_bem);
      hlib_bem = -1;

      specify_command_line_flag("-numerical-int-bem", &numerical_int_bem);
      numerical_int_bem = -1;

      specify_command_line_flag("-mallinson", &mallinson);
      mallinson = -1;

      specify_command_line_flag("-ms-method", &ms_method);
      ms_method = "implicit";

      specify_command_line_flag("-check-angles", &check_angles);
      check_angles = -1;
    }

    bool is_decoupled(const std::string& ms_method) const
    {
      return to_lower(ms_method) == "decoupled"
        || to_lower(ms_method) == "decoupled-no-extrapolation";
    }


    virtual void run_factories()
    {
      using namespace SemiImplicitFactories;
      using namespace LLGFactories;
      using namespace Factories;

      // Figure out how to build meshes
      if(is_decoupled(ms_method))
        {
          mesh_factory_pt = &llg_mesh_factory;
        }
      else
        {
          mesh_factory_pt = &mesh_factory;
        }


      MyCliArgs::run_factories();

      initial_m_name = to_lower(initial_m_name);
      h_app_name = to_lower(h_app_name);
      mag_params_name = to_lower(mag_params_name);

      initial_condition_fpt = InitialM::initial_m_factory(initial_m_name);
      h_app_fct_pt = HApp::h_app_factory(h_app_name);
      mag_params_pt = magnetic_parameters_factory(mag_params_name);

      if(command_line_flag_has_been_set("-damping"))
        {
          mag_params_pt->Gilbert_damping = damping;
        }

      // Copy flags into bools in this class
      pin_boundary_m = command_line_flag_has_been_set("-pin-boundary-m");

      if(is_decoupled(ms_method))
        {
          // Pick the factory function for creating the phi 1 surface mesh
          phi_1_flux_mesh_factory_fct_pt = phi_1_flux_mesh_factory_factory
            (phi_1_mesh_pts[0]->finite_element_pt(0));

          // Pick the factory function for creating the BEM elements
          bem_element_factory_fct_pt = bem_element_factory_factory
            (mesh_pts[0]->finite_element_pt(0));
        }

    }

    void build_meshes()
      {
        // Build the main mesh(es)
        MyCliArgs::build_meshes();

        if(is_decoupled(ms_method))
          {
            // Time stepper for phi, store history values for derivative
            // calculations and extrapolation.
            TimeStepper* phi_time_stepper_pt = new Steady<2>;

            // Time stepper for phi1, don't store history values.
            TimeStepper* phi1_time_stepper_pt = new Steady<0>;


            // Also build separate poisson meshes if needed
            using namespace SemiImplicitFactories;
            phi_mesh_pts = build_meshes_helper(phi_mesh_factory,
                                               phi_time_stepper_pt);
            phi_1_mesh_pts = build_meshes_helper(phi_mesh_factory,
                                                 phi1_time_stepper_pt);
          }
      }

    virtual void assign_specific_parameters(MyProblem* problem_pt) const
    {
      LLGProblem* llg_pt = checked_dynamic_cast<LLGProblem*>(problem_pt);
      llg_pt->applied_field_fct_pt() = h_app_fct_pt;
      llg_pt->set_mag_parameters_pt(mag_params_pt);
      llg_pt->Renormalise_each_time_step = renormalise;
      llg_pt->Pin_boundary_m = pin_boundary_m;

      if(to_lower(ms_method) == "implicit")
        {
          llg_pt->Decoupled_ms = false;
          llg_pt->Disable_ms = false;
          llg_pt->Extrapolate_decoupled_ms = false;
        }
      else if(to_lower(ms_method) == "decoupled-no-extrapolation")
        {
          llg_pt->Decoupled_ms = true;
          llg_pt->Disable_ms = false;
          llg_pt->Extrapolate_decoupled_ms = false;
        }
      else if(to_lower(ms_method) == "decoupled")
        {
          llg_pt->Decoupled_ms = true;
          llg_pt->Disable_ms = false;
          llg_pt->Extrapolate_decoupled_ms = true;
        }
      else if(to_lower(ms_method) == "disabled")
        {
          llg_pt->Decoupled_ms = false;
          llg_pt->Disable_ms = true;
          llg_pt->Extrapolate_decoupled_ms = false;
        }
      else
        {
          llg_pt->Decoupled_ms = false;
          llg_pt->Disable_ms = true;
          llg_pt->Extrapolate_decoupled_ms = false;
          llg_pt->Analytic_ms_fct_pt =
            MagnetostaticFieldFunctions::ms_factory(to_lower(ms_method));
        }

      // ??ds this should maybe be a general one?
      llg_pt->Use_fd_jacobian = use_fd_jacobian;

      // Set exact solution if we have one
      if(((h_app_name == "minus_z")
          && (initial_m_name == "z")
          && (mag_params_pt->gilbert_damping() != 0.0)
          && (to_lower(ms_method) == "disabled"))
         || mallinson == 1)
        {
          llg_pt->Compare_with_mallinson = true;
        }

      if(check_angles != -1)
        {
          llg_pt->Check_angles = bool(check_angles);
        }

      llg_pt->Phi_1_flux_mesh_factory_fct_pt = phi_1_flux_mesh_factory_fct_pt;

      llg_pt->Bem_element_factory_pt = bem_element_factory_fct_pt;
    }

    HApp::HAppFctPt h_app_fct_pt;
    MagneticParameters* mag_params_pt;

    // Strings for input to factory functions
    std::string initial_m_name;
    std::string h_app_name;
    std::string mag_params_name;


    /// Flag to control renormalisation of |m| after each step. -1 =
    /// default for timestepper, 0 = off, 1 = on.
    int renormalise;

    double damping;

    int numerical_int_bem;
    int hlib_bem;
    int mallinson;
    int check_angles;

    std::string ms_method;

    bool pin_boundary_m;

    Vector<Mesh*> phi_1_mesh_pts;
    Vector<Mesh*> phi_mesh_pts;

    GenericPoissonProblem::FluxMeshFactoryFctPt phi_1_flux_mesh_factory_fct_pt;
    BEMElementFactoryFctPt bem_element_factory_fct_pt;

  };


}

#endif
