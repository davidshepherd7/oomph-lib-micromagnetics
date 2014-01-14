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

namespace oomph
{

  // Function types
  // ============================================================

  /// \short Function pointer type for function which returns a BEM
  /// element.
  typedef MicromagBEMElementEquations*
  (*BEMElementFactoryFctPt)(FiniteElement* const, const int&);

  /// Function pointer type for function to create flux meshes
  typedef Mesh* (*FluxMeshFactoryFctPt)(Mesh* bulk_mesh_pt,
                                        const Vector<unsigned> &boundaries);



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
      Renormalise_each_time_step(false),
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
      Disable_ms = false;
      Phi_problem_pt = 0;
      Phi_1_problem_pt = 0;

      // Debugging switches
      Pin_boundary_m = false;
      Use_fd_jacobian = false;
      Use_hlib = false;
    }

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
                ->node_to_global(nd);
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

      oomph_info << std::endl
                << "Newton step " << Nnewton_iter_taken + 1 << std::endl
                << "---------------------------------------" << std::endl;
    }


    virtual void actions_before_newton_solve()
    {
      // Call base class version
      MyProblem::actions_before_newton_solve();

      // Update BEM magnetostatics boundary conditions (if we are doing
      // them fully implicitly).
      maybe_update_bem_boundary_conditions();
    }

    virtual void actions_after_explicit_timestep()
      {
        MyProblem::actions_after_explicit_timestep();

        // We need to keep M normalised...
        oomph_info << "Renormalising nodal magnetisations." << std::endl;
        renormalise_magnetisation();
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

#ifdef PARANOID

      // From the nmag user manual:
      // [Begin quote M Donahue]
      // * if the spin angle is approaching 180 degrees, then the results are completely bogus.
      // * over 90 degrees the results are highly questionable.
      // * Under 30 degrees the results are probably reliable.
      // [end quote]
      // (the spin angle is the angle between two neighbouring magnetisations).

      Vector<double> a = elemental_max_m_angle_variations();
      double max_angle_var = *std::max_element(a.begin(), a.end());
      if(max_angle_var > MathematicalConstants::Pi/4)
        {
          std::string error_msg
            = "Large angle variations of " + to_string(max_angle_var)
            + " > " + to_string(MathematicalConstants::Pi/4)
            + " across a single element,\n";
          error_msg += "this often means that your mesh is not sufficiently refined.";
          OomphLibWarning(error_msg, OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
        }
#endif

      // Calculate and store the energy (ready to be output)
      calculate_energies();
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

    /// Integrate a function given by func_pt over every element in a mesh
    /// and return the total. This should probably be in the mesh class but
    /// that's core oomph-lib so I'll leave it here.
    double integrate_over_mesh(const ElementalFunction* func_pt,
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
    double integrate_over_problem(const ElementalFunction* func_pt) const
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

    /// \short Calculate energies and store them for easy reference
    /// (e.g. for output).
    void calculate_energies(bool calculate_effective_damping=true)
      {
        // If you want to turn off energy calculations (e.g. for speed)
        // this is the place to do it. Replace values with
        // MyProblem::Dummy_doc_data.

        // Calculate and store new values
        Exchange_energy = exchange_energy();
        Zeeman_energy = zeeman_energy();
        Crystalline_anisotropy_energy = crystalline_anisotropy_energy();
        Magnetostatic_energy = magnetostatic_energy();

        // Store energy for damping calculations
        Previous_energies.push_front(micromagnetic_energy());

        // Keep the list of previous energies reasonably small (we only
        // need N for any bdf<N> calculation).
        if(Previous_energies.size() > 5) Previous_energies.pop_back();

        // Calculate and store effective damping if not disabled.
        if(calculate_effective_damping)
          {
            Effective_damping_constant = effective_damping_used();
            Alt_eff_damp = alt_effective_damping_used();
          }
      }


    /// \short Calculate the total (micromagnetic) energy for all meshes in
    /// the problem.
    double micromagnetic_energy() const
    {
      return Exchange_energy + Zeeman_energy +
        Crystalline_anisotropy_energy + Magnetostatic_energy;
    }


    double exchange_energy() const
    {
      ExchangeEnergyFunction f;
      return integrate_over_problem(&f);
    }


    double zeeman_energy() const
    {
      ZeemanEnergyFunction f;
      return integrate_over_problem(&f);
    }

    double crystalline_anisotropy_energy() const
    {
      CrystallineAnisotropyEnergyFunction f;
      return integrate_over_problem(&f);
    }


    double magnetostatic_energy() const
    {
      MagnetostaticEnergyFunction f;
      return integrate_over_problem(&f);
    }

    double integral_of_dmdt_squared() const
    {
      DmdtSquaredFunction f;
      return integrate_over_problem(&f);
    }

    double dEnergydt() const
    {
      dExchangeEnergydtFunction de_exdt;
      double I_de_exdt = integrate_over_problem(&de_exdt);

      dZeemanEnergydtFunction de_zeedt;
      double I_de_zeedt = integrate_over_problem(&de_zeedt);

      dCrystallineAnisotropydtEnergyFunction de_cadt;
      double I_de_cadt = integrate_over_problem(&de_cadt);

      dMagnetostaticEnergydtFunction de_ms;
      double I_de_ms = integrate_over_problem(&de_ms);

      return I_de_exdt + I_de_zeedt + I_de_cadt + I_de_ms;
    }

    double alt_dEnergydt() const
    {
      // Make a BDF2 time stepper to look up weights from (because I'm
      // lazy...)
      BDF<2> bdf;
      TimeStepper* node_ts_pt = mesh_pt()->finite_element_pt(0)->node_pt(0)
        ->time_stepper_pt();
      bdf.time_pt() = node_ts_pt->time_pt();
      bdf.set_weights();

      // Calculate first derivative
      double deriv = 0.0;
      for(unsigned t=0;t<bdf.ntstorage();t++)
        {
          deriv += bdf.weight(1,t) * Previous_energies[t];
        }

      return deriv;
    }

    /// \short Compute the effective damping constant (alpha) for the
    /// previous time step (see Albuquerque2001).
    double alt_effective_damping_used() const
      {
        // Integral over all space of (dm/dt)^2 used in last step
        double dmdt_squared = integral_of_dmdt_squared(); //??ds

        // If no change then damping is undefined
        if(dmdt_squared  == 0) return nan("");

        // Forumla from Albuquerque2001 & dAquino2005
        double dEdt = alt_dEnergydt();
        double effective_alpha = - dEdt / dmdt_squared;

        return effective_alpha;
      }


    /// \short Compute the effective damping constant (alpha) for the
    /// previous time step (see Albuquerque2001).
    double effective_damping_used() const
      {
        // Integral over all space of (dm/dt)^2 used in last step
        double dmdt_squared = integral_of_dmdt_squared();

        // If no change then damping is undefined
        if(dmdt_squared  == 0) return nan("");

        // Forumla from Albuquerque2001 & dAquino2005
        double dEdt = dEnergydt();
        double effective_alpha = - dEdt / dmdt_squared;

        return effective_alpha;
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
        << Trace_seperator << "alt_effective_damping";

    }

    void write_additional_trace_data(std::ofstream& trace_file) const
    {

      // Get average (and standard deviation) of |m| - 1
      double m_error_avg(0), m_error_stddev(0);
      norm_m_error(m_error_avg, m_error_stddev);

      Vector<double> angle_variations = elemental_max_m_angle_variations();
      Vector<double> mean_m = mean_magnetisation();

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
        << Trace_seperator << Alt_eff_damp;

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
          }

        calculate_energies(false);
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
          return -1;
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
              SemiImplicitMicromagEquations* ele_pt
                = checked_dynamic_cast<SemiImplicitMicromagEquations*>
                (mesh_pt()->element_pt(e));

              // Interpolate
              Vector<double> ms;
              MMInterpolator intp(ele_pt, s);
              ele_pt->get_magnetostatic_field(&intp, ms);

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

    /// \short Non-const access function for Renormalise_each_time_step.
    bool& renormalise_each_time_step() {return Renormalise_each_time_step;}

    /// \short Const access function for Renormalise_each_time_step.
    bool renormalise_each_time_step() const {return Renormalise_each_time_step;}

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
        return (!Decoupled_ms) && (!Disable_ms);
      }

    /// Can we check the solution using Mallinson's exact time + phi
    /// solutions?
    bool Compare_with_mallinson;

    bool Disable_ms;

    bool Pin_boundary_m;
    bool Use_fd_jacobian;
    bool Use_hlib;

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

    /// Normalise magnetisation problem after each step?
    bool Renormalise_each_time_step;

    /// \short Exchange_energy, computed after previous Newton solve.
    double Exchange_energy;

    /// \short Zeeman energy, computed after previous Newton solve.
    double Zeeman_energy;

    /// \short Crystalline anisotropy energy, computed after previous Newton
    /// solve.
    double Crystalline_anisotropy_energy;

    /// \short Magnetostatic energy, computed after previous Newton solve.
    double Magnetostatic_energy;

    /// \short Pointer to class for handling BEM
    BoundaryElementHandler* Bem_handler_pt;

    /// \short Mesh for flux elements to impose boundary condition on phi1.
    Mesh* Flux_mesh_pt;

public:

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


  public:

    /// Are we solving for ms "properly" or using a separate solve?
    bool Decoupled_ms;

    BEMElementFactoryFctPt Bem_element_factory_pt;

    GenericPoissonProblem::FluxMeshFactoryFctPt Phi_1_flux_mesh_factory_fct_pt;

    void build_decoupled_ms(Vector<Mesh*>& llg_mesh_pts,
                            Vector<Mesh*>& phi_mesh_pts,
                            Vector<Mesh*>& phi_1_mesh_pts);

    /// \short Solve for the magnetostatic field.
    void magnetostatics_solve()
    {
      if(Decoupled_ms)
        {
          oomph_info << std::endl
                    << "BEM solve" << std::endl
                    << "--------------------------" <<std::endl;

          // solve for phi1
          oomph_info << "solving phi1" << std::endl;
          phi_1_problem_pt()->newton_solve();

          // update boundary values of phi
          oomph_info << "solving BEM" << std::endl;
          double t_start = TimingHelpers::timer();
          Bem_handler_pt->get_bem_values(Phi_boundary_values_pts);
          double t_end = TimingHelpers::timer();
          oomph_info << "BEM time taken: " << t_end - t_start << std::endl;

          // push old phi values back in time (so that we can use them later to
          // get time derivatives of the field). Note that we don't use the
          // problem's shift time values function because we don't want to
          // shift the timestepper (that has been done by the llg problem
          // already) and we don't have any external data to shift.
          phi_problem_pt()->mesh_pt()->shift_time_values();

          // solve for phi
          oomph_info << "solving phi" << std::endl;
          phi_problem_pt()->newton_solve();

          oomph_info << "mean field is " << average_magnetostatic_field() << std::endl;
        }
      else
        {
          std::string err = "requested magnetostatics solve but problem is not ";
          err += "decoupled, so we can't do it, doing nothing instead.";
          throw OomphLibWarning(err, OOMPH_EXCEPTION_LOCATION,
                                OOMPH_CURRENT_FUNCTION);
        }
    }

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
    MMArgs() : h_app_fct_pt(0),
               magnetic_parameters_pt(0) {}

    virtual void set_flags()
    {
      MyCliArgs::set_flags();

      specify_command_line_flag("-initm", &initial_m_name);
      initial_m_name = "z";

      specify_command_line_flag("-happ", &h_app_name);
      h_app_name = "minus_z";

      specify_command_line_flag("-mag-params", &magnetic_parameters_name);
      magnetic_parameters_name = "simple-llg";

      specify_command_line_flag("-renorm_m", &Renormalise);
      Renormalise = -1;

      specify_command_line_flag("-dampc", &dampc);
      dampc = -10;

      // Flags automatically default to false
      specify_command_line_flag("-numerical-BEM");
      specify_command_line_flag("-decoupled-ms");
      specify_command_line_flag("-disable-ms");
      specify_command_line_flag("-pin-boundary-m");
    }


    virtual void run_factories()
    {
      using namespace SemiImplicitFactories;
      using namespace LLGFactories;
      using namespace Factories;

      decoupled_ms = command_line_flag_has_been_set("-decoupled-ms");

      // Figure out how to build meshes
      if(decoupled_ms)
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
      magnetic_parameters_name = to_lower(magnetic_parameters_name);

      initial_condition_fpt = InitialM::initial_m_factory(initial_m_name);
      h_app_fct_pt = HApp::h_app_factory(h_app_name);
      magnetic_parameters_pt =
        magnetic_parameters_factory(magnetic_parameters_name);

      if(command_line_flag_has_been_set("-dampc"))
        {
          magnetic_parameters_pt->gilbert_damping() = dampc;
        }

      // Copy flags into bools in this class
      use_numerical_integration_bem = command_line_flag_has_been_set("-numerical-BEM");
      disable_ms = command_line_flag_has_been_set("-disable-ms");
      pin_boundary_m = command_line_flag_has_been_set("-pin-boundary-m");

      // Only want one of these to be true at once ??ds enumeration instead?
      if(decoupled_ms && disable_ms)
        {
          std::string err = "Requested decoupled ms and disabled ms, ";
          err += "I'm just going to disabled it.";
            throw OomphLibWarning(err, OOMPH_EXCEPTION_LOCATION,
                                  OOMPH_CURRENT_FUNCTION);
          decoupled_ms = false;
        }


      if(decoupled_ms)
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

        if(decoupled_ms)
          {
            // Also build separate poisson meshes if needed
            using namespace SemiImplicitFactories;
            phi_mesh_pts = build_meshes_helper(phi_mesh_factory);
            phi_1_mesh_pts = build_meshes_helper(phi_mesh_factory);
          }
      }

    virtual void assign_specific_parameters(MyProblem* problem_pt) const
    {
      LLGProblem* llg_pt = checked_dynamic_cast<LLGProblem*>(problem_pt);
      llg_pt->applied_field_fct_pt() = h_app_fct_pt;
      llg_pt->set_mag_parameters_pt(magnetic_parameters_pt);
      llg_pt->renormalise_each_time_step() = renormalise_flag();
      llg_pt->Pin_boundary_m = pin_boundary_m;
      llg_pt->Decoupled_ms = decoupled_ms;
      llg_pt->Disable_ms = disable_ms;

      // ??ds this should maybe be a general one?
      llg_pt->Use_fd_jacobian = use_fd_jacobian; //??ds

      // Set exact solution if we have one
      if((h_app_name == "minus_z")
         && (initial_m_name == "z")
         && (magnetic_parameters_pt->gilbert_damping() != 0.0)
         && disable_ms)
        {
          llg_pt->Compare_with_mallinson = true;
        }

      llg_pt->Phi_1_flux_mesh_factory_fct_pt = phi_1_flux_mesh_factory_fct_pt;

      llg_pt->Bem_element_factory_pt = bem_element_factory_fct_pt;

    }

    /// Write out all args (in a parseable format) to a stream.
    virtual void dump_args(std::ostream& out_stream) const
    {
      MyCliArgs::dump_args(out_stream);

      out_stream
        << "initial_m " << initial_m_name << std::endl
        << "h_app " << h_app_name << std::endl
        << "mag_params " << magnetic_parameters_name << std::endl
        << "Renormalise " << Renormalise << std::endl
        << "damping_parameter_override " << dampc << std::endl
        << "numerical-BEM " << use_numerical_integration_bem << std::endl
        << "decoupled_ms " << decoupled_ms << std::endl
        << "disable_ms " << disable_ms << std::endl
        << "pin_boundary_m " << pin_boundary_m << std::endl
        ;
    }


    bool renormalise_flag() const
    {
      // If flag not set then only do it for bdf timesteppers
      if(Renormalise == -1)
        {
          return (time_stepper_name == "bdf2")
            || (time_stepper_name == "bdf1");
        }
      // Otherwise do what the flag says
      else
        {
          return bool(Renormalise);
        }
    }

    HApp::HAppFctPt h_app_fct_pt;
    MagneticParameters* magnetic_parameters_pt;

    // Strings for input to factory functions
    std::string initial_m_name;
    std::string h_app_name;
    std::string magnetic_parameters_name;


    /// Flag to control renormalisation of |m| after each step. -1 =
    /// default for timestepper, 0 = off, 1 = on.
    int Renormalise;

    double dampc;

    bool use_numerical_integration_bem;


    bool decoupled_ms;
    bool disable_ms;
    bool pin_boundary_m;

    Vector<Mesh*> phi_1_mesh_pts;
    Vector<Mesh*> phi_mesh_pts;

    GenericPoissonProblem::FluxMeshFactoryFctPt phi_1_flux_mesh_factory_fct_pt;
    BEMElementFactoryFctPt bem_element_factory_fct_pt;

  };


}

#endif
