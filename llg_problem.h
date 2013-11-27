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


namespace oomph
{

  class ResidualCalculator;

  // ============================================================
  ///
  // ============================================================
  class LLGProblem : public MyProblem
  {
  public:


    /// Function type which creates flux meshes
    typedef Mesh* (*FluxMeshFactoryFctPt)(Mesh* bulk_mesh_pt,
                                          const Vector<unsigned> &boundaries);

    /// Default constructor - do nothing except nulling pointers.
    LLGProblem() :
      Compare_with_mallinson(false),
      Swap_solver_large_dt(false),
      Use_fd_jacobian(false),
      Use_implicit_ms(false),
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

      // Bem stuff
      Bem_handler_pt = 0;
      Flux_mesh_pt = 0;
      Flux_mesh_factory_pt = 0;
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

    void actions_before_newton_step()
    {
      std::cout << std::endl
                << "Newton step " << Nnewton_iter_taken + 1 << std::endl
                << "---------------------------------------" << std::endl;

      if(Swap_solver_large_dt)
        {
          if(this->time_pt()->dt() > 1e-2)
            {
              linear_solver_pt() = Super_LU_solver_pt;
            }
          else
            {
              linear_solver_pt() = My_linear_solver_pt;
            }
        }
    }


    virtual void actions_before_newton_solve()
    {
      // Call base class version
      MyProblem::actions_after_newton_solve();

      // Update BEM magnetostatics boundary conditions (if we are doing
      // them fully implicitly).
      maybe_update_bem_boundary_conditions();
    }

    void actions_after_newton_solve()
    {

      std::cout << std::endl
                << "Finalising" << std::endl
                << "-------------------------------------------" << std::endl;

      // If we're using BDF we need to keep M normalised.
      if(renormalise_each_time_step())
        {
          std::cout << "Renormalising nodal magnetisations." << std::endl;
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

    virtual void actions_after_newton_step()
    {
      // Call base class actions function
      MyProblem::actions_after_newton_step();

      // Update BEM magnetostatics boundary conditions (if we are doing them
      // fully implicitly).
      maybe_update_bem_boundary_conditions();
    }

    void maybe_update_bem_boundary_conditions()
      {
        if(Use_implicit_ms)
          {
            throw OomphLibError("Not yet implemented.",
                                OOMPH_EXCEPTION_LOCATION,
                                OOMPH_CURRENT_FUNCTION);

            Bem_handler_pt->get_bem_values_and_copy_into_values(phi_index());
          }
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

        // Keep the lsit of previous energies reasonably small (we only
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

    /// Set up an initial M
    void set_initial_condition(const InitialM::InitialMFctPt initial_m_pt);


    /// \short Return a vector of the maximum angle variation in each
    /// element.
    Vector<double> elemental_max_m_angle_variations() const
    {
      Vector<double> angles;

      for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
        {
          // Skip non-bulk meshes
          if(mesh_pt(msh)->node_pt(0)->ndim() != Dim) continue;

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

    // /// Elementwise calculation of m . dm/dn on the boundaries. This gives a good measure of
    // /// orthogonality of m and dmdn.
    // double mean_orthogonality_m_error() const
    // {
    //   // unsigned nele = surface_exchange_mesh_pt()->nelement();
    //   // double temp_sum = 0.0;

    //   // // Calculate m . dm/dn for each element
    //   // for(unsigned e=0; e < nele; e++)
    //   //   {
    //   //     MicromagFluxElement<MicromagEquations>* ele_pt
    //   //       = checked_dynamic_cast<MicromagFluxElement<MicromagEquations>*>
    //   //       (surface_exchange_mesh_pt()->element_pt(e));
    //   //     Vector<double> s(ele_pt->dim(), 0.5); // middle of element
    //   //     temp_sum += ele_pt->interpolated_mdotdmdn_micromag(s);
    //   //   }

    //   // // Take average
    //   // return temp_sum / double(nele);

    //   std::ostringstream error_msg;
    //   error_msg << "Not implemented.";
    //   throw OomphLibError(error_msg.str(),
    //                       OOMPH_CURRENT_FUNCTION,
    //                       OOMPH_EXCEPTION_LOCATION);

    // }

    // void orthogonality_m_error(double &orthogonality_error_avg,
    //                            double &orthogonality_error_stddev) const
    // {
    //   orthogonality_error_avg = mean_orthogonality_m_error();

    //   unsigned nele = surface_exchange_mesh_pt()->nelement();
    //   double temp_sum = 0.0;

    //   // Calculate m . dm/dn for each element
    //   for(unsigned e=0; e < nele; e++)
    //     {
    //       MicromagFluxElement<MicromagEquations>* ele_pt
    //         = checked_dynamic_cast<MicromagFluxElement<MicromagEquations>*>
    //         (surface_exchange_mesh_pt()->element_pt(e));
    //       Vector<double> s(ele_pt->dim(), 0.5); // middle of element
    //       temp_sum += pow( ele_pt->interpolated_mdotdmdn_micromag(s)
    //                        - orthogonality_error_avg, 2);
    //     }

    //   // Take stddev
    //   orthogonality_error_stddev = std::sqrt(temp_sum / double(nele));
    // }


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

    // Access functions
    // ============================================================

    unsigned m_index(const unsigned &j) const
    {return this->ele_pt()->m_index_micromag(j);}

    unsigned phi_index() const {return this->ele_pt()->phi_index_micromag();}

    unsigned phi_1_index() const
    {return this->ele_pt()->phi_1_index_micromag();}

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

    //     /// \short Non-const access function for Surface_exchange_mesh_pt.
    //     Mesh*& surface_exchange_mesh_pt() {return Surface_exchange_mesh_pt;}

    //     /// \short Const access function for Surface_exchange_mesh_pt.
    //     Mesh* surface_exchange_mesh_pt() const
    //     {
    // #ifdef PARANOID
    //       if(Surface_exchange_mesh_pt == 0)
    //         {
    //           std::ostringstream error_msg;
    //           error_msg << "Surface exchange mesh pointer not set.";
    //           throw OomphLibError(error_msg.str(),
    //                               OOMPH_CURRENT_FUNCTION,
    //                               OOMPH_EXCEPTION_LOCATION);
    //         }
    // #endif

    //       return Surface_exchange_mesh_pt;
    //     }

    /// \short Get a pointer to a MicromagEquations element.
    MicromagEquations* ele_pt() const
    {
      return checked_dynamic_cast<MicromagEquations*>
        (mesh_pt(0)->element_pt(0));
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

    /// Can we check the solution using Mallinson's exact time + phi
    /// solutions?
    bool Compare_with_mallinson;

    /// \short Should we swap to superlu for large dt solves?
    bool Swap_solver_large_dt;

    bool Use_fd_jacobian;
    bool Use_implicit_ms;

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

    /// \short Pointer to function for creating the flux mesh (can't be
    /// hard coded becuase it depends on the element type, which depends on
    /// dimension etc.)
    FluxMeshFactoryFctPt Flux_mesh_factory_pt;

public:
    /// \short Recomputed effective damping constant for the last time step
    /// (based on actual change in energy).
    double Effective_damping_constant;
    double Alt_eff_damp;
    std::deque<double> Previous_energies;

  private:

    /// \short ??ds
    void create_surface_exchange_elements(const unsigned& b);

    /// Inaccessible copy constructor
    LLGProblem(const LLGProblem & dummy)
    {BrokenCopy::broken_copy("LLGProblem");}

    /// Inaccessible assignment operator
    void operator=(const LLGProblem &dummy)
    {BrokenCopy::broken_assign("LLGProblem");}

  };


  /// \short Functions for constructing things needed for LLG problems
  namespace LLGFactories
  {
    /// \short Make a mesh as specified by an input argument. Refined
    /// according to the given refinement level (in some way appropriate
    /// for that mesh type). Assumption: this will be passed into a
    /// problem, which will delete the pointer when it's done.
    Vector<Mesh*> mesh_factory(const std::string& _mesh_name,
                               int refinement_level,
                               TimeStepper* time_stepper_pt,
                               unsigned nnode1d = 2);

    LLGResidualCalculator* residual_calculator_factory(const std::string& residual);

    /// \short Create a variable order quadrature object based on the
    /// dimension and shape of the element. Only works for some element
    /// types.
    Integral* variable_order_integrator_factory(const FiniteElement* const el_pt);


    /// \short Function pointer type for function which returns a BEM
    /// element.
    typedef MicromagBEMElementEquations*
    (*BEMElementFactoryFctPt)(FiniteElement* const, const int&);

    typedef Mesh* (*FluxMeshFactoryFctPt)(Mesh* bulk_mesh_pt,
                                          const Vector<unsigned> &boundaries);


    /// \short very simple function: create a new face element of type
    /// ELEMENT.
    template<class ELEMENT>
    MicromagBEMElementEquations* bem_element_factory(FiniteElement* ele,
                                                     const int& face)
    {
      return new ELEMENT(ele, face);
    }


    /// \short Create a function to create bem elements based on the
    /// elements used in the bulk mesh.
    BEMElementFactoryFctPt bem_element_factory_factory
    (const FiniteElement* bulk_ele_pt);

    FluxMeshFactoryFctPt
    mm_flux_mesh_factory_factory(const FiniteElement* bulk_ele_pt);
  }


  /// Base command line args processing class for pure llg and semi
  /// implicit (and any other magnetism problems).
  class MMArgs : public MyCliArgs
  {
  public:
    /// Constructor: Initialise pointers to null.
    MMArgs() : initial_m_fct_pt(0), h_app_fct_pt(0),
               magnetic_parameters_pt(0),
               residual_calculator_pt(0) {}

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

      specify_command_line_flag("-resi", &residual_to_use);
      residual_to_use = "llg";
    }


    virtual void run_factories()
    {
      MyCliArgs::run_factories();

      initial_m_name = to_lower(initial_m_name);
      h_app_name = to_lower(h_app_name);
      magnetic_parameters_name = to_lower(magnetic_parameters_name);

      initial_m_fct_pt = InitialM::initial_m_factory(initial_m_name);
      h_app_fct_pt = HApp::h_app_factory(h_app_name);
      magnetic_parameters_pt =
        Factories::magnetic_parameters_factory(magnetic_parameters_name);

      if(command_line_flag_has_been_set("-dampc"))
        {
          magnetic_parameters_pt->gilbert_damping() = dampc;
        }

      residual_calculator_pt = LLGFactories::
        residual_calculator_factory(residual_to_use);

      // If we are using llg residual with midpoint method then we need to
      // swap residuals over in the explicit predictor time steps. Put in
      // the class to do this here.
      if(residual_to_use == "llg" && time_stepper_name == "midpoint")
        {
          MidpointMethod* mp_pt = checked_dynamic_cast<MidpointMethod*>
            (time_stepper_pt);

          // We've already run the base classes factories so we have the
          // timestepper ready to play with. Create and set up our rs
          // timestepper.
          ResidualSwappingExplicitTimestepper* rsts_pt
            = new ResidualSwappingExplicitTimestepper;
          rsts_pt->underlying_time_stepper_pt = mp_pt->predictor_pt();
          rsts_pt->residual_pt = checked_dynamic_cast<LLGResidualCalculator*>
            (residual_calculator_pt);

          mp_pt->set_predictor_pt(rsts_pt);
        }
    }

    /// Write out all args (in a parseable format) to a stream.
    virtual void dump_args(std::ostream& out_stream) const
    {
      MyCliArgs::dump_args(out_stream);

      out_stream
        << "initial_m " << initial_m_name << std::endl
        << "h_app " << h_app_name << std::endl
        << "mag_params " << magnetic_parameters_name << std::endl
        << "damping_parameter_override " << dampc << std::endl;
    }


    bool renormalise_flag()
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

    InitialM::InitialMFctPt initial_m_fct_pt;
    HApp::HAppFctPt h_app_fct_pt;
    MagneticParameters* magnetic_parameters_pt;
    LLGResidualCalculator* residual_calculator_pt;

    // Strings for input to factory functions
    std::string initial_m_name;
    std::string h_app_name;
    std::string magnetic_parameters_name;
    std::string residual_to_use;


    /// Flag to control renormalisation of |m| after each step. -1 =
    /// default for timestepper, 0 = off, 1 = on.
    int Renormalise;

    double dampc;

  };


  /// Command line args class for llg problems. Just add the mesh
  /// stuff. ??ds refactor to combine with MMArgs?
  class LLGArgs : public MMArgs
  {
  public:

    /// Constructor
    LLGArgs() {}

    virtual void set_flags()
    {
      MMArgs::set_flags();

      specify_command_line_flag("-mesh", &mesh_name);
      mesh_name = "sq_square";

      specify_command_line_flag("-nnode1d", &nnode1d);
      nnode1d = 2;

      specify_command_line_flag("-implicit-ms");
    }


    virtual void run_factories()
    {
      MMArgs::run_factories();

      // Do the mesh last of all because it can be slow
      mesh_name = to_lower(mesh_name);
      mesh_pts = LLGFactories::mesh_factory(mesh_name, refinement,
                                            time_stepper_pt, nnode1d);

      use_implicit_ms = command_line_flag_has_been_set("-implicit-ms");
    }

    /// Write out all args (in a parseable format) to a stream.
    virtual void dump_args(std::ostream& out_stream) const
    {
      MMArgs::dump_args(out_stream);
      out_stream << "mesh " << mesh_name << std::endl;
      out_stream << "nnode1d " << nnode1d << std::endl;
      out_stream << "use_implicit_ms " << use_implicit_ms << std::endl;
    }


    Vector<Mesh*> mesh_pts;

    unsigned nnode1d;

    bool use_implicit_ms;

    // Strings for input to factory functions
    std::string mesh_name;
  };

}

#endif
