#ifndef OOMPH_ODE_PROBLEM_H
#define OOMPH_ODE_PROBLEM_H

#include "my_cli.h"
#include "my_generic_problem.h"
#include "micromag_types.h"
#include "mallinson_solution.h"

// For creating happ function for llg ode
#include "llg_factories.h"

namespace oomph
{
  using namespace MathematicalConstants;
  using namespace StringConversion;
  using namespace VectorOps;


  typedef TimeSpaceToDoubleVectFctPt TimeValueToDoubleVectFctPt;

  namespace deriv_functions
  {
    class LLODESolution;

    inline Vector<double> cos(const double& time, const Vector<double>& x)
    {
      Vector<double> values(1);
      values[0] = std::cos(time);
      return values;
    }
    inline Vector<double> dcos(const double& t, const Vector<double>& x,
                               const Vector<double>& u)
    {
      Vector<double> deriv(1, 0.0);
      deriv[0] = -1*std::sin(t);
      return deriv;
    }


    inline Vector<double> sin(const double& time, const Vector<double>& x)
    {
      Vector<double> values(1);
      values[0] = std::sin(time);
      return values;
    }
    inline Vector<double> dsin(const double& t, const Vector<double>& x,
                               const Vector<double>& u)
    {
      Vector<double> deriv(1, 0.0);
      deriv[0] = std::cos(t);
      return deriv;
    }


    inline Vector<double> tanh(const double& time, const Vector<double>& x)
    {
      double a = 10;
      double t_step = 1;
      Vector<double> values(1);
      values[0] = std::tanh(a * (time - t_step));
      return values;
    }
    inline Vector<double> dtanh(const double& t, const Vector<double>& x,
                                const Vector<double>& u)
    {
      double a = 10;
      double t_step = 1;
      Vector<double> deriv(1, 0.0);
      double th = std::tanh(a * (t - t_step));
      deriv[0] = a * (-th * th  + 1);
      return deriv;
    }


    inline Vector<double> exp(const double& time, const Vector<double>& x)
    {
      Vector<double> values(1);
      values[0] = std::exp(time);
      return values;
    }
    inline Vector<double> dexp(const double& t, const Vector<double>& x,
                               const Vector<double>& u)
    {
      Vector<double> deriv(1, 0.0);
      deriv[0] = u[0];
      return deriv;
    }


    // A polynomial of degree 2
    inline Vector<double> poly2(const double& time, const Vector<double>& x)
    {
      double b0 = 0.5, b1 = 0, b2 = 1;

      Vector<double> values(1);
      values[0] =  b2*time*time + b1*time +b0;
      return values;
    }
    inline Vector<double> dpoly2(const double& t, const Vector<double>& x,
                                 const Vector<double>& u)
    {
      double b1 = 0, b2 = 1;

      Vector<double> deriv(1, 0.0);
      deriv[0] = 2*t*b2 + b1;
      return deriv;
    }


    // A polynomial of degree 3
    inline Vector<double> poly3(const double& time, const Vector<double>& x)
    {
      double a0 = 0.5, a1 = 0, a2 = 0, a3 = 1;

      Vector<double> values(1);
      values[0] = a3*time*time*time + a2*time*time + a1*time +a0;
      return values;
    }
    inline Vector<double> dpoly3(const double& t, const Vector<double>& x,
                                 const Vector<double>& u)
    {
      double a1 = 0, a2 = 0, a3 = 1;

      Vector<double> deriv(1, 0.0);
      deriv[0] = 3*t*t*a3 + 2*t*a2 + a1;
      return deriv;
    }

    // stiff ode, example from Iserles pg. 54
    inline Vector<double> stiff_test(const double& time, const Vector<double>& x)
    {
      Vector<double> x1(2), x2(2);
      x1[0] = 0; x1[1] = 0;
      x2[0] = 1.0; x2[1] = 999.0/10;

      Vector<double> values(2);
      values[0] = std::exp(-100*time)*x1[0] + std::exp(-0.1*time)*x2[0];
      values[1] = std::exp(-100*time)*x1[1] + std::exp(-0.1*time)*x2[1];
      return values;
    }
    inline Vector<double> dstiff_test(const double& t, const Vector<double>& x,
                                      const Vector<double>& u)
    {
      Vector<double> deriv(2, 0.0);
      deriv[0] = -100*u[0] + u[1];
      deriv[1] = 0*u[0] - 0.1*u[1];
      return deriv;
    }

    /// A damped oscillation solution
    class DampedOscillation : public SolutionFunctorBase
    {
    public:
      /// Constructor
      DampedOscillation()
      {
        Beta = 0.1;
        Omega = 2;
      }

      /// Virtual destructor
      virtual ~DampedOscillation() {}

      /// Function call
      Vector<double> operator()(const double& t, const Vector<double>& x) const
      {
        Vector<double> values(1);
        values[0] = std::exp(-Beta*t) * std::sin(Omega*t);
        return values;
      }

      /// Derivative call
      Vector<double> derivative(const double& t, const Vector<double>& x,
                                const Vector<double>& u) const
      {
        Vector<double> deriv(1, 0.0);
        deriv[0] = -Beta * std::exp(-Beta*t) * std::sin(Omega*t)
          + Omega * std::exp(-Beta*t) * std::cos(Omega*t);
        return deriv;
      }

      double Beta;
      double Omega;
    };

    /// Another stiff solution: Atkinson equation (8.1) pg 128
    class SimpleStiffTest : public SolutionFunctorBase
    {
    public:
      /// Constructor
      SimpleStiffTest()
      {
        Lambda = 100;
        Y_intial = 1;
      }

      /// Virtual destructor
      virtual ~SimpleStiffTest() {}

      /// Function call
      Vector<double> operator()(const double& t, const Vector<double>& x) const
      {
        Vector<double> values(1);
        values[0] = std::exp(-Lambda * t) * Y_intial;
        return values;
      }

      /// Derivative call
      Vector<double> derivative(const double& t, const Vector<double>& x,
                                const Vector<double>& u) const
      {
        Vector<double> deriv(1, 0.0);
        deriv[0] = -Lambda * u[0];
        return deriv;
      }

      double Lambda;
      double Y_intial;
    };

    /// Another stiff solution: Atkinson pg. 158, also example 8.2 pg 129.
    class OrderReductionTest : public SolutionFunctorBase
    {
    public:
      /// Constructor
      OrderReductionTest()
      {
        Lambda = -100;
      }

      /// Virtual destructor
      virtual ~OrderReductionTest() {}

      /// Function call
      Vector<double> operator()(const double& t, const Vector<double>& x) const
      {
        Vector<double> values(1);
        values[0] = std::sin(t);
        return values;
      }

      /// Derivative call
      Vector<double> derivative(const double& t, const Vector<double>& x,
                                const Vector<double>& u) const
      {
        Vector<double> deriv(1, 0.0);
        deriv[0] = Lambda*u[0] - Lambda*std::sin(t) + std::cos(t);
        return deriv;
      }

      double Lambda;
    };

  }

  namespace ODEFactories
  {
    // Pick an exact solution using a name
    SolutionFunctorBase* exact_solutions_factory(const std::string& exact_name,
                                                 const std::string& initial_m_name,
                                                 double sc_beta=0.0);
  }



  //=====================================================================
  /// Element that integrates single ode with bdf scheme
  //=====================================================================
  class ODEElement : public GeneralisedElement
  {

  public:


    /// Constructor: Pass timestepper
    ODEElement(TimeStepper* timestepper_pt,
               SolutionFunctorBase* exact_solution_pt)
    {
      Exact_solution_pt = exact_solution_pt;

      Vector<double> exact = this->exact_solution(0);
      unsigned nvalue = exact.size();

      add_internal_data(new Data(timestepper_pt, nvalue));

      Use_fd_jacobian = false;
    }

    virtual ~ODEElement() {}

    unsigned nvalue() const
    {
      return internal_data_pt(0)->nvalue();
    }

    /// Get residuals
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Get pointer to one-and-only internal data object
      Data* dat_pt = internal_data_pt(0);

      // Get it's values
      Vector<double> u(nvalue(), 0.0);
      dat_pt->value(u);

      // Get timestepper
      TimeStepper* timestepper_pt = dat_pt->time_stepper_pt();

      // Get continuous time
      double t = timestepper_pt->time();

      Vector<double> deriv = derivative_function(t, u);
      for(unsigned j=0, nj=deriv.size(); j<nj; j++)
        {
          // Get dudt approximation from timestepper
          double dudt = timestepper_pt->time_derivative(1, dat_pt, j);

          // Residual is difference between the exact derivative and our
          // timestepper's derivative estimate.
          residuals[j] = deriv[j] - dudt;
        }
    }

    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Get residuals
      fill_in_contribution_to_residuals(residuals);

      if(Exact_solution_pt->have_jacobian() && !Use_fd_jacobian)
        {
          // get df/du jacobian
          double t = internal_data_pt(0)->time_stepper_pt()->time();
          Vector<double> dummy, u(nvalue(), 0.0);
          internal_data_pt(0)->value(u);
          Exact_solution_pt->jacobian(t, dummy, u, jacobian);

          // We need jacobian of residual = f - dudt so subtract diagonal
          // (dudt)/du term.
          const double a = internal_data_pt(0)->time_stepper_pt()->weight(1,0);
          const unsigned n = nvalue();
          for(unsigned i=0; i<n; i++)
            {
              jacobian(i, i) -= a;
            }
        }
      else
        {
          // Use FD for jacobian
          GeneralisedElement::fill_in_jacobian_from_internal_by_fd
            (residuals, jacobian, true);
        }
    }

    void fill_in_contribution_to_mass_matrix(Vector<double>& residuals,
                                             DenseMatrix<double>& mm)
    {
      fill_in_contribution_to_residuals(residuals);
      for(unsigned j=0, nj=nvalue(); j<nj; j++)
        {
          mm(j, j) = 1;
        }
    }

    /// Exact solution
    Vector<double> exact_solution(const double& t) const
    {
#ifdef PARANOID
      if(Exact_solution_pt == 0)
        {
          throw OomphLibError("No exact solution function",
                              OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif
      Vector<double> dummy_x;
      return (*Exact_solution_pt)(t, dummy_x);
    }

    /// Exact solution
    Vector<double> derivative_function(const double& t,
                                       const Vector<double>& u)
    {
#ifdef PARANOID
      if(Exact_solution_pt == 0)
        {
          throw OomphLibError("No derivative function",
                              OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif
      Vector<double> dummy_x;
      return Exact_solution_pt->derivative(t, dummy_x, u);
    }

    SolutionFunctorBase* Exact_solution_pt;

    bool Use_fd_jacobian;
  };


  class ODEProblem : public MyProblem
  {
  public:

    /// constructor
    ODEProblem()
    {
      // Don't output to trace file every step, often too many steps
      Always_write_trace = false;

      Use_fd_jacobian = false;
    }

    virtual ~ODEProblem() {}

    virtual void build(Vector<Mesh*>& bulk_mesh_pts)
    {
      // Call the underlying build
      MyProblem::build(bulk_mesh_pts);

      // Set up the global mesh
      build_global_mesh();

      element_pt()->Use_fd_jacobian = Use_fd_jacobian;

      // assign equation numbers
      this->assign_eqn_numbers();
      oomph_info << "Number of equations: " << ndof() << std::endl;
    }

    void set_initial_condition(const InitialConditionFct& ic)
    {
      // Loop over current & previous timesteps
      const unsigned nprev_values = time_stepper_pt()->nprev_values();
      for(unsigned t=0; t<nprev_values+1; t++)
        {
          double time = time_pt()->time(t);

          std::cout << "setting IC at time =" << time << std::endl;

          // Get + set the (only) value
          Vector<double> dummy(nvalue(), 1.0);
          Vector<double> values = ic(time, dummy);

          for(unsigned j=0, nj=nvalue(); j<nj; j++)
            {
              mesh_pt()->element_pt(0)->internal_data_pt(0)
                ->set_value(t, j, values[j]);
            }
        }

      actions_after_set_initial_condition();
    }

    virtual void write_additional_trace_headers(std::ofstream& trace_file) const
    {
      trace_file << Trace_seperator << "exact";
    }

    virtual void write_additional_trace_data(const unsigned& t_hist,
                                             std::ofstream& trace_file) const
    {
      trace_file << Trace_seperator << exact_solution(time_pt()->time(t_hist));
    }

    virtual double get_error_norm(const unsigned& t_hist=0) const
    {
      Vector<double> val = trace_values(t_hist);
      Vector<double> exact = exact_solution(time_pt()->time(t_hist));

      return VectorOps::two_norm_diff(val, exact);
    }

    Vector<double> exact_solution(const double& time) const
    {
      ODEElement* el_pt = checked_dynamic_cast<ODEElement*>
        (mesh_pt()->element_pt(0));
      return el_pt->exact_solution(time);
    }

    /// Error norm: use abs(error in data).
    double global_temporal_error_norm()
    {
      Data* dat_pt=mesh_pt()->element_pt(0)->internal_data_pt(0);

      return std::abs(ts_pt()->temporal_error_in_value(dat_pt, 0));
    }

    Vector<double> trace_values(const unsigned& t_hist=0) const
    {
      return solution(t_hist);
    }

    ODEElement* element_pt()
    {return checked_dynamic_cast<ODEElement*>(mesh_pt()->element_pt(0));}

    const ODEElement* element_pt() const
    {return checked_dynamic_cast<ODEElement*>(mesh_pt()->element_pt(0));}

    TimeStepper* ts_pt() const
    {
      return element_pt()->internal_data_pt(0)->time_stepper_pt();
    }

    unsigned nvalue() const
    {
      return element_pt()->nvalue();
    }

    // Output solution
    void output_solution(const unsigned& t, std::ostream& outstream,
                         const unsigned& npoints=2) const
    {
      std::cout << solution(t) << std::endl;
      outstream << solution(t) << std::endl;

      // npoints is ignored
    }

    Vector<double> solution(const unsigned& timestep=0) const
    {
      Data* dat_pt=mesh_pt()->element_pt(0)->internal_data_pt(0);
      Vector<double> solution(nvalue(), 0.0);
      dat_pt->value(timestep, solution);
      return solution;
    }

    bool Use_fd_jacobian;

  };


  class LLGODEProblem : public ODEProblem
  {

  public:

    LLGODEProblem()
    {
      Exchange_energy = MyProblem::Dummy_doc_data;
      Zeeman_energy = MyProblem::Dummy_doc_data;
      Crystalline_anisotropy_energy = MyProblem::Dummy_doc_data;
      Magnetostatic_energy = MyProblem::Dummy_doc_data;
      Effective_damping_constant = MyProblem::Dummy_doc_data;
      Alt_eff_damp = MyProblem::Dummy_doc_data;

      Mallinson_applicable = false;
      Magnetic_parameters_pt = 0;

      // Default to effectively infinity => never renormalise
      M_length_error_tol = 1e200;

      Renormalise_history = false;
    }

    virtual ~LLGODEProblem()
    {
      delete Magnetic_parameters_pt; Magnetic_parameters_pt = 0;
    }

    bool should_renormalise()
    {
#ifdef PARANOID
      if(M_length_error_tol < 0.0)
	{
	  std::string err = "Error tol is negative";
	  throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);
	}
#endif

      return m_length_error() > M_length_error_tol;
    }

    void renormalise_magnetisation()
    {
      oomph_info << "Renormalising nodal magnetisations." << std::endl;
      unsigned tn=1;
      if(Renormalise_history)
	{
	  tn = time_stepper_pt()->nprev_values();
	}

      for(unsigned ti=0; ti<tn; ti++)
	{
	  Vector<double> values(3, 0.0);
	  mesh_pt()->element_pt(0)->internal_data_pt(0)->value(ti, values);
	  VectorOps::normalise(values);

	  for(unsigned j=0; j<3; j++)
	    {
	      mesh_pt()->element_pt(0)->internal_data_pt(0)
		->set_value(ti, j, values[j]);
	    }
	}
    }

    virtual void actions_after_newton_solve()
    {
      ODEProblem::actions_after_newton_solve();

      // Renormalise if needed
      if(should_renormalise())
	{
	  renormalise_magnetisation();
	}

      // Update energies etc.
      Vector<double> dummy;

      double prev_energy = Zeeman_energy + Crystalline_anisotropy_energy;

      Zeeman_energy = -VectorOps::dot(solution(),
				      Magnetic_parameters_pt->h_app(time(), dummy));
      Crystalline_anisotropy_energy = 0.0; //??ds not implemented
      double energy = Zeeman_energy + Crystalline_anisotropy_energy;

      if(Magnetic_parameters_pt->gilbert_damping() != 0)
	{
	  Alt_eff_damp = effective_damping(energy, prev_energy);
	}
    }

    virtual void actions_after_set_initial_condition()
    {
      ODEProblem::actions_after_set_initial_condition();

      // Get energies
      Vector<double> dummy;
      Zeeman_energy = -VectorOps::dot(solution(),
				      Magnetic_parameters_pt->h_app(time(), dummy));
      Crystalline_anisotropy_energy = 0.0; //??ds not implemented
      Magnetostatic_energy = 0.0;
      Exchange_energy = 0.0;

      Initial_energy = Zeeman_energy + Crystalline_anisotropy_energy;
    }

    /// Calculate effective damping from midpoint
    double effective_damping(const double& energy, const double& prev_energy)
    {
      // Get dt, make sure it's the one we just used not the one that
      // adaptivity has chosen as the next dt by calculating it from the
      // times.
      const double dt = time_pt()->time(0) - time_pt()->time(1);

      Data* dat_pt = mesh_pt()->element_pt(0)->internal_data_pt(0);

      Vector<double> dmdt(3, 0.0);
      for(unsigned i=0; i<3; i++)
	{
	  dmdt[i] = (dat_pt->value(0, i) - dat_pt->value(1, i))/dt;
	}
      const double dmdtnorm = two_norm(dmdt);

      return -(energy - prev_energy)/(dt*dmdtnorm*dmdtnorm);
    }

    virtual void actions_after_explicit_timestep()
    {
      ODEProblem::actions_after_explicit_timestep();

      // We need to keep M normalised...
      if(should_renormalise())
	{
	  renormalise_magnetisation();
	}
    }

    virtual void actions_before_time_integration()
    {
      // Check that initial condition is unit length
#ifdef PARANOID
      if(m_length_error() > 1e-14)
        {
          std::string err = "Initial condition has wrong |m|.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
    }



    virtual void build(Vector<Mesh*>& bulk_mesh_pts)
    {
      ODEProblem::build(bulk_mesh_pts);

#ifdef PARANOID
      if(Magnetic_parameters_pt == 0)
        {
          std::string err = "Magnetic_parameters_pt is null!";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      // Set parameters in solution
      element_pt()->Exact_solution_pt->initialise_from_problem(this);

      // Rough check for if we can use Mallinson: h is constant in time and
      // axis aligned and non-zero damping.
      double ftol = 1e-12;
      Vector<double> dummy;
      Vector<double> h0 = Magnetic_parameters_pt->h_app(0, dummy);
      Vector<double> h1 = Magnetic_parameters_pt->h_app(0.1, dummy);
      Vector<double> h2 = Magnetic_parameters_pt->h_app(1000000, dummy);
      Mallinson_applicable =
        (two_norm_diff(h0, h1) < ftol)
        && (two_norm_diff(h1, h2) < ftol)
        && (std::abs(h0[0]) < ftol) // no x component
        && (std::abs(h0[1]) < ftol) // no y component
        && (Magnetic_parameters_pt->gilbert_damping() > 1e-12) // some damping
        ;
    }

    double m_length_error() const
    {
      return std::abs(1 - VectorOps::two_norm(solution()));
    }

    virtual void write_additional_trace_headers(std::ofstream& trace_file) const
    {
      trace_file
        << Trace_seperator << "m_length_error_means"
        << Trace_seperator << "m_length_error_std_devs"
        << Trace_seperator << "m_length_error_maxes"
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
        << Trace_seperator << "switching_time_error"
        << Trace_seperator << "energy_change"
        ;
    }

    virtual void write_additional_trace_data(const unsigned& t_hist,
                                             std::ofstream& trace_file) const
    {

      Vector<double> m = solution();

      trace_file
        << Trace_seperator << m_length_error()
        << Trace_seperator << 0
        << Trace_seperator << m_length_error()
        << Trace_seperator << 0
        << Trace_seperator << m[0]
        << Trace_seperator << m[1]
        << Trace_seperator << m[2];

      if(t_hist == 0)
        {
          double total_energy = Exchange_energy + Zeeman_energy
            + Crystalline_anisotropy_energy + Magnetostatic_energy;
          double energy_change = total_energy - Initial_energy;

          double sw_time_error = 0;
            if(Mallinson_applicable)
              {
                sw_time_error = get_switching_time_error_norm(t_hist);
              }
            else
              {
                sw_time_error = MyProblem::Dummy_doc_data;
              }


          trace_file
            << Trace_seperator << Exchange_energy
            << Trace_seperator << Zeeman_energy
            << Trace_seperator << Crystalline_anisotropy_energy
            << Trace_seperator << Magnetostatic_energy
            << Trace_seperator << MyProblem::Dummy_doc_data
            << Trace_seperator << Effective_damping_constant
            << Trace_seperator << Alt_eff_damp
            << Trace_seperator << sw_time_error
            << Trace_seperator << energy_change
            ;
        }
      else
        {
          trace_file
            << Trace_seperator << MyProblem::Dummy_doc_data
            << Trace_seperator << MyProblem::Dummy_doc_data
            << Trace_seperator << MyProblem::Dummy_doc_data
            << Trace_seperator << MyProblem::Dummy_doc_data
            << Trace_seperator << MyProblem::Dummy_doc_data
            << Trace_seperator << MyProblem::Dummy_doc_data
            << Trace_seperator << MyProblem::Dummy_doc_data
            << Trace_seperator << MyProblem::Dummy_doc_data
            << Trace_seperator << MyProblem::Dummy_doc_data;
        }
    }


    virtual double get_error_norm(const unsigned& t_hist=0) const;

    double get_switching_time_error_norm(const unsigned& t_hist) const;

    /// Storage for magnetic parameters object
    MagneticParameters* Magnetic_parameters_pt;

    /// Can we use mallinson solution?
    bool Mallinson_applicable;

    double M_length_error_tol;

    /// When we renormalise m, should we also renormalise any history
    /// values?
    bool Renormalise_history;

  private:

    double Initial_energy;
    double Exchange_energy;
    double Zeeman_energy;
    double Crystalline_anisotropy_energy;
    double Magnetostatic_energy;
    double Effective_damping_constant;
    double Alt_eff_damp;
  };


  class ODECliArgs : public MyCliArgs
  {
  public:

    ODECliArgs() {}
    virtual ~ODECliArgs() {}

    virtual void set_flags()
    {
      MyCliArgs::set_flags();

      specify_command_line_flag("-exact", &exact_name);
      exact_name = "sin";

      specify_command_line_flag("-initial-m", &initial_m_name);
      initial_m_name = "ode_z";

      specify_command_line_flag("-fd-jacobian", &fd_jacobian,
                                "Finite difference the Jacobian, default: -1.");
      fd_jacobian = -1;

      specify_command_line_flag("-sc-beta", &sc_beta, "Self correcting parameter (only for LL equation)");
      sc_beta = 0.0;
    }

    virtual void run_factories()
    {

      exact_name = to_lower(exact_name);
      initial_condition_pt = ODEFactories::exact_solutions_factory(exact_name,
                                                                   initial_m_name,
                                                                   sc_beta);

      initial_is_exact = true;

      MyCliArgs::run_factories();
    }

    /// Overloaded to just create a single ode element
    virtual void build_meshes()
    {
      mesh_pts.push_back(new Mesh);
      mesh_pts[0]->
        add_element_pt(new ODEElement(time_stepper_pt, exact_solution_pt()));
    }

    virtual void assign_specific_parameters(MyProblem* problem_pt) const
    {
      ODEProblem* ode_pt = checked_dynamic_cast<ODEProblem*>(problem_pt);
      if(fd_jacobian != -1)
        {
          ode_pt->Use_fd_jacobian = bool(fd_jacobian);
        }

    }

    std::string exact_name;
    std::string initial_m_name;

    int fd_jacobian;

    double sc_beta;
  };


  class LLGODECliArgs : public ODECliArgs
  {
  public:

    LLGODECliArgs()
    {
      mag_parameters_pt = 0;
    }

    virtual ~LLGODECliArgs() {}

    virtual void set_flags()
    {
      ODECliArgs::set_flags();

      specify_command_line_flag("-damping", &damping, "LLG damping constant.");
      damping = 0.01;

      specify_command_line_flag("-h-app", &h_app_name);
      h_app_name =  "minus_z";

      specify_command_line_flag("-renorm-tol", &renorm_tol);
      renorm_tol = 1e200; // By default never renormalise

      specify_command_line_flag("-also-renorm-history", &also_renorm_history,
				"Should we also renormalise the history values, default: -1.");
      also_renorm_history = -1;
    }

    virtual void run_factories()
    {
      ODECliArgs::run_factories();

      // Create magnetic parameters
      mag_parameters_pt = Factories::magnetic_parameters_factory("simple-llg");
      mag_parameters_pt->Gilbert_damping = damping;
      mag_parameters_pt->Applied_field_fct_pt = Factories::h_app_factory(h_app_name);
    }

    virtual void assign_specific_parameters(MyProblem* problem_pt) const
    {
      ODECliArgs::assign_specific_parameters(problem_pt);

      // Assign magnetic parameters pointer
      LLGODEProblem* llg_ode_pt = checked_dynamic_cast<LLGODEProblem*>(problem_pt);
      llg_ode_pt->Magnetic_parameters_pt = mag_parameters_pt;
      llg_ode_pt->M_length_error_tol = renorm_tol;

      if(also_renorm_history != -1)
	{
	  llg_ode_pt->Renormalise_history = bool(also_renorm_history);
	}
    }

    double damping;
    std::string h_app_name;
    MagneticParameters* mag_parameters_pt;

    double renorm_tol;
    int also_renorm_history;
  };


} // End of oomph namespace

#endif
