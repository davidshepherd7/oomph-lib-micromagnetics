#ifndef OOMPH_ODE_PROBLEM_H
#define OOMPH_ODE_PROBLEM_H

#include "my_general_header.h"
#include "my_generic_problem.h"
#include "micromag_types.h"

namespace oomph
{
  using namespace MathematicalConstants;
  using namespace StringConversion;
  using namespace VectorOps;


  typedef TimeSpaceToDoubleVectFctPt TimeValueToDoubleVectFctPt;

  namespace deriv_functions
  {

    inline Vector<double> llg(const double& time, const Vector<double>& x)
    {
      Vector<double> values(3, 0.0);
      return values;
    }
    inline Vector<double> dllg(const double& t, const Vector<double>& u)
    {
      throw OomphLibError("Function not yet implemented",
                          OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
      Vector<double> deriv(3, 0.0);
      return deriv;
    }


    inline Vector<double> ll(const double& time, const Vector<double>& x)
    {
      Vector<double> m(3, 0.0);
      m[2] = 1.0;
      m[0] = 0.2;
      normalise(m);
      return m;
    }
    inline Vector<double> dll(const double& t, const Vector<double>& m)
    {
      Vector<double> deriv(3, 0.0), dummy;
      Vector<double> happ = HApp::minus_z(t, dummy);

      Vector<double> mxh = cross(m, happ);
      Vector<double> mxmxh = cross(m, mxh);

      double damping = 1.0;

      for(unsigned j=0; j<3; j++)
        {
          deriv[j] = -(1/(1+damping*damping))*mxh[j] - damping* mxmxh[j];
        }

      return deriv;
    }


    inline Vector<double> cos(const double& time, const Vector<double>& x)
    {
      Vector<double> values(1);
      values[0] = std::cos(time);
      return values;
    }
    inline Vector<double> dcos(const double& t, const Vector<double>& u)
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
    inline Vector<double> dsin(const double& t, const Vector<double>& u)
    {
      Vector<double> deriv(1, 0.0);
      deriv[0] = std::cos(t);
      return deriv;
    }


    inline Vector<double> exp(const double& time, const Vector<double>& x)
    {
      Vector<double> values(1);
      values[0] = std::exp(time);
      return values;
    }
    inline Vector<double> dexp(const double& t, const Vector<double>& u)
    {
      Vector<double> deriv(1, 0.0);
      deriv[0] = u[0];
      return deriv;
    }


    // A polynomial of degree 2
    double b0 = 0.5, b1 = 0, b2 = 1;
    inline Vector<double> poly2(const double& time, const Vector<double>& x)
    {
      Vector<double> values(1);
      values[0] =  b2*time*time + b1*time +b0;
      return values;
    }
    inline Vector<double> dpoly2(const double& t, const Vector<double>& u)
    {
      Vector<double> deriv(1, 0.0);
      deriv[0] = 2*t*b2 + b1;
      return deriv;
    }


    // A polynomial of degree 3
    double a0 = 0.5, a1 = 0, a2 = 0, a3 = 1;
    inline Vector<double> poly3(const double& time, const Vector<double>& x)
    {
      Vector<double> values(1);
      values[0] = a3*time*time*time + a2*time*time + a1*time +a0;
      return values;
    }
    inline Vector<double> dpoly3(const double& t, const Vector<double>& u)
    {
      Vector<double> deriv(1, 0.0);
      deriv[0] = 3*t*t*a3 + 2*t*a2 + a1;
      return deriv;
    }
  }

  namespace ODEFactories
  {

    // Pick an exact solution using a name
    inline TimeSpaceToDoubleVectFctPt exact_solutions_factory
    (const std::string& exact_name)
    {
      if(exact_name == "sin") return &deriv_functions::sin;
      else if(exact_name == "cos") return &deriv_functions::cos;
      else if(exact_name == "exp") return &deriv_functions::exp;
      else if(exact_name == "poly3") return &deriv_functions::poly3;
      else if(exact_name == "poly2") return &deriv_functions::poly2;
      else if(exact_name == "ll") return &deriv_functions::ll;
      else if(exact_name == "llg") return &deriv_functions::llg;
      else
        {
          throw OomphLibError("Unrecognised exact solution " + exact_name,
                              OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
    }

    // Pick a derivative function using a name
    inline TimeValueToDoubleVectFctPt derivative_function_factory
    (const std::string& exact_name)
    {
      if(exact_name == "sin") return &deriv_functions::dsin;
      else if(exact_name == "cos") return &deriv_functions::dcos;
      else if(exact_name == "exp") return &deriv_functions::dexp;
      else if(exact_name == "poly3") return &deriv_functions::dpoly3;
      else if(exact_name == "poly2") return &deriv_functions::dpoly2;
      else if(exact_name == "ll") return &deriv_functions::dll;
      else if(exact_name == "llg") return &deriv_functions::dllg;
      else
        {
          throw OomphLibError("Unrecognised exact solution " + exact_name,
                              OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
    }

  }



  //=====================================================================
  /// Element that integrates single ode with bdf scheme
  //=====================================================================
  class ODEElement : public GeneralisedElement
  {

  public:


    /// Constructor: Pass timestepper
    ODEElement(TimeStepper* timestepper_pt,
               TimeSpaceToDoubleVectFctPt exact_solution_pt,
               TimeValueToDoubleVectFctPt derivative_function_pt)
    {
      Exact_solution_pt = exact_solution_pt;
      Derivative_function_pt = derivative_function_pt;

      Vector<double> dummy(3, 1.0);
      Vector<double> deriv = derivative_function(0, dummy);
      unsigned nvalue = deriv.size();

      add_internal_data(new Data(timestepper_pt, nvalue));
    }

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

      // Use FD for jacobian
      GeneralisedElement::fill_in_jacobian_from_internal_by_fd
        (residuals, jacobian, true);

      // // Or we can use this for problems where f(t, u) = f(t), if better
      // // convergence is needed
      // #warning "jacobian assume no direct u dependence in residual"
      //     jacobian(0, 0) = internal_data_pt(0)->time_stepper_pt()->weight(1, 0);

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
    Vector<double> exact_solution(const double& t)
    {
#ifdef PARANOID
      if(Exact_solution_pt == 0)
        {
          throw OomphLibError("No exact solution function",
                              OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif
      Vector<double> dummy;
      return Exact_solution_pt(t, dummy);
    }

    /// Exact solution
    Vector<double> derivative_function(const double& t,
                                       const Vector<double>& u)
    {
#ifdef PARANOID
      if(Derivative_function_pt == 0)
        {
          throw OomphLibError("No derivative function",
                              OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif
      return Derivative_function_pt(t, u);
    }


  private:

    TimeSpaceToDoubleVectFctPt Exact_solution_pt;
    TimeValueToDoubleVectFctPt Derivative_function_pt;
  };


  class ODECliArgs : public MyCliArgs
  {
  public:

    ODECliArgs() : exact_solution_pt(0), derivative_function_pt(0) {}
    virtual ~ODECliArgs() {}

    virtual void set_flags()
    {
      MyCliArgs::set_flags();

      specify_command_line_flag("-exact", &exact_name);
      exact_name = "sin";

      // specify_command_line_flag("-initial-m", &initial_m_name);
      // initial_m_name = "z";

    }

    virtual void run_factories()
    {

      exact_name = to_lower(exact_name);
      exact_solution_pt = ODEFactories::exact_solutions_factory(exact_name);
      derivative_function_pt = ODEFactories::derivative_function_factory(exact_name);

      MyCliArgs::run_factories();

      initial_condition_fpt = exact_solution_pt;
    }

    /// Overloaded to just create a single ode element
    virtual void build_meshes()
    {
      mesh_pts.push_back(new Mesh);
      mesh_pts[0]->
        add_element_pt(new ODEElement(time_stepper_pt, exact_solution_pt,
                                      derivative_function_pt));
    }

    std::string exact_name;
    TimeSpaceToDoubleVectFctPt exact_solution_pt;
    TimeValueToDoubleVectFctPt derivative_function_pt;

    // for micromag odes only
    std::string initial_m_name;
  };



  class ODEProblem : public MyProblem
  {
  public:

    /// constructor
    ODEProblem() {}

    void build(Vector<Mesh*>& bulk_mesh_pts)
    {
      // Call the underlying build
      MyProblem::build(bulk_mesh_pts);

      // Set up the global mesh and assign equation numbers
      build_global_mesh();
      this->assign_eqn_numbers();
      oomph_info << "Number of equations: " << ndof() << std::endl;
    }

    void set_initial_condition(InitialConditionFctPt ic_fpt)
    {
#ifdef PARANOID
      if(ic_fpt == 0)
        {
          std::string err = "Null initial condition function pointer.";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif

      // Loop over current & previous timesteps
      for (int t=time_stepper_pt()->nprev_values(); t>=0; t--)
        {
          double time = time_pt()->time(t);

          std::cout << "setting IC at time =" << time << std::endl;

          // Get + set the (only) value
          Vector<double> dummy(nvalue(), 1.0);
          Vector<double> values = ic_fpt(time, dummy);

          for(unsigned j=0, nj=nvalue(); j<nj; j++)
            {
              mesh_pt()->element_pt(0)->internal_data_pt(0)
                ->set_value(t, j, values[j]);
            }

        }
    }

    void write_additional_trace_headers(std::ofstream& trace_file) const
    {
      trace_file << Trace_seperator << "exact";
    }

    void write_additional_trace_data(std::ofstream& trace_file) const
    {
      trace_file << Trace_seperator << exact_solution(time());
    }

    double get_error_norm() const
    {
      if(nvalue() == 1)
        {
          return std::abs(trace_value() - exact_solution(time())[0]);
        }
      else
        {
          return MyProblem::Dummy_doc_data;
        }
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

    /// Get trace
    double trace_value() const
    {
      return mesh_pt()->element_pt(0)->internal_data_pt(0)->value(0);
    }

    TimeStepper* ts_pt() const
    {
      return mesh_pt()->element_pt(0)->internal_data_pt(0)->time_stepper_pt();
    }

    unsigned nvalue() const
    {
      return checked_dynamic_cast<ODEElement*>(mesh_pt()->element_pt(0))
        ->nvalue();
    }

    // Output solution
    void doc_solution_additional(std::ofstream& soln_file) const
    {
      Data* dat_pt=mesh_pt()->element_pt(0)->internal_data_pt(0);
        Vector<double> solution(nvalue(), 0.0);
      dat_pt->value(solution);

      std::cout << solution << std::endl;
      soln_file << solution << std::endl;
    }

  };

} // End of oomph namespace

#endif
