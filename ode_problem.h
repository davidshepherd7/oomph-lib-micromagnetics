#ifndef OOMPH_ODE_PROBLEM_H
#define OOMPH_ODE_PROBLEM_H

#include "my_general_header.h"
#include "my_generic_problem.h"

namespace oomph
{
  using namespace MathematicalConstants;
  using namespace StringConversion;



  // Our derivative function type
  typedef double (*FunctionOfTimeValPt)(const double& t, const double& u);

  // Our exact solution function type
  typedef InitialConditionFctPt FunctionOfTimePt;


  namespace deriv_functions
  {
    inline Vector<double> cos(const double& time, const Vector<double>& x)
    {
      Vector<double> values(1);
      values[0] = std::cos(time);
      return values;
    }
    inline double dcos(const double& t, const double& u) {return -1*std::sin(t);}


    inline Vector<double> sin(const double& time, const Vector<double>& x)
    {
      Vector<double> values(1);
      values[0] = std::sin(time);
      return values;
    }
    inline double dsin(const double& t, const double& u) {return std::cos(t);}


    inline Vector<double> exp(const double& time, const Vector<double>& x)
    {
      Vector<double> values(1);
      values[0] = std::exp(time);
      return values;
    }
    inline double dexp(const double& t, const double& u) {return u;}


    // A polynomial of degree 2
    double b0 = 0.5, b1 = 0, b2 = 1;
    inline Vector<double> poly2(const double& time, const Vector<double>& x)
    {
      Vector<double> values(1);
      values[0] =  b2*time*time + b1*time +b0;
      return values;
    }
    inline double dpoly2(const double& t, const double& u) {return 2*t*b2 + b1;}


    // A polynomial of degree 3
    double a0 = 0.5, a1 = 0, a2 = 0, a3 = 1;
    inline Vector<double> poly3(const double& time, const Vector<double>& x)
    {
      Vector<double> values(1);
      values[0] = a3*time*time*time + a2*time*time + a1*time +a0;
      return values;
    }
    inline double dpoly3(const double& t, const double& u) {return 3*t*t*a3 + 2*t*a2 + a1;}
  }

  namespace ODEFactories
  {

    // Pick an exact solution using a name
    inline FunctionOfTimePt exact_solutions_factory(const std::string& exact_name)
    {
      if(exact_name == "sin") return &deriv_functions::sin;
      else if(exact_name == "cos") return &deriv_functions::cos;
      else if(exact_name == "exp") return &deriv_functions::exp;
      else if(exact_name == "poly3") return &deriv_functions::poly3;
      else if(exact_name == "poly2") return &deriv_functions::poly2;
      else
        {
          throw OomphLibError("Unrecognised exact solution " + exact_name,
                              OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
    }

    // Pick a derivative function using a name
    inline FunctionOfTimeValPt derivative_function_factory(const std::string& exact_name)
    {
      if(exact_name == "sin") return &deriv_functions::dsin;
      else if(exact_name == "cos") return &deriv_functions::dcos;
      else if(exact_name == "exp") return &deriv_functions::dexp;
      else if(exact_name == "poly3") return &deriv_functions::dpoly3;
      else if(exact_name == "poly2") return &deriv_functions::dpoly2;
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
               FunctionOfTimePt exact_solution_pt,
               FunctionOfTimeValPt derivative_function_pt)
    {
      add_internal_data(new Data(timestepper_pt, 1));
      Exact_solution_pt = exact_solution_pt;
      Derivative_function_pt = derivative_function_pt;
    }

    /// Get residuals
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Get pointer to one-and-only internal data object
      Data* dat_pt = internal_data_pt(0);

      // Get it's value
      double u = dat_pt->value(0);

      // Get timestepper
      TimeStepper* timestepper_pt = dat_pt->time_stepper_pt();

      // Get continuous time
      double t = timestepper_pt->time();

      // Get dudt approximation from timestepper: 1st deriv of 0th value
      double dudt = timestepper_pt->time_derivative(1, dat_pt, 0);

      // Residual is difference between the exact derivative and our
      // timestepper's derivative estimate.
      residuals[0] = derivative_function(t, u) - dudt;
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
      mm(0,0) = 1;
    }

    /// Exact solution
    double exact_solution(const double& t)
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
      return Exact_solution_pt(t, dummy)[0];
    }

    /// Exact solution
    double derivative_function(const double& t,
                               const double& u)
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

    FunctionOfTimePt Exact_solution_pt;
    FunctionOfTimeValPt Derivative_function_pt;
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
    }

    virtual void run_factories()
    {
      // Need exact solution for element constructor so this needs to go
      // first
      exact_name = to_lower(exact_name);
      exact_solution_pt = ODEFactories::exact_solutions_factory(exact_name);
      derivative_function_pt = ODEFactories::derivative_function_factory(exact_name);

      MyCliArgs::run_factories();

      initial_condition_fpt = exact_solution_pt;
    }

    virtual void dump_args(std::ostream& out_stream) const
    {
      MyCliArgs::dump_args(out_stream);
      out_stream << "exact_name " << exact_name << std::endl;
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
    FunctionOfTimePt exact_solution_pt;
    FunctionOfTimeValPt derivative_function_pt;
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
          Vector<double> dummy;
          Vector<double> values = ic_fpt(time, dummy);
          mesh_pt()->element_pt(0)->internal_data_pt(0)
            ->set_value(t, 0, values[0]);
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

    double exact_solution(const double& time) const
    {
      ODEElement* el_pt = checked_dynamic_cast<ODEElement*>
        (mesh_pt()->element_pt(0));
      return el_pt->exact_solution(time);
    }

    /// Error norm: use abs(error in data).
    double global_temporal_error_norm()
    {
      Data* dat_pt=mesh_pt()->element_pt(0)->internal_data_pt(0);

      std::cout << "corrector: " << dat_pt->value(0, 0)
                << " predictor: " << dat_pt->value(4, 0)
                << std::endl;

      return std::abs(ts_pt()->temporal_error_in_value(dat_pt, 0));
    }

    // Output solution
    void doc_solution_additional(std::ofstream &some_file) const
    {
      // Number of plot points
      unsigned npts = 2;

      // Output solution with specified number of plot points per element
      mesh_pt()->output(some_file, npts);
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

  };

} // End of oomph namespace

#endif
