//Generic routines
#include "generic.h"

#include "micromag.h"

// Floating point error checks
#include <fenv.h>

// Signal handling
#include <signal.h>


using namespace oomph;
using namespace MathematicalConstants;
using namespace StringConversion;

// Our exact solution function type
typedef double (*FunctionOfTimePt)(double t);

// Our derivative function type
typedef double (*FunctionOfTimeValPt)(double t, double u);


namespace deriv_functions
{
  inline double dcos(double t, double u) {return -1*std::sin(t);}
  inline double dsin(double t, double u) {return std::cos(t);}
  inline double dexp(double t, double u) {return u;}

  double a0 = 0.5, a1 = 0, a2 = 0, a3 = 1;
  inline double poly3(double t) {return a3*t*t*t + a2*t*t + a1*t +a0;}
  inline double dpoly3(double t, double u) {return 3*t*t*a3 + 2*t*a2 + a1;}

  double b0 = 0.5, b1 = 0, b2 = 1;
  inline double poly2(double t) {return b2*t*t + b1*t +b0;}
  inline double dpoly2(double t, double u) {return 2*t*b2 + b1;}
}

namespace ODEFactories
{

  // Pick an exact solution using a name
  inline FunctionOfTimePt exact_solutions_factory(const std::string& exact_name)
  {
    if(exact_name == "sin") return &std::sin;
    else if(exact_name == "cos") return &std::cos;
    else if(exact_name == "exp") return &std::exp;
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
    MyCliArgs::run_factories();

    exact_name = to_lower(exact_name);

    exact_solution_pt = ODEFactories::exact_solutions_factory(exact_name);
    derivative_function_pt = ODEFactories::derivative_function_factory(exact_name);
  }

  virtual void dump_args(std::ostream& out_stream) const
  {
    MyCliArgs::dump_args(out_stream);
    out_stream << "exact_name " << exact_name << std::endl;
  }

  std::string exact_name;
  FunctionOfTimePt exact_solution_pt;
  FunctionOfTimeValPt derivative_function_pt;
};


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
    return Exact_solution_pt(t);
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
      oomph_info << "LLG Number of equations: " << ndof() << std::endl;
    }

  void set_initial_condition()
  {
    // Past history needs to be established for t=time0-deltat, ...
    // Then provide current values (at t=time0) which will also form
    // the initial guess for the first solve at t=time0+deltat

    // Loop over current & previous timesteps
    for (int t=time_stepper_pt()->nprev_values(); t>=0; t--)
      {
        double time = time_pt()->time(t);

        std::cout << "setting IC at time =" << time << std::endl;

        // Get + set the (only) value
        double val = exact_solution(time);
        mesh_pt()->element_pt(0)->internal_data_pt(0)->set_value(t, 0, val);
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


int main(int argc, char* argv[])
{

#ifdef PARANOID
  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
#endif

  // Parse args
  ODECliArgs args;
  args.parse(argc, argv);

  // Create problem
  ODEProblem problem;

  // Assign from args
  problem.linear_solver_pt() = args.solver_pt;
  problem.newton_solver_tolerance() = args.newton_tol;
  // problem.Use_fd_jacobian = args.use_fd_jacobian;
  problem.Use_time_adaptive_newton = args.adaptive_flag();

  // Assign time steppers. Only one of the timesteppers will be a real
  // timestepper, the other will be null or a dummy instead.
  problem.add_time_stepper_pt(args.time_stepper_pt);
  problem.set_explicit_time_stepper_pt(args.explicit_time_stepper_pt);

  // Create mesh and build problem
  Vector<Mesh*> mesh_pts;
  mesh_pts.push_back(new Mesh);
  mesh_pts[0]->
    add_element_pt(new ODEElement(args.time_stepper_pt, args.exact_solution_pt,
                                  args.derivative_function_pt));

  problem.build(mesh_pts);

  // Initialise problem and output
  problem.initialise_dt(args.dt);
  problem.set_initial_condition();

  problem.Doc_info.set_directory(args.outdir);
  problem.Doc_info.Args_pt = &args;
  problem.Doc_info.output_jacobian = args.output_jacobian;
  problem.set_doc_times(args.doc_times);
  problem.initial_doc();

  // All ready: step until completion
  double dt = args.dt;
  unsigned time_step_number = 0;
  while (problem.time() < args.tmax)
    {
      time_step_number++;

      std::cout
        << std::endl
        << std::endl
        << "Time step " << time_step_number << std::endl
        << "==========================" << std::endl
        << "time = " << problem.time()
        << ", dt = " << dt
        << std::endl;

      // Do the newton solve (different ones depending flags set)
      dt = problem.smart_newton_solve(dt, args.tol);

      // Output
      problem.doc_solution();
    }


} // end_of_main
