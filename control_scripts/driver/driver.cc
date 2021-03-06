
#include "generic.h"
#include "micromag.h"

#include <string>

// Floating point error checks
#include <fenv.h>

// For getcwd
#include <unistd.h>
#include <stdio.h>

#include <memory> // unique_ptr

using namespace oomph;
using namespace MathematicalConstants;
using namespace StringConversion;

namespace oomph
{

  namespace DriverFactories
  {

    /// Get the first argument (not the
    std::string get_problem_name(int& argc, char *argv[])
      {
        if(argc < 2)
          {
            std::string err = "No problem type argument found!";
            throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                OOMPH_CURRENT_FUNCTION);
          }

        // Get the problem name and check it
        std::string problem_name(argv[1]);

        if(problem_name[0] == '-')
          {
            std::string err = "Problem name seems to start with -, you probably forgot it.";
            throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                OOMPH_CURRENT_FUNCTION);
          }

        // In place modify argc and argv to remove this argument so that
        // other parsers can read the argument list.
        for(unsigned j=1; j<unsigned(argc); j++)
          {
            argv[j] = argv[j+1];
          }
        argc--;

        return problem_name;
      }


    MyProblem* problem_factory(const std::string& _problem_name)
    {
       const std::string problem_name = to_lower(_problem_name);
       MyProblem* problem_pt = 0;

       if(problem_name == "llg")
         {
           LLGProblem* llgp_pt = new LLGProblem;
           llgp_pt->Residual_calculator_pt
             = Factories::residual_calculator_factory("llg");
           problem_pt = llgp_pt;
         }
       else if(problem_name == "ll")
         {
           LLGProblem* llgp_pt = new LLGProblem;
           llgp_pt->Residual_calculator_pt
             = Factories::residual_calculator_factory("ll");
           problem_pt = llgp_pt;
         }
       else if(problem_name == "ode")
         {
           problem_pt = new ODEProblem;
         }
       else if(problem_name == "llgode") //??ds temp hack...
         {
           problem_pt = new LLGODEProblem;
         }
       else if(problem_name == "poisson")
         {
           problem_pt = new GenericPoissonProblem;
         }

       else
         {
           Factories::unrecognised_name(problem_name, OOMPH_CURRENT_FUNCTION);
         }

       return problem_pt;
    }


    MyCliArgs* args_factory(const std::string& _cli_args_name)
    {
      const std::string cli_args_name = to_lower(_cli_args_name);
      MyCliArgs* cli_args_pt = 0;

      if(cli_args_name == "llg" || cli_args_name == "ll")
        {
          cli_args_pt = new LLGArgs;
        }
      else if(cli_args_name == "ode")
        {
          cli_args_pt = new ODECliArgs;
        }
      else if(cli_args_name == "llgode") //??ds temp hack...
        {
          cli_args_pt = new LLGODECliArgs;
        }
      else if(cli_args_name == "poisson")
        {
          cli_args_pt = new MyCliArgs;
        }
       else
         {
           Factories::unrecognised_name(cli_args_name, OOMPH_CURRENT_FUNCTION);
         }

      return cli_args_pt;
    }

  }

  std::string get_working_path()
  {
    char temp [ PATH_MAX ];

    if( getcwd(temp, PATH_MAX) != 0)
      {
        return std::string ( temp );
      }
    else
      {
        std::string err = "Failed to get cwd.";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
  }

}

using namespace DriverFactories;


int main(int argc, char *argv[])
{
  // Take first argument out of list: it's the problem type name
  std::string problem_name = get_problem_name(argc, argv);

  // Start MPI if necessary
  MPI_Helpers::init(argc, argv);

// #ifdef PARANOID
//   // Enable some floating point error checkers
//   feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
// #endif

  std::cout << "cwd is:" << get_working_path() << std::endl;

  // Create problem class and argument parser
  oomph_info << "Making " << problem_name << " problem." << std::endl;
  std::auto_ptr<MyProblem> problem_pt(problem_factory(problem_name));
  std::auto_ptr<MyCliArgs> args_pt(args_factory(problem_name));
  args_pt->parse(argc, argv);

  oomph_info << "With the following arguments:" << std::endl;
  CommandLineArgs::doc_specified_flags();
  oomph_info << std::endl;


  // Assign time steppers to problem. At most one of these two will be a
  // real timestepper, the other will be null or a dummy.
  problem_pt->set_explicit_time_stepper_pt(args_pt->explicit_time_stepper_pt);

  // Assign general purpose parameters
  problem_pt->linear_solver_pt() = args_pt->solver_pt;
  problem_pt->newton_solver_tolerance() = args_pt->newton_tol;
  problem_pt->max_residuals() = args_pt->newton_max_residual;
  problem_pt->max_newton_iterations() = args_pt->newton_max_iterations;
  problem_pt->maximum_dt() = args_pt->max_dt;
  problem_pt->Exact_solution_pt = args_pt->exact_solution_pt();

  if(args_pt->crash_newton_fail != -1)
    {
      problem_pt->time_adaptive_newton_crash_on_solve_fail()
        = args_pt->crash_newton_fail;
    }

  problem_pt->Error_norm_limit = args_pt->error_norm_limit;
  problem_pt->Solution_norm_limit = args_pt->solution_norm_limit;

  if(args_pt->disable_mm_opt != -1)
    {
      problem_pt->Disable_mass_matrix_solver_optimisations =
        bool(args_pt->disable_mm_opt);
    }

  if(args_pt->dump != -1)
    {
      problem_pt->Dump = bool(args_pt->dump);
    }

  if(args_pt->output_ltes != -1)
    {
      problem_pt->Output_ltes = bool(args_pt->output_ltes);
    }

  if(args_pt->output_predictor_values != -1)
    {
      problem_pt->Output_predictor_values =
        bool(args_pt->output_predictor_values);
    }

  if(args_pt->should_doc_boundaries != -1)
    {
      problem_pt->Should_doc_boundaries = bool(args_pt->should_doc_boundaries);
    }

  // Assign doc parameters
  problem_pt->Doc_info.copy_args_string();
  problem_pt->Doc_info.set_directory(args_pt->outdir);
  problem_pt->Doc_info.output_jacobian = args_pt->output_jacobian;
  problem_pt->set_doc_times(args_pt->doc_times);

  if(args_pt->output_initial_conditions != -1)
    {
      problem_pt->Output_initial_condition = bool(args_pt->output_initial_conditions);
    }

  if(args_pt->doc_exact != -1)
    {
      problem_pt->Want_doc_exact = bool(args_pt->doc_exact);
    }

  if(args_pt->always_write_trace != -1)
    {
      problem_pt->Always_write_trace = args_pt->always_write_trace;
    }

  if(args_pt->predictor_as_initial_guess != -1)
    {
      problem_pt->use_predictor_values_as_initial_guess()
        = bool(args_pt->predictor_as_initial_guess);
    }

  // Assign anything problem specific
  args_pt->assign_specific_parameters(problem_pt.get());


  // Build and initialise the problem
  problem_pt->build(args_pt->mesh_pts);

  //??ds temp hack for mm problems
  LLGArgs* mm_args_pt = dynamic_cast<LLGArgs*>(args_pt.get());
  if(mm_args_pt !=0)
    {
      LLGProblem* llg_pt = checked_dynamic_cast<LLGProblem*>(problem_pt.get());

      if(llg_pt->implicit_ms_flag() || llg_pt->Decoupled_ms)
        {
          // Create the bem handler
          llg_pt->Bem_handler_pt = Factories::bem_handler_factory
            (args_pt->mesh_pts, 0,
             mm_args_pt->hlib_bem,
             false,
             mm_args_pt->numerical_int_bem,
             llg_pt->pin_any_phi_1() || llg_pt->pin_a_boundary_phi_1(),
             mm_args_pt->cached_bem_matrix_filename,
             mm_args_pt->cached_bem_matrix_filename_out
             );
        }
      else
        {
          llg_pt->Bem_handler_pt = 0;
        }
    }

  // Choose which dt to use (calculate from refinement level for convergence
  // test, or use the dtinitial for adaptive or just use dt?).
  double dt;
  if(args_pt->convergence_test != -1
     && bool(args_pt->convergence_test))
    {
      // Assuming that mesh is reasonably uniform...
      double Nx = std::pow(problem_pt->mesh_pt()->nnode(),
                           1.0/double(problem_pt->dim()));

      // Then choose dt as in Jeong2014
      dt = 0.32*(1/Nx);
    }
  else if(args_pt->tol != 0)
    {
      dt = args_pt->dt_initial;
    }
  else
    {
      dt = args_pt->dt;
    }


  // Give the initial condition functor any problem info it needs
  args_pt->initial_condition_pt->initialise_from_problem(problem_pt.get());

  // Get initial condition from either a function pt or a restart file
  if(args_pt->restart_file == "")
    {
      // Set all dts to the value given in args
      problem_pt->initialise_dt(dt);

      // Set values useing the initial condition function
      problem_pt->set_initial_condition(*args_pt->initial_condition_pt);
    }
  else
    {
      std::ifstream restart_file(args_pt->restart_file.c_str());
      problem_pt->read(restart_file);

      if(args_pt->impulsive_restart == 1)
        {
          problem_pt->time_pt()->time() = 0.0;
          problem_pt->initialise_dt(dt);
          problem_pt->set_up_impulsive_initial_condition();
        }
    }

  problem_pt->initial_doc();

  // Do any pre-time integration actions (like relaxing the magnetisation).
  problem_pt->actions_before_time_integration();

  // Solve
  if(problem_pt->is_steady())
    {
      problem_pt->newton_solve();
    }
  else
    {
      // Initialise loop variables
      double tmax = args_pt->tmax;
      double tol = args_pt->tol;
      unsigned max_steps = args_pt->max_steps;


      // Time step to end or to max number of steps
      while((problem_pt->time() < tmax)
            && (problem_pt->N_steps_taken < max_steps)
            && (!problem_pt->finished()))
        {
          // Output some basic info
          oomph_info
            << std::endl
            << std::endl
            << "Time step " << problem_pt->N_steps_taken
            << ", time = " << problem_pt->time()
            << ", dt = " << dt << std::endl
            << "=============================================" << std::endl
            << std::endl;

          // Do the newton solve (different ones depending flags set)
          dt = problem_pt->smart_time_step(dt, tol);

          // Output
          problem_pt->doc_solution();
        }
    }


  problem_pt->final_doc();

  // Maybe have to clean up bem handler, I wish we could use c++11 so we
  // could not have to do this stuff by hand...
  LLGProblem* llg_pt = dynamic_cast<LLGProblem*>(problem_pt.get());
  if(llg_pt != 0)
    {
      delete llg_pt->Bem_handler_pt; llg_pt->Bem_handler_pt = 0;
      delete llg_pt->Residual_calculator_pt; llg_pt->Residual_calculator_pt = 0;
    }

  // Shut down oomph-lib's MPI
  MPI_Helpers::finalize();

  return 0;
}
