
#include "generic.h"
#include "micromag.h"

#include <string>

// Floating point error checks
#include <fenv.h>

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
             = LLGFactories::residual_calculator_factory("llg");
           problem_pt = llgp_pt;
         }
       else if(problem_name == "ll")
         {
           LLGProblem* llgp_pt = new LLGProblem;
           llgp_pt->Residual_calculator_pt
             = LLGFactories::residual_calculator_factory("ll");
           problem_pt = llgp_pt;
         }
       else if(problem_name == "ode")
         {
           problem_pt = new ODEProblem;
         }
       else if(problem_name == "poisson")
         {
           problem_pt = new GenericPoissonProblem;
         }
       else if(problem_name == "heat")
         {
           problem_pt = new UnsteadyHeatProblem;
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
          cli_args_pt = new MMArgs;
        }
      else if(cli_args_name == "ode")
        {
          cli_args_pt = new ODECliArgs;
        }
      else if(cli_args_name == "poisson")
        {
          cli_args_pt = new MyCliArgs;
        }
      else if(cli_args_name == "heat")
        {
          cli_args_pt = new UnsteadyHeatArgs;
        }

       else
         {
           Factories::unrecognised_name(cli_args_name, OOMPH_CURRENT_FUNCTION);
         }

      return cli_args_pt;
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

#ifdef PARANOID
  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
#endif

  // Create problem class and argument parser
  oomph_info << "Making " << problem_name << " problem." << std::endl;
  MyProblem* problem_pt = problem_factory(problem_name);
  MyCliArgs* args_pt = args_factory(problem_name);
  args_pt->parse(argc, argv);

  oomph_info << "With the following arguments:" << std::endl;
  args_pt->dump_args(*oomph_info.stream_pt());
  oomph_info << std::endl;


  // Assign time steppers to problem. At most one of these two will be a
  // real timestepper, the other will be null or a dummy.
  problem_pt->add_time_stepper_pt(args_pt->time_stepper_pt);
  problem_pt->set_explicit_time_stepper_pt(args_pt->explicit_time_stepper_pt);

  // Assign general purpose parameters
  problem_pt->linear_solver_pt() = args_pt->solver_pt;
  problem_pt->newton_solver_tolerance() = args_pt->newton_tol;
  // problem_pt->Use_fd_jacobian = args_pt->use_fd_jacobian; //??ds
  problem_pt->Error_norm_limit = args_pt->error_norm_limit;

  // Assign doc parameters
  problem_pt->Doc_info.copy_args_string(args_pt);
  problem_pt->Doc_info.set_directory(args_pt->outdir);
  problem_pt->Doc_info.output_jacobian = args_pt->output_jacobian;
  problem_pt->set_doc_times(args_pt->doc_times);

  // Assign anything problem specific
  args_pt->assign_specific_parameters(problem_pt);


  // Build and initialise the problem
  MMArgs* mm_args_pt = dynamic_cast<MMArgs*>(args_pt);
  if(mm_args_pt !=0 && mm_args_pt->decoupled_ms)
    {
      //?? clean this up
      LLGProblem* llg_pt = checked_dynamic_cast<LLGProblem*>(problem_pt);
      llg_pt->build_decoupled_ms(mm_args_pt->mesh_pts,
                                 mm_args_pt->phi_mesh_pts,
                                 mm_args_pt->phi_1_mesh_pts);
    }
  else
    {
      problem_pt->build(args_pt->mesh_pts);
    }
  problem_pt->initialise_dt(args_pt->dt); //??ds is this ok for steady state prob?
  problem_pt->set_initial_condition(args_pt->initial_condition_fpt);
  problem_pt->initial_doc();

  // Initialise loop variables
  double dt = args_pt->dt, tmax = args_pt->tmax;
  double tol = args_pt->tol;
  unsigned time_step_number = 0, max_steps = args_pt->max_steps;

  // Solve
  if(problem_pt->is_steady())
    {
      problem_pt->newton_solve();
    }
  else
    {
      // Time step to end or to max max number of steps
      while((problem_pt->time() < tmax) && (time_step_number < max_steps))
        {
          time_step_number++;

          // Output some basic info
          oomph_info
            << std::endl
            << std::endl
            << "Time step " << time_step_number
            << ", time = " << problem_pt->time()
            << ", dt = " << dt << std::endl
            << "=============================================" << std::endl
            << std::endl;

          // Do the newton solve (different ones depending flags set)
          dt = problem_pt->smart_newton_solve(dt, tol);

          // Output
          problem_pt->doc_solution();
        }
    }


  problem_pt->final_doc();

  // Done with problem and args pointers
  delete problem_pt; problem_pt = 0;
  delete args_pt; args_pt = 0;

  // Shut down oomph-lib's MPI
  MPI_Helpers::finalize();

  return 0;
}
