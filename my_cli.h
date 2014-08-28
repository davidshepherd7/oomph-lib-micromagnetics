#ifndef OOMPH_MY_CLI_H
#define OOMPH_MY_CLI_H

#include "magnetics_helpers.h"
#include "micromag_types.h"
#include "oomph_factories.h"

// solver factories
#include "../../src/generic/linear_solver.h"
#include "../../src/generic/iterative_linear_solver.h"

// Preconditioner factories
#include "../../src/generic/preconditioner.h"
#include "../../src/generic/general_purpose_preconditioners.h"
#include "../../src/generic/general_purpose_block_preconditioners.h"
#include "sum_of_matrices_preconditioner.h"

namespace oomph
{

  class MyProblem;

  using namespace CommandLineArgs;
  using namespace StringConversion;
  using namespace Factories;



  /// \short Parse inputs and store in a struct-like format. The objects
  /// specified are created using factory functions. Extension to specific
  /// problems can be done by inheriting and overloading set_flags and
  /// run_factories as appropriate.
  class MyCliArgs
  {

  public:

    /// Constructor: Initialise pointers to null.
    MyCliArgs() : time_stepper_pt(0), solver_pt(0), prec_pt(0)
    {
      explicit_time_stepper_pt = 0;
      mesh_factory_pt = 0;
      initial_condition_pt = 0;
    }

    /// Destructor: clean up everything we made in the factories.
    virtual ~MyCliArgs()
      {
        delete prec_pt; prec_pt = 0;
        delete solver_pt; solver_pt = 0;
        delete time_stepper_pt; time_stepper_pt = 0;
        delete explicit_time_stepper_pt; explicit_time_stepper_pt = 0;
        delete initial_condition_pt; initial_condition_pt = 0;
      }

    virtual void set_flags()
      {
        specify_command_line_flag("-dt", &dt, "Time step size");
        dt = 0.1;

        specify_command_line_flag("-tmax", &tmax, "Time at which to stop");
        tmax = 1.0;

        specify_command_line_flag("-tol", &tol,
                                  "Adaptive time step tolerance (default 0 = fixed step)");
        tol = 0.0;

        specify_command_line_flag("-dt-initial", &dt_initial,
                                  "Initial dt for adaptive time step selection, default: 1e-5.");
        dt_initial = 1e-5;

        specify_command_line_flag("-ref", &refinement, "Spatial refinement level of mesh");
        refinement = 1;

        specify_command_line_flag("-newton-tol", &newton_tol, "Newton solver tolerance");
        newton_tol = 1e-8;

        specify_command_line_flag("-newton-max-residual", &newton_max_residual,
                                  "Maximum Newton residual before we give up");
        newton_max_residual = 10.0;

        specify_command_line_flag("-newton-max-iterations", &newton_max_iterations,
                                  "Maximum Newton iterations before we give up");
        newton_max_iterations = 10;

        specify_command_line_flag("-crash-newton-fail", &crash_newton_fail,
                                  "If set then always throw an error on Newton solve fails (normally with adaptive time stepping we try again with a reduced step)");
        crash_newton_fail = -1;

        specify_command_line_flag("-fd-jac", "Use finite differences to calculate the elemental Jacobians");

        specify_command_line_flag("-outdir", &outdir, "Directory to write output to, default is results");
        outdir = "results";

        specify_command_line_flag("-output-jac", &output_jacobian,
                                  "When should we write out Jacobian/mass matrix/residuals, valid settings are 'never' (default), 'always', 'at_start' or 'at_end");
        output_jacobian = "never";

        specify_command_line_flag("-ts", &ts_name, "The time stepper to use, default is bdf2.");
        ts_name = "bdf2";

        specify_command_line_flag("-mp-pred", &mp_pred_name, "The explicit time stepper to use as a predictor for implicit midpoint rule, default is ebdf3.");
        mp_pred_name = "ebdf3";

        specify_command_line_flag("-solver", &solver_name, "The linear solver to use, default is SuperLU.");
        solver_name = "superlu";

        specify_command_line_flag("-matrix-type", &matrix_type,
                                  "The type of Jacobian matrix to use, default: \"cr\".");
        matrix_type = "cr";


        specify_command_line_flag("-prec", &prec_name, "The preconditioner to use for iterative solves.");
        prec_name = "none";

        specify_command_line_flag("-doc-interval", &doc_times_interval,
                                  "The amount of (simulated) time to allow between writing out all data. Default is 0.1. Set to 0 to output at every step. Trace file is usually updated after every step regardless (unless -always-write-trace is set to false).");
        doc_times_interval = 0.1;

        specify_command_line_flag("-mesh", &mesh_name, "The name of the mesh to use.");
        mesh_name = "sq_square";

        specify_command_line_flag("-nnode1d", &nnode1d,
                                  "The number of nodes that will lie along one axis of each element, default is 2 (linear elements).");
        nnode1d = 2;

        specify_command_line_flag("-xshift", &xshift);
        xshift = 1.5;

        specify_command_line_flag("-yshift", &yshift);
        yshift = 1.5;

        specify_command_line_flag("-rotate-xy", &rotate_xy_angle,
                                  "Rotate 2d mesh by angle (in degrees).");
        rotate_xy_angle = 0.0;

        specify_command_line_flag("-scale", &scale);
        scale = 1.0;

        specify_command_line_flag("-blocking", &blocking_name);
        blocking_name = "none";

        specify_command_line_flag("-max-steps", &max_steps);
        max_steps = UINT_MAX; // can't get bigger than this or we overflow

        specify_command_line_flag("-error-norm-limit", &error_norm_limit);
        error_norm_limit = -1.0;

        specify_command_line_flag("-solution-norm-limit", &solution_norm_limit);
        solution_norm_limit = -1.0;

        specify_command_line_flag("-disable-mm-opt");

        specify_command_line_flag("-predictor-as-initial-guess",
                                  &predictor_as_initial_guess);
        predictor_as_initial_guess = -1;

        specify_command_line_flag("-dump", &dump);
        dump = -1;

        specify_command_line_flag("-restart", &restart_file);
        restart_file = "";

        specify_command_line_flag("-impulsive-restart", &impulsive_restart);
        impulsive_restart = -1;

        specify_command_line_flag("-dummy-adaptivity", &dummy_adaptivity);
        dummy_adaptivity = -1;

        specify_command_line_flag("-mp-update-pinned", &mp_update_pinned);
        mp_update_pinned = -1;

        specify_command_line_flag("-krylov-tol", &krylov_tol);
        krylov_tol = -1;

        specify_command_line_flag("-initial-is-exact", &initial_is_exact);
        initial_is_exact = 0;

        specify_command_line_flag("-always-write-trace", &always_write_trace);
        always_write_trace = -1;

        specify_command_line_flag("-output-ltes", &output_ltes,
                                  "Output local truncation errors at each node to ltes*.csv files");
        output_ltes = -1;

        specify_command_line_flag("-output-predictor-values",
                                  &output_predictor_values,
                                  "After each step also output the values calculated by the predictor (to a separate file).");
        output_predictor_values = -1;

        specify_command_line_flag("-output-intial-conditions",
                                  &output_initial_conditions,
                                  "Output the initial conditions as though they were a time step (or as much as possible). Note that this shifts the doc info numbers.");
        output_initial_conditions = -1;

        specify_command_line_flag("-output-boundary-node-positions", &should_doc_boundaries,
                          "Output the positions of boundary nodes?");
        should_doc_boundaries = -1;

        specify_command_line_flag("-convergence-test",
                                  &convergence_test,
                                  "Automatically select dt depending on given refinement to obtain convergence plots (as in Jeong2014).");
        convergence_test = -1;

        specify_command_line_flag("-doc-exact",
                                  &doc_exact,
                                  "Output exact solution if available.");
        doc_exact = -1;
       }

    void parse(int argc, char *argv[])
    {
      // Store command line args
      setup(argc,argv);
      this->set_flags();
      parse_and_assign(argc, argv, true);

      // if(command_line_flag_has_been_set("-dt") &&
      //    command_line_flag_has_been_set("-tol"))
      //   {
      //     std::string err = "Cannot set both dt and adaptive dt tolerance, use -initial-dt to set initial dt for adaptive time step.";
      //     throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
      //                         OOMPH_EXCEPTION_LOCATION);
      //   }

      // Run any processing that needs to be done on the arguments
      // (e.g. creating meshes).
      this->run_factories();
    }

    /// Function to overload to assign any additional parameters which
    /// can't be dealt with here.
    virtual void assign_specific_parameters(MyProblem* problem_pt) const {}

    virtual void run_factories()
    {
      // Make sure all strings are lower case
      ts_name = to_lower(ts_name);
      mp_pred_name = to_lower(mp_pred_name);

      // Build the time stepper
      time_stepper_pt = time_stepper_factory(ts_name, mp_pred_name,
                                             mp_update_pinned);

      if(time_stepper_pt == 0) // failed, so try explicit
        {
          explicit_time_stepper_pt =
            explicit_time_stepper_factory(ts_name);

          // Still need a dummy timestepper to get everything set up
          // without segfaults (ts determines how much storage to create in
          // nodes).
          time_stepper_pt = time_stepper_factory("steady");
        }

      solver_pt = Factories::linear_solver_factory(solver_name, matrix_type,
                                                   krylov_tol);


      // Create and set preconditioner pointer if our solver is iterative.
      IterativeLinearSolver* its_pt
        = dynamic_cast<IterativeLinearSolver*>(solver_pt);
      if(its_pt == 0)
        {
#ifdef PARANOID
          if(prec_name != "none")
            {
              std::string error_msg
                = "Cannot use a preconditioner with a non-iterative solver of type "
                + to_string(solver_name);
              throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
#endif
        }
      else
        {
          its_pt->preconditioner_pt() = this->preconditioner_factory(prec_name);
        }

      doc_times = doc_times_factory(doc_times_interval, tmax);

      // Store boolean flags
      use_fd_jacobian = command_line_flag_has_been_set("-fd-jac");
      disable_mass_matrix_solver_optimisations =
        command_line_flag_has_been_set("-disable-mm-opt");

      // Build the meshes using whatever function the sub class defines
      build_meshes();
    }

    virtual Preconditioner* preconditioner_factory(const std::string& name) const
      {
        return Factories::preconditioner_factory(name);
      }

    // Adaptive if a tolerance has been set
    bool adaptive_flag() {return tol != 0.0;}

    /// Use initial condition if we can, otherwise none so far.
    InitialConditionFct* exact_solution_pt() const
      {
        if(initial_is_exact)
          {
            return initial_condition_pt;
          }
        else
          {
            return 0;
          }
      }

    /// Explcit if an explicit_time_stepper_pt has been set
    bool explicit_flag() const
    {
      // Check that we don't have two "real" timesteppers
      if((!time_stepper_pt->is_steady()) && (explicit_time_stepper_pt != 0))
        {
          std::string err = "Somehow implicit timestepper is steady but explicit one is set!";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
      else
        {
          return explicit_time_stepper_pt != 0;
        }
    }

    virtual void build_meshes()
      {
        this->mesh_pts = build_meshes_helper(mesh_factory_pt,
                                             this->time_stepper_pt);
      }

    Vector<Mesh*> build_meshes_helper(MeshFactoryFctPt mesh_factory_pt,
                                      TimeStepper* time_stepper_pt)
      {

#ifdef PARANOID
        if(mesh_factory_pt == 0)
          {
            std::string err = "Mesh factory pointer is null!";
            throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                OOMPH_CURRENT_FUNCTION);
          }
#endif

        Vector<Mesh*> mesh_pts;

        mesh_name = to_lower(mesh_name);

        // If it has this prefix then make multiple meshes
        if(has_prefix("multi_", mesh_name))
          {
            std::string base_name = rest_of_name("multi_", mesh_name);
            mesh_pts = simple_multimesh_factory(mesh_factory_pt,
                                                base_name, refinement,
                                                time_stepper_pt, xshift,
                                                scale,
                                                nnode1d);
          }

        // Or with "many" prefix make a load of meshes
        else if(has_prefix("many_", mesh_name))
          {
            std::string base_name = rest_of_name("many_", mesh_name);

            mesh_pts = simple_many_multimesh_factory
              (mesh_factory_pt, base_name, refinement,
               time_stepper_pt, xshift, yshift, scale, nnode1d);
          }

        // Otherwise just make a single mesh
        else
          {
            mesh_pts.push_back
              (mesh_factory_pt(mesh_name, refinement, time_stepper_pt,
                                          scale, rotate_xy_angle, nnode1d));
          }

        return mesh_pts;
      }


    // Variables
    double dt;
    double dt_initial;
    double tmax;
    double tol;
    int dummy_adaptivity;
    int convergence_test;

    int refinement;
    bool use_fd_jacobian;
    unsigned max_steps;
    double error_norm_limit;
    double solution_norm_limit;
    bool disable_mass_matrix_solver_optimisations;
    int predictor_as_initial_guess;
    int mp_update_pinned;

    double newton_tol;
    double newton_max_residual;
    unsigned newton_max_iterations;
    int crash_newton_fail;
    double krylov_tol;

    std::string outdir;
    std::string output_jacobian;
    int output_ltes;
    int output_predictor_values;
    int doc_exact;
    int output_initial_conditions;

    InitialConditionFct* initial_condition_pt;
    int initial_is_exact;

    TimeStepper* time_stepper_pt;
    ExplicitTimeStepper* explicit_time_stepper_pt;
    LinearSolver* solver_pt;
    Preconditioner* prec_pt;
    Vector<unsigned> dof_to_block_map;

    // Strings for input to factory functions
    std::string ts_name;
    std::string mp_pred_name;
    std::string solver_name;
    std::string matrix_type;
    std::string prec_name;
    std::string blocking_name;
    std::string mesh_name;


    double doc_times_interval;
    Vector<double> doc_times;
    int always_write_trace;
    int should_doc_boundaries;

    int dump;
    std::string restart_file;
    int impulsive_restart;

    // Mesh parameters
    MeshFactoryFctPt mesh_factory_pt;
    Vector<Mesh*> mesh_pts;
    unsigned nnode1d;
    double xshift;
    double yshift;
    double scale;
    double rotate_xy_angle;

  };



}

#endif
