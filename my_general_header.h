#ifndef OOMPH_MY_GENERAL_HEADER_H
#define OOMPH_MY_GENERAL_HEADER_H

/*
  A header for all my debug and output stuff.
*/


// Include the appropriate version of the pretty print header depending on if we
// are using c++11 or not
#ifdef __GXX_EXPERIMENTAL_CXX0X__
#include "prettyprint.hpp"
#else
#include "prettyprint98.hpp"
#endif

#include <ostream>
#include <climits>

#include "../../src/generic/Vector.h"

#include "../../src/generic/preconditioner.h"


#include "magnetics_helpers.h"
#include "micromag_types.h"
#include "micromag_factories.h"



//??ds tidy up this list

// mesh factories
#include "../../src/generic/mesh.h"

// solver factories
#include "../../src/generic/linear_solver.h"
#include "../../src/generic/iterative_linear_solver.h"

// timestepper factories
#include "../../src/generic/timesteppers.h"
#include "../../src/generic/midpoint_method.h"
#include "../../src/generic/explicit_timesteppers.h"

// Preconditioner factories
#include "../../src/generic/preconditioner.h"
#include "../../src/generic/general_purpose_preconditioners.h"
#include "../../src/generic/general_purpose_block_preconditioners.h"
#include "sum_of_matrices_preconditioner.h"

namespace oomph
{

  class MyProblem;

  /// Given a preconditioner:
  /// 1) if it's a the right type of preconditioner return it
  /// 2) otherwise if its a som preconditioner containing" the right type
  /// of preconditioner then return a pointer to the underlying
  /// preconditioner.
  /// 3) otherwise return null
  template<class T>
  T smart_cast_preconditioner(Preconditioner* prec_pt)
  {
    T bp_pt = dynamic_cast<T> (prec_pt);
    if(bp_pt != 0)
      {
        return bp_pt;
      }
    else
      {
        MainMatrixOnlyPreconditioner* som_main_prec_pt
          = dynamic_cast<MainMatrixOnlyPreconditioner*>(prec_pt);
        if(som_main_prec_pt != 0)
          {
            T ul_bp_pt = dynamic_cast<T>
              (som_main_prec_pt->underlying_preconditioner_pt());
            if(ul_bp_pt != 0)
              {
                return ul_bp_pt;
              }
          }
        else
          {
            return 0;
          }
      }

    // Never get here?
    return 0;
  }


  using namespace CommandLineArgs;
  using namespace StringConversion;
  using namespace Factories;



  inline bool small(const double& test_double)
  {
    return std::abs(test_double) < 1e-5;
  }

  struct RowColVal
  {
  public:
    RowColVal(int row_, int col_, double val_)
      : row(row_), col(col_), val(val_)
    {}

    int row;
    int col;
    double val;

    bool operator<(const RowColVal& other) const
    {
      if (this->row == other.row)
        return (this->col < other.col);
      else
        return (this->row < other.row);
    }
  };


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
      initial_condition_fpt = 0;
    }

    /// Destructor: clean up everything we made in the factories.
    virtual ~MyCliArgs()
      {
        delete prec_pt; prec_pt = 0;
        delete solver_pt; solver_pt = 0;
        delete time_stepper_pt; time_stepper_pt = 0;
        delete explicit_time_stepper_pt; explicit_time_stepper_pt = 0;
      }

    virtual void set_flags()
      {
        specify_command_line_flag("-dt", &dt);
        dt = 1e-6;

        specify_command_line_flag("-tmax", &tmax);
        tmax = 1.0;

        specify_command_line_flag("-tol", &tol);
        tol = 0.0;

        specify_command_line_flag("-ref", &refinement);
        refinement = 1;

        specify_command_line_flag("-newton-tol", &newton_tol);
        newton_tol = 1e-8;

        specify_command_line_flag("-newton-max-residual", &newton_max_residual);
        newton_max_residual = 10.0;

        specify_command_line_flag("-newton-max-iterations", &newton_max_iterations);
        newton_max_iterations = 10;

        specify_command_line_flag("-fd-jac");

        specify_command_line_flag("-outdir", &outdir);
        outdir = "results";

        specify_command_line_flag("-output-jac", &output_jacobian);
        output_jacobian = "never";

        specify_command_line_flag("-ts", &time_stepper_name);
        time_stepper_name = "bdf2";

        specify_command_line_flag("-mp-pred", &mp_pred_name);
        mp_pred_name = "ebdf3";

        specify_command_line_flag("-solver", &solver_name);
        solver_name = "superlu";

        specify_command_line_flag("-prec", &prec_name);
        prec_name = "none";

        specify_command_line_flag("-doc-interval", &doc_times_interval);
        doc_times_interval = 0.1;

        specify_command_line_flag("-mesh", &mesh_name);
        mesh_name = "sq_square";

        specify_command_line_flag("-nnode1d", &nnode1d);
        nnode1d = 2;

        specify_command_line_flag("-xshift", &xshift);
        xshift = 1.5;

        specify_command_line_flag("-yshift", &yshift);
        yshift = 1.5;

        specify_command_line_flag("-scale", &scale);
        scale = 1.0;

        specify_command_line_flag("-blocking", &blocking_name);
        blocking_name = "none";

        specify_command_line_flag("-max-steps", &max_steps);
        max_steps = UINT_MAX; // can't get bigger than this or we overflow

        specify_command_line_flag("-error-norm-limit", &error_norm_limit);
        error_norm_limit = -1.0;

        specify_command_line_flag("-disable-explicit-solver-optimisations");

      }

    void parse(int argc, char *argv[])
    {
      // Store command line args
      setup(argc,argv);
      this->set_flags();
      parse_and_assign(argc, argv, true);
      doc_specified_flags();

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
      time_stepper_name = to_lower(time_stepper_name);
      mp_pred_name = to_lower(mp_pred_name);

      // Build all the pointers to stuff
      time_stepper_pt = time_stepper_factory(time_stepper_name,
                                             mp_pred_name);
      if(time_stepper_pt == 0) // failed, so try explicit
        {
          explicit_time_stepper_pt =
            explicit_time_stepper_factory(time_stepper_name);

          // Still need a dummy timestepper to get everything set up
          // without segfaults (ts determines how much storage to create in
          // nodes).
          time_stepper_pt = time_stepper_factory("steady");
        }

      solver_pt = Factories::linear_solver_factory(solver_name);


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
          // Maybe make a preconditioner which only acts on the main matrix of
          // a sum of matrices.
          if(has_prefix("som-main-", prec_name))
            {
              std::string ul_prec_name = rest_of_name("som-main-", prec_name);
              Preconditioner* ul_prec = Factories::preconditioner_factory(ul_prec_name);
              prec_pt = new MainMatrixOnlyPreconditioner(ul_prec);

              its_pt->preconditioner_pt() = prec_pt;
            }
          // Otherwise just make a normal preconditioner
          else if(prec_name != "none")
            {
              prec_pt = Factories::preconditioner_factory(prec_name);

              its_pt->preconditioner_pt() = prec_pt;
            }
        }


      // If possible then set the dof to block mapping, otherwise check
      // that we didn't want to set one.
      GeneralPurposeBlockPreconditioner<CRDoubleMatrix>* gpbp_pt
        = smart_cast_preconditioner<GeneralPurposeBlockPreconditioner<CRDoubleMatrix>*>
        (prec_pt);

      // Test for case where the preconditioner itself is a block preconditioner
      if(gpbp_pt != 0)
        {
          dof_to_block_map = dof_to_block_factory(blocking_name);
          gpbp_pt->set_dof_to_block_map(dof_to_block_map);
        }
      else
        {
          if(blocking_name != "none")
            {
              std::string err = "Dof to block map specified but cannot be set";
              err += " because preconditioner is not a GeneralPurposeBlockPreconditioner.";
              throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                  OOMPH_CURRENT_FUNCTION);
            }
        }

      doc_times = doc_times_factory(doc_times_interval, tmax);

      // Store boolean flags
      use_fd_jacobian = command_line_flag_has_been_set("-fd-jac");
      disable_explicit_solver_optimisations =
        command_line_flag_has_been_set("-disable-explicit-solver-optimisations");

      // Build the meshes using whatever function the sub class defines
      build_meshes();
    }

    /// Write out all args (in a parseable format) to a stream.
    virtual void dump_args(std::ostream& out_stream) const
    {
      out_stream
        << "initial_dt " << dt << std::endl
        << "tmax " << tmax << std::endl
        << "max_steps " << max_steps << std::endl
        << "tol " << tol << std::endl
        << "refinement " << refinement << std::endl
        << "error-norm-limit " << error_norm_limit << std::endl
        << "disable-explicit-solver-optimisations "
        << disable_explicit_solver_optimisations << std::endl

        << "newton-tol " << newton_tol << std::endl
        << "newton-max-residual " << newton_max_residual << std::endl
        << "newton-max-iterations " << newton_max_iterations << std::endl

        << "outdir " << outdir << std::endl
        << "output_jacobian " << output_jacobian << std::endl

        << "time_stepper " << time_stepper_name << std::endl
        << "mp_pred_name " << mp_pred_name << std::endl
        << "solver_name " << solver_name << std::endl
        << "preconditioner_name " << prec_name <<std::endl
        << "blocking_name " << blocking_name << std::endl

        << "doc_times_interval " << doc_times_interval << std::endl

        << "mesh " << mesh_name << std::endl
        << "nnode1d " << nnode1d << std::endl
        << "xshift " << xshift << std::endl
        << "yshift " << yshift << std::endl
        << "scale " << scale << std::endl
        ;
    }

    // Adaptive if a tolerance has been set
    bool adaptive_flag() {return tol != 0.0;}

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
        this->mesh_pts = build_meshes_helper(mesh_factory_pt);
      }

    Vector<Mesh*> build_meshes_helper(MeshFactoryFctPt mesh_factory_pt)
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
                                          scale, nnode1d));
          }

        return mesh_pts;
      }


    // Variables
    double dt;
    double tmax;
    double tol;
    int refinement;
    bool use_fd_jacobian;
    unsigned max_steps;
    double error_norm_limit;
    bool disable_explicit_solver_optimisations;

    double newton_tol;
    double newton_max_residual;
    unsigned newton_max_iterations;

    std::string outdir;
    std::string output_jacobian;

    InitialConditionFctPt initial_condition_fpt;
    TimeStepper* time_stepper_pt;
    ExplicitTimeStepper* explicit_time_stepper_pt;
    LinearSolver* solver_pt;
    Preconditioner* prec_pt;
    Vector<unsigned> dof_to_block_map;

    // Strings for input to factory functions
    std::string time_stepper_name;
    std::string mp_pred_name;
    std::string solver_name;
    std::string prec_name;
    std::string blocking_name;
    std::string mesh_name;


    double doc_times_interval;
    Vector<double> doc_times;

    // Mesh parameters
    MeshFactoryFctPt mesh_factory_pt;
    Vector<Mesh*> mesh_pts;
    unsigned nnode1d;
    double xshift;
    double yshift;
    double scale;

  };



}

#endif
