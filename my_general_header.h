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

#include "../../src/generic/Vector.h"
#include "../../src/generic/timesteppers.h"
#include "../../src/generic/midpoint_method.h"
#include "../../src/generic/explicit_timesteppers.h"

// meshes for factory
#include "../../src/generic/mesh.h"
#include "multi_mesh.h"

// Solvers for solver factory
#include "../../src/generic/linear_solver.h"
#include "../../src/generic/iterative_linear_solver.h"
#include "../../src/generic/hypre_solver.h"
#include "../../src/generic/sum_of_matrices.h"

// Preconditioners for factory
#include "../../src/generic/general_purpose_block_preconditioners.h"


#include "./magnetics_helpers.h"
#include "./sum_of_matrices_preconditioner.h"

namespace oomph
{

  using namespace CommandLineArgs;
  using namespace StringConversion;

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


  namespace Factories
  {

    inline bool has_prefix(std::string prefix, std::string test_string)
      {
        return test_string.find(prefix) == 0;
      }

    inline std::string rest_of_name(std::string prefix, std::string test_string)
      {
#ifdef PARANOID
        if(!has_prefix(prefix, test_string))
          {
            std::string err = "This string doesn't have the prefix to remove!";
            throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                OOMPH_CURRENT_FUNCTION);
          }
#endif
        return std::string(test_string, prefix.size(), string::npos);
      }

    /// \short Make a timestepper from an input argument. Assumption: this
    /// will be passed into a problem, which will delete the pointer when
    /// it's done.
    inline TimeStepper* time_stepper_factory
    (const std::string& ts_name, const std::string& mp_pred_name="rk4")
    {

      // Always make timestepper adaptive, we can control adaptivity by
      // calling adaptive or non adaptive newton solve.
      bool adaptive_flag = true;

      if(ts_name == "bdf1")
        {
          return new BDF<1>(adaptive_flag);
        }
      else if(ts_name == "bdf2")
        {
          return new BDF<2>(adaptive_flag);
        }
      else if(ts_name == "midpoint")
        {
          MidpointMethod* mp_pt = new MidpointMethod(adaptive_flag);

          if(mp_pred_name == "rk4")
            {
              mp_pt->set_predictor_pt(new RungeKutta<4>);
            }
          else if(mp_pred_name == "lrk4")
            {
              mp_pt->set_predictor_pt(new LowStorageRungeKutta<4>);
            }
          else if (mp_pred_name == "ebdf3")
            {
              mp_pt->set_predictor_pt(new EBDF3);
            }
          else
            {
              throw OomphLibError("Unrecognised predictor name "
                                  + mp_pred_name,
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
          return mp_pt;
        }
      else
        {
          throw OomphLibError("Unrecognised timestepper name " + ts_name,
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }

    /// Simple structure to contain information about how to create and
    /// combine multiple meshes.
    struct ShiftedMeshDetails
    {
      ShiftedMeshDetails()
      {
        mesh_name = "";
        xshift = 0;
        yshift = 0;
        zshift = 0;
      }

      std::string mesh_name;
      double xshift;
      double yshift;
      double zshift;
    };

    /// function pointer type to create a mesh
    typedef Mesh* (MeshFactoryFctPt)(const std::string&, int,
                                     TimeStepper*, unsigned);


    /// Create a vector of meshes with the names given in mesh details and
    /// shifted as also specified in mesh_details.
    inline Vector<Mesh*> multimesh_factory(MeshFactoryFctPt* underlying_factory,
                                           Vector<ShiftedMeshDetails>& mesh_details,
                                           int refinement_level,
                                           TimeStepper* time_stepper_pt,
                                           unsigned nnode1d=2)
    {
      Vector<Mesh*> mesh_pts;

      const unsigned nj = mesh_details.size();
      for(unsigned j=0; j<nj; j++)
        {
          mesh_pts.push_back(underlying_factory(mesh_details[j].mesh_name,
                                                refinement_level,
                                                time_stepper_pt,
                                                nnode1d));

          shift_mesh(mesh_details[j].xshift,
                     mesh_details[j].yshift,
                     mesh_details[j].zshift,
                     mesh_pts[j]);
        }

      return mesh_pts;
    }


    /// Create a pair of meshes near to each other (shifted along x).
    inline Vector<Mesh*> simple_multimesh_factory(MeshFactoryFctPt* underlying_factory,
                                                  const std::string& mesh_name,
                                                  int refinement_level,
                                                  TimeStepper* time_stepper_pt,
                                                  double xshift,
                                                  unsigned nnode1d = 2)
    {
      Vector<ShiftedMeshDetails> inputs(2);
      inputs[0].mesh_name = mesh_name;
      inputs[0].xshift = -xshift;

      inputs[1].mesh_name = mesh_name;
      inputs[1].xshift = +xshift;

      return multimesh_factory(underlying_factory,
                               inputs, refinement_level,
                               time_stepper_pt, nnode1d);
    }


    template<class ELEMENT>
    inline Mesh* surface_mesh_factory(Mesh* bulk_mesh_pt,
                               const Vector<unsigned> &boundaries)
    {
      Mesh* flux_mesh_pt = new Mesh;

      // Loop over boundaries which need a surface mesh
      for(unsigned i=0, ni=boundaries.size(); i<ni; i++)
        {
          // Get boundary number
          unsigned b = boundaries[i];

          // Loop over the bulk elements adjacent to boundary b
          for(unsigned e=0, ne=bulk_mesh_pt->nboundary_element(b); e<ne; e++)
            {
              // Get pointer to the bulk element that is adjacent to boundary b
              FiniteElement* bulk_elem_pt = bulk_mesh_pt->boundary_element_pt(b, e);

              // Get the index of the face of the bulk element e on bondary b
              int face_index = bulk_mesh_pt->face_index_at_boundary(b, e);

              // Build the corresponding prescribed-flux element
              ELEMENT* flux_element_pt = new ELEMENT(bulk_elem_pt, face_index);

              // Add the prescribed-flux element to the surface mesh
              flux_mesh_pt->add_element_pt(flux_element_pt);
            }
        }

      return flux_mesh_pt;
    }


    inline LinearSolver* linear_solver_factory(const std::string& _solver_name)
    {
      const std::string solver_name = to_lower(_solver_name);

      LinearSolver* solver_pt;

      if(solver_name == "superlu")
        { solver_pt = new SuperLUSolver; }
      else if(solver_name == "gmres")
        {
          IterativeLinearSolver* its_pt = new GMRES<CRDoubleMatrix>;
          its_pt->max_iter() = 200;
          solver_pt = its_pt;
        }
      else if(solver_name == "cg")
        {
          IterativeLinearSolver* its_pt = new CG<CRDoubleMatrix>;
          its_pt->max_iter() = 200;
          solver_pt = its_pt;
        }
      else if(solver_name == "fdlu")
        { solver_pt = new FD_LU; }
      else if(solver_name == "som-gmres")
        {
          IterativeLinearSolver* its_pt = new GMRES<SumOfMatrices>;
          its_pt->max_iter() = 200;
          solver_pt = its_pt;
        }
      else
        {
          std::string err("Unrecognised solver name ");
          err += solver_name;
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      return solver_pt;
    }



    inline Preconditioner* block_llg_factory(const std::string &_prec_name)
    {
      // Parse the parameter string
      Vector<std::string> parameters = split_string(to_lower(_prec_name), '-');

#ifdef PARANOID
      if(parameters[0] != "blockllg")
        {
          std::string error_msg = _prec_name + " is not a block llg preconditioner";
          error_msg += "(should begin with blockllg-).";
            throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }
      if(parameters.size() != 5)
        {
          std::string error_msg = "Not enough parameters in llg block string.";
            throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }
#endif

      std::cout << "block llg parameters = " << parameters << std::endl;

      // Pick the basic block structure
      GeneralPurposeBlockPreconditioner<CRDoubleMatrix>* bp_pt = 0;
      if(parameters[1] == "blockexact")
        {
          bp_pt = new ExactBlockPreconditioner<CRDoubleMatrix>;
        }
      else if(parameters[1] == "uppertriangular")
        {
          BlockTriangularPreconditioner<CRDoubleMatrix>* tri_prec_pt
            = new BlockTriangularPreconditioner<CRDoubleMatrix>;
          tri_prec_pt->upper_triangular();
          bp_pt = tri_prec_pt;
        }
      else if(parameters[1] == "lowertriangular")
        {
          BlockTriangularPreconditioner<CRDoubleMatrix>* tri_prec_pt
            = new BlockTriangularPreconditioner<CRDoubleMatrix>;
          tri_prec_pt->lower_triangular();
          bp_pt = tri_prec_pt;
        }
      else if(parameters[1] == "blockdiagonal")
        {
          bp_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;
        }
      else
        {
          throw OomphLibError("Unrecognised block structure setting "
                              + parameters[1], OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }


      // Create + set reordering of blocks: the two dofs given (x = m_x
      // etc.) are put into the first block, the other magnetisation dof is
      // put into the second block. Poisson dofs are in their own blocks
      // (empty for now).
      Vector<unsigned> a(5);
      if(parameters[3] == "xy")
        {
          a[0] = 2; a[1] = 3; // poisson blocks
          a[2] = 0; a[3] = 0; // first m block
          a[4] = 1; // second m block
        }
      else if(parameters[3] == "xz")
        {
          a[0] = 2; a[1] = 3; // poisson blocks
          a[2] = 0; a[4] = 0; // first m block
          a[3] = 1; // second m block
        }
      else if(parameters[3] == "yz")
        {
          a[0] = 2; a[1] = 3; // poisson blocks
          a[4] = 0; a[4] = 0; // first m block
          a[2] = 1; // second m block
        }
      else
        {
          throw OomphLibError("Unrecognised block swapping setting "
                              + parameters[3], OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      bp_pt->set_dof_to_block_map(a);


      // Pick the solver to use on the individual blocks (or on the entire
      // thing if we are using "blockexact").
      if(parameters[4] == "exact")
        {
          // Do nothing--default GeneralPurposeBlockPreconditioner solver
          // is SuperLU.
        }
      // else if(parameters[4] == "amg")
      //   {
      //     bp_pt->set_subsidiary_preconditioner_function();
      //   }
      else
        {
          throw OomphLibError("Unrecognised block structure setting "
                              + parameters[4], OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }


      // Now create the upper-left block sub preconditioner.
      // ============================================================
      GeneralPurposeBlockPreconditioner<CRDoubleMatrix>* sub_bp_pt = 0;
      if(parameters[2] == "blockexact")
        {
          sub_bp_pt = new ExactBlockPreconditioner<CRDoubleMatrix>;
        }
      else if(parameters[2] == "uppertriangular")
        {
          BlockTriangularPreconditioner<CRDoubleMatrix>* tri_prec_pt
            = new BlockTriangularPreconditioner<CRDoubleMatrix>;
          tri_prec_pt->upper_triangular();
          sub_bp_pt = tri_prec_pt;
        }
      else if(parameters[2] == "lowertriangular")
        {
          BlockTriangularPreconditioner<CRDoubleMatrix>* tri_prec_pt
            = new BlockTriangularPreconditioner<CRDoubleMatrix>;
          tri_prec_pt->lower_triangular();
          sub_bp_pt = tri_prec_pt;
        }
      else if(parameters[2] == "blockdiagonal")
        {
          sub_bp_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;
        }
      else if(parameters[2] == "blockantidiagonal")
        {
          sub_bp_pt = new BlockAntiDiagonalPreconditioner<CRDoubleMatrix>;
        }
      else
        {
          std::string err = "Unrecognised sub-preconditioner block structure";
          err += " setting " + parameters[2];
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      // The blocks we want to use it on are just the ones which are in the
      // zeroth block of the master. Find out which ones these are.
      Vector<unsigned> sub_bp_mapping;
      for(unsigned j=0; j<a.size(); j++)
        {if(a[j] == 0) sub_bp_mapping.push_back(j);}

      // And set the master pointer along with this mapping
      sub_bp_pt->turn_into_subsidiary_block_preconditioner(bp_pt,
                                                           sub_bp_mapping);

      // Provide our subsidiary preconditioner to the master for use on
      // block 0, which is the first m block.
      bp_pt->set_subsidiary_preconditioner_pt(sub_bp_pt, 0);

      return bp_pt;
    }


    inline Preconditioner* preconditioner_factory(const std::string &_prec_name)
    {
      const std::string prec_name = to_lower(_prec_name);
      Preconditioner* prec_pt = 0;

      // AMG with optimal poisson paramters
      // ============================================================
      if(prec_name == "poisson-amg")
        {
#ifdef OOMPH_HAS_HYPRE
          HyprePreconditioner* amg_pt = new HyprePreconditioner;
          amg_pt->hypre_method() = HyprePreconditioner::BoomerAMG;

          // Use good Poisson settings for 3D (2D problems should be ok
          // with the same ones I hope...).
          Hypre_default_settings::set_defaults_for_3D_poisson_problem(amg_pt);

          prec_pt = amg_pt;

#else // If no Hypre then give a warning and use exact
          OomphLibWarning("Don't have Hypre, using exact preconditioner.",
                          OOMPH_CURRENT_FUNCTION,OOMPH_EXCEPTION_LOCATION);
          prec_pt = preconditioner_factory("exact");
#endif
        }

        // General purpose AMG (default parameters)
        // ============================================================
        else if(prec_name == "amg")
          {
#ifdef OOMPH_HAS_HYPRE
          HyprePreconditioner* amg_pt = new HyprePreconditioner;
          amg_pt->hypre_method() = HyprePreconditioner::BoomerAMG;
          prec_pt = amg_pt;

#else // If no Hypre then give a warning and use exact
          OomphLibWarning("Don't have Hypre, using exact preconditioner.",
            OOMPH_CURRENT_FUNCTION,OOMPH_EXCEPTION_LOCATION);
          prec_pt = preconditioner_factory("exact");
#endif
        }

        else if(prec_name == "identity")
          { prec_pt = new IdentityPreconditioner; }

        else if(prec_name == "none")
          { prec_pt = 0; }

       else if(prec_name == "exact")
          { prec_pt = new SuperLUPreconditioner; }

       else if(prec_name == "blockexact")
          { prec_pt = new ExactBlockPreconditioner<CRDoubleMatrix>; }


      // if it starts with blockllg then call that factory
       else if(split_string(prec_name, '-')[0] == "blockllg")
         {
           prec_pt = block_llg_factory(prec_name);
         }

        else
          {
            std::string err("Unrecognised preconditioner name ");
            err += prec_name;
            throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }

      return prec_pt;

    }


    /// \short Construct a list of times to output the full solution at
    /// based on command line input in label.
    inline Vector<double> doc_times_factory(const double& doc_interval,
                                            const double &t_max)
      {
        Vector<double> doc_times;

        if(doc_interval == 0)
          {
            // Do nothing: empty vector = output at every step.
          }

        // Otherwise we have a number giving the interval between doc
        // times.
        else
          {
            // Add an output time every "doc_interval" time units until we get
            // to t_max.
            double doc_t = 0.0;
            while(doc_t < t_max)
              {
                doc_times.push_back(doc_t);
                doc_t += doc_interval;
              }
          }

        return doc_times;
      }

  }


  /// \short Parse inputs and store in a struct-like format. The objects
  /// specified are created using factory functions. Extension to specific
  /// problems can be done by inheriting and overloading set_flags and
  /// run_factories as appropriate.
  class MyCliArgs
  {
  public:

    /// Constructor: Initialise pointers to null.
    MyCliArgs() : time_stepper_pt(0), solver_pt(0), prec_pt(0) {}

    /// Destructor: clean up everything we made in the factories.
    virtual ~MyCliArgs()
      {
        delete prec_pt; prec_pt = 0;
        delete solver_pt; solver_pt = 0;
        delete time_stepper_pt; time_stepper_pt = 0;
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

        specify_command_line_flag("-preconditioner", &prec_name);
        prec_name = "none";

        specify_command_line_flag("-doc-interval", &doc_times_interval);
        doc_times_interval = 0.1;
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

    virtual void run_factories()
    {
      // Make sure all strings are lower case
      time_stepper_name = to_lower(time_stepper_name);
      mp_pred_name = to_lower(mp_pred_name);

      // Build all the pointers to stuff
      time_stepper_pt = Factories::time_stepper_factory(time_stepper_name,
                                                        mp_pred_name);
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
      else if(prec_name != "none")
        {
          prec_pt = Factories::preconditioner_factory(prec_name);
          its_pt->preconditioner_pt() = prec_pt;
        }

      doc_times = Factories::doc_times_factory(doc_times_interval, tmax);

      use_fd_jacobian = command_line_flag_has_been_set("-fd-jac");
    }

    /// Write out all args (in a parseable format) to a stream.
    virtual void dump_args(std::ostream& out_stream) const
    {
      out_stream
        << "initial_dt " << dt << std::endl
        << "tmax " << tmax << std::endl
        << "tol " << tol << std::endl
        << "refinement " << refinement << std::endl
        << "newton-tol " << newton_tol << std::endl

        << "outdir " << outdir << std::endl
        << "output_jacobian " << output_jacobian << std::endl

        << "time_stepper " << time_stepper_name << std::endl
        << "solver_name " << solver_name << std::endl
        << "preconditioner_name " << prec_name <<std::endl

        << "doc_times_interval " << doc_times_interval << std::endl;
    }

    // Adaptive if a tolerance has been set
    bool adaptive_flag() {return tol != 0.0;}


    // Variables
    double dt;
    double tmax;
    double tol;
    int refinement;
    double newton_tol;
    bool use_fd_jacobian;

    std::string outdir;
    std::string output_jacobian;

    TimeStepper* time_stepper_pt;
    LinearSolver* solver_pt;
    Preconditioner* prec_pt;

    // Strings for input to factory functions
    std::string time_stepper_name;
    std::string mp_pred_name;
    std::string solver_name;
    std::string prec_name;

    double doc_times_interval;
    Vector<double> doc_times;

  };



}

#endif
