#ifndef OOMPH_MY_GENERAL_HEADER_H
#define OOMPH_MY_GENERAL_HEADER_H

/*
  A header for all my debug and output stuff.
*/


// Include the appropriate version of the pretty print header depending on if we
// are using c++0x or not
#ifdef __GXX_EXPERIMENTAL_CXX0X__
#include "prettyprint.hpp"
#else
#include "prettyprint98.hpp"
#endif

#include <ostream>

// We need full definition of elements so that we can create meshes of
// them.
#include "./micromag.h"
// #include "./micromagnetics_element.h"

#include "./magnetics_helpers.h"

// Timesteppers
#include "./midpoint_method.h"

// Meshes
#include "meshes/simple_rectangular_quadmesh.h"
#include "meshes/simple_rectangular_tri_mesh.h"
#include "meshes/simple_cubic_tet_mesh.h"
#include "meshes/simple_cubic_mesh.h"
#include "meshes/tetgen_mesh.h"
#include "meshes/triangle_mesh.h"

// Variable order integrators
#include "./variable_order_quadrature.h"

using namespace oomph;

namespace oomph
{


  bool small(const double& test_double)
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

    /// \short Make a timestepper from an input argument. Assumption: this will be
    /// passed into a problem, which will delete the pointer when it's
    /// done.
    TimeStepper* time_stepper_factory(const std::string& ts_name,
                                      double tol)
    {
      bool adaptive_flag = (tol != 0.0);

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
          return new MidpointMethod(adaptive_flag, 2, 0.1);
          //??ds add access to interp points, fudge factor?
        }
      else
        {
          throw OomphLibError("Unrecognised timestepper name " + ts_name,
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }


    /// \short Make a mesh as specified by an input argument. Refined
    /// according to the given refinement level (in some way appropriate
    /// for that mesh type). Assumption: this will be passed into a
    /// problem, which will delete the pointer when it's done.
    Mesh* mesh_factory(const std::string& _mesh_name,
                       int refinement_level,
                       TimeStepper* time_stepper_pt,
                       unsigned nnode1d = 2)
    {
      // Ignore case in mesh names
      const std::string mesh_name = to_lower(_mesh_name);

      // Make the mesh and store a pointer to it
      Mesh* mesh_pt = 0;
      if(mesh_name == "sq_square" && nnode1d == 2)
        {
          double lx = 1.0;
          unsigned nx = 5 * std::pow(2, refinement_level);
          mesh_pt = new SimpleRectangularQuadMesh<QMicromagElement<2,2> >
            (nx, nx, lx, lx, time_stepper_pt);
        }
      else if(mesh_name == "ut_square" && nnode1d == 2)
        {
          mesh_pt = new TriangleMesh<TMicromagElement<2, 2> >
            ("../meshes/square." + to_string(refinement_level) + ".node",
             "../meshes/square." + to_string(refinement_level) + ".ele",
             "../meshes/square." + to_string(refinement_level) + ".poly",
             time_stepper_pt);
        }
      else if(mesh_name == "st_cubeoid" && nnode1d == 2)
        {
          // nmag cubeoid
          double lx = 30, ly = lx, lz = 100;
          unsigned nx = 2 * std::pow(2, refinement_level);
          unsigned ny = nx, nz = std::ceil(lz/lx) * nx;
          mesh_pt = new SimpleCubicTetMesh<TMicromagElement<3, 2> >
            (nx, ny, nz, lx, ly, lz, time_stepper_pt);
        }
      else if(mesh_name == "ut_cubeoid" && nnode1d == 2)
        {
          mesh_pt = new TetgenMesh<TMicromagElement<3, 2> >
            ("../meshes/cubeoid." + to_string(refinement_level) + ".node",
             "../meshes/cubeoid." + to_string(refinement_level) + ".ele",
             "../meshes/cubeoid." + to_string(refinement_level) + ".face",
             time_stepper_pt);
        }
      else if(mesh_name == "st_cubeoid" && nnode1d == 2)
        {
          double lx = 30, ly = lx, lz = 100;
          unsigned nx = std::pow(2, refinement_level);
          mesh_pt = new SimpleCubicTetMesh<TMicromagElement<3, 2> >
            (nx, nx, int(lz/lx)*nx, lx, ly, lz, time_stepper_pt);
        }
      else if(mesh_name == "sq_cubeoid" && nnode1d == 2)
        {
          double lx = 30, ly = lx, lz = 100;
          unsigned nx = std::pow(2, refinement_level);
          mesh_pt = new SimpleCubicMesh<QMicromagElement<3, 2> >
            (nx, nx, int(lz/lx)*nx, lx, ly, lz, time_stepper_pt);
        }
      else if(mesh_name == "ut_sphere" && nnode1d == 2)
        {
          mesh_pt = new TetgenMesh<TMicromagElement<3, 2> >
            ("../meshes/sphere." + to_string(refinement_level) + ".node",
             "../meshes/sphere." + to_string(refinement_level) + ".ele",
             "../meshes/sphere." + to_string(refinement_level) + ".face",
             time_stepper_pt);
        }
      else
        {
          throw OomphLibError("Unrecognised mesh name " + mesh_name,
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      // For some reason we have to call this manually...
      mesh_pt->setup_boundary_element_info();

      // Done: pass out the mesh pointer
      return mesh_pt;
    }

    /// \short Create a variable order quadrature object based on the
    /// dimension and shape of the element. Only works for
    Integral* variable_order_integrator_factory(const FiniteElement* const el_pt)
    {
      if((el_pt->nodal_dimension() == 2) && (el_pt->nvertex_node() == 3))
        {
          return new TVariableOrderGaussLegendre<1>;
        }
      else if((el_pt->nodal_dimension() == 2) && (el_pt->nvertex_node() == 4))
        {
          return new QVariableOrderGaussLegendre<1>;
        }
      else if((el_pt->nodal_dimension() == 3) && (el_pt->nvertex_node() == 4))
        {
          return new TVariableOrderGaussLegendre<2>;
        }
      else if((el_pt->nodal_dimension() == 3) && (el_pt->nvertex_node() == 8))
        {
          return new QVariableOrderGaussLegendre<2>;
        }
      else
        {
          std::string err("Cannot determine element type.\n");
          err += "Maybe it is a higher order element (NNODE_1D > 2)?\n";
          err += "Variable order quadratures are not supported for this case.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }

    template<class ELEMENT>
    Mesh* surface_mesh_factory(Mesh* bulk_mesh_pt,
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


    LinearSolver* linear_solver_factory(const std::string& _solver_name)
    {
      const std::string solver_name = to_lower(_solver_name);

      LinearSolver* solver_pt;
      if(solver_name == "superlu")
        {
          solver_pt = new SuperLUSolver;
        }

      else if(solver_name == "gmres-ilu")
        {
          IterativeLinearSolver* gmres_pt = new GMRES<CRDoubleMatrix>;
          gmres_pt->preconditioner_pt()
            = new ILUZeroPreconditioner<CRDoubleMatrix>;

          solver_pt = gmres_pt;
        }

      else if(solver_name == "cg-ilu")
        {
          IterativeLinearSolver* cg_pt = new CG<CRDoubleMatrix>;
          cg_pt->preconditioner_pt()
            = new ILUZeroPreconditioner<CRDoubleMatrix>;

          solver_pt = cg_pt;
        }

      else if(solver_name == "gmres-amg")
        {
#ifdef OOMPH_HAS_HYPRE
          HyprePreconditioner* amg_pt = new HyprePreconditioner;
          amg_pt->hypre_method() = HyprePreconditioner::BoomerAMG;

          IterativeLinearSolver* gmres_pt = new GMRES<CRDoubleMatrix>;
          gmres_pt->preconditioner_pt() = amg_pt;

          solver_pt = gmres_pt;
#else // If no Hypre then give an error
          throw OomphLibError("Don't have Hypre.",
                              OOMPH_CURRENT_FUNCTION,OOMPH_EXCEPTION_LOCATION);
#endif
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

  }


  /// \short Parse inputs and store in a struct-like format. The objects
  /// specified are created using factory functions. Extension to specific
  /// problems can be done by inheriting and overloading set_flags and
  /// run_factories as appropriate.
  class MyCliArgs
  {
  public:

    /// Constructor: Initialise pointers to null.
    MyCliArgs() : time_stepper_pt(0) {}

    virtual void set_flags()
      {
        CommandLineArgs::specify_command_line_flag("-dt", &dt);
        dt = 1e-6;

        CommandLineArgs::specify_command_line_flag("-tmax", &tmax);
        tmax = 1.0;

        CommandLineArgs::specify_command_line_flag("-tol", &tol);
        tol = 0.0;

        CommandLineArgs::specify_command_line_flag("-ref", &refinement);
        refinement = 1;

        CommandLineArgs::specify_command_line_flag("-outdir", &outdir);
        outdir = "results";

        CommandLineArgs::specify_command_line_flag("-output_jac", &output_jacobian);
        output_jacobian = "never";

        CommandLineArgs::specify_command_line_flag("-ts", &time_stepper_name);
        time_stepper_name = "bdf2";

        CommandLineArgs::specify_command_line_flag("-solver", &solver_name);
        solver_name = "superlu";
      }

    void parse(int argc, char *argv[])
    {
      // Store command line args
      CommandLineArgs::setup(argc,argv);
      this->set_flags();
      CommandLineArgs::parse_and_assign();
      CommandLineArgs::doc_specified_flags();

      // Run any processing that needs to be done on the arguments
      // (e.g. creating meshes).
      this->run_factories();
    }

    virtual void run_factories()
    {
      // Make sure all strings are lower case
      time_stepper_name = to_lower(time_stepper_name);

      // Build all the pointers to stuff
      time_stepper_pt = Factories::time_stepper_factory(time_stepper_name, tol);
      solver_pt = Factories::linear_solver_factory(solver_name);
      // etc. for precond?
    }

    /// Write out all args (in a parseable format) to a stream.
    virtual void dump_args(std::ostream& out_stream) const
    {
      out_stream
        << "initial_dt " << dt << std::endl
        << "tmax " << tmax << std::endl
        << "tol " << tol << std::endl
        << "refinement " << refinement << std::endl

        << "outdir " << outdir << std::endl
        << "output_jacobian " << output_jacobian << std::endl

        << "time_stepper " << time_stepper_name << std::endl
        << "solver_name " << solver_name << std::endl;
    }

    // Adaptive if a tolerance has been set
    bool adaptive_flag() {return tol != 0.0;}

    // Variables
    double dt;
    double tmax;
    double tol;
    int refinement;

    std::string outdir;
    std::string output_jacobian;

    TimeStepper* time_stepper_pt;
    LinearSolver* solver_pt;

    // Strings for input to factory functions
    std::string time_stepper_name;
    std::string solver_name;

  };


  /// Command line args class for llg problems.
  class LLGArgs : public MyCliArgs
  {
    public:

    /// Constructor: Initialise pointers to null.
    LLGArgs() : mesh_pt(0), initial_m_fct_pt(0), h_app_fct_pt(0) {}


    virtual void set_flags()
    {
      MyCliArgs::set_flags();


      CommandLineArgs::specify_command_line_flag("-mesh", &mesh_name);
      mesh_name = "sq_square";

      CommandLineArgs::specify_command_line_flag("-initm", &initial_m_name);
      initial_m_name = "z";

      CommandLineArgs::specify_command_line_flag("-happ", &h_app_name);
      h_app_name = "minus_z";
    }


    virtual void run_factories()
    {
      MyCliArgs::run_factories();

      mesh_name = to_lower(mesh_name);
      initial_m_name = to_lower(initial_m_name);
      h_app_name = to_lower(h_app_name);

      mesh_pt = Factories::mesh_factory(mesh_name, refinement, time_stepper_pt);
      initial_m_fct_pt = InitialM::initial_m_factory(initial_m_name);
      h_app_fct_pt = HApp::h_app_factory(h_app_name);
    }

    /// Write out all args (in a parseable format) to a stream.
    virtual void dump_args(std::ostream& out_stream) const
    {
      MyCliArgs::dump_args(out_stream);

      out_stream
        << "mesh " << mesh_name << std::endl
        << "initial_m " << initial_m_name << std::endl
        << "h_app " << h_app_name << std::endl;
    }

    Mesh* mesh_pt;
    InitialM::InitialMFctPt initial_m_fct_pt;
    HApp::HAppFctPt h_app_fct_pt;

    // Strings for input to factory functions
    std::string mesh_name;
    std::string initial_m_name;
    std::string h_app_name;
  };


}

#endif
