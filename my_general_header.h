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


// // Basic oomph-lib headers
// #include "generic.h"

#include <ostream>
// #include <utility>


// // Floating point debugging
// #include <fenv.h>

// // Quicker to use vector functions
// #include "./vector_helpers.h"


#include "./micromagnetics_element.h"
#include "./magnetics_helpers.h"

// Problems
// #include "./my_generic_problem.h"
// #include "./implicit_llg_problem.h"

// Timesteppers
#include "./midpoint_method.h"

// Meshes
#include "meshes/simple_rectangular_quadmesh.h"
#include "meshes/simple_rectangular_tri_mesh.h"
#include "meshes/simple_cubic_tet_mesh.h"
#include "meshes/tetgen_mesh.h"
#include "meshes/triangle_mesh.h"

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
    Mesh* mesh_factory(const std::string& mesh_name,
                       int refinement_level,
                       TimeStepper* time_stepper_pt,
                       unsigned nnode1d = 2)
    {
      // Lower case mesh names!

      if(mesh_name == "sq_square" && nnode1d == 2)
        {
          double lx = 1.0;
          unsigned nx = 5 * std::pow(2, refinement_level);
          return new SimpleRectangularQuadMesh<QMicromagElement<2,2> >
            (nx, nx, lx, lx, time_stepper_pt);
        }
      else if(mesh_name == "ut_square" && nnode1d == 2)
        {
          return new TriangleMesh<TMicromagElement<2, 2> >
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
          return new SimpleCubicTetMesh<TMicromagElement<3, 2> >
            (nx, ny, nz, lx, ly, lz, time_stepper_pt);
        }
      else if(mesh_name == "ut_cubeoid" && nnode1d == 2)
        {
          return new TetgenMesh<TMicromagElement<3, 2> >
            ("../meshes/cubeoid." + to_string(refinement_level) + ".node",
             "../meshes/cubeoid." + to_string(refinement_level) + ".ele",
             "../meshes/cubeoid." + to_string(refinement_level) + ".face",
             time_stepper_pt);
        }
      else if(mesh_name == "st_cubeoid" && nnode1d == 2)
        {
          double lx = 30, ly = lx, lz = 100;
          unsigned nx = std::pow(2, refinement_level);
          return new SimpleCubicTetMesh<TMicromagElement<3, 2> >
            (nx, nx, int(lz/lx)*nx, lx, ly, lz, time_stepper_pt);
        }
      else if(mesh_name == "ut_sphere" && nnode1d == 2)
        {
          return new TetgenMesh<TMicromagElement<3, 2> >
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
    }

  }


  /// \short Parse inputs and store in a struct-like format. The objects
  /// specified are created using factory functions.
  // ??ds Should be two classes: CliArgs and MMCliArgs...
  class MyCliArgs
  {
  public:

    // Initialise pointers to null
    MyCliArgs()
      : initial_m_fct_pt(0), h_app_fct_pt(0), time_stepper_pt(0), mesh_pt(0) {}

    /// \short Fill in defaults
    void assign_default_values()
    {
      dt = 1e-6;
      tmax = 1.0;
      tol = 0.0;
      refinement = 1;

      outdir = "results";
      output_jacobian = "never";

      initial_m_fct_pt = 0;
      h_app_fct_pt = 0;
      time_stepper_pt = 0;
      mesh_pt = 0;

      time_stepper_name = "bdf2";
      initial_m_name = "z";
      h_app_name = "minus_z";
      mesh_name = "sq_square";
    }

    void parse(int argc, char *argv[])
    {
      // Fill in the default values
      assign_default_values();

      // Store command line args
      CommandLineArgs::setup(argc,argv);
      CommandLineArgs::specify_command_line_flag("-dt", &dt);
      CommandLineArgs::specify_command_line_flag("-tmax", &tmax);
      CommandLineArgs::specify_command_line_flag("-tol", &tol);
      CommandLineArgs::specify_command_line_flag("-ref", &refinement);

      CommandLineArgs::specify_command_line_flag("-outdir", &outdir);
      CommandLineArgs::specify_command_line_flag("-output_jac", &output_jacobian);

      CommandLineArgs::specify_command_line_flag("-ts", &time_stepper_name);

      CommandLineArgs::specify_command_line_flag("-initm", &initial_m_name);
      CommandLineArgs::specify_command_line_flag("-happ", &h_app_name);
      CommandLineArgs::specify_command_line_flag("-mesh", &mesh_name);

      CommandLineArgs::parse_and_assign();
      CommandLineArgs::output();
      CommandLineArgs::doc_specified_flags();

      // Make sure all strings are lower case
      time_stepper_name = to_lower(time_stepper_name);
      initial_m_name = to_lower(initial_m_name);
      h_app_name = to_lower(h_app_name);
      mesh_name = to_lower(mesh_name);

      // Build all the pointers to stuff
      time_stepper_pt = Factories::time_stepper_factory(time_stepper_name, tol);
      initial_m_fct_pt = InitialM::initial_m_factory(initial_m_name);
      h_app_fct_pt = HApp::h_app_factory(h_app_name);
      mesh_pt = Factories::mesh_factory(mesh_name, refinement, time_stepper_pt);

      // etc. for precond?
    }

    /// Write out all args (in a parseable format) to a stream.
    void dump_args(std::ostream& out_stream) const
    {
      out_stream
        << "initial_dt " << dt << std::endl
        << "tmax " << tmax << std::endl
        << "tol " << tol << std::endl
        << "refinement " << refinement << std::endl

        << "outdir " << outdir << std::endl
        << "output_jacobian " << output_jacobian << std::endl

        << "time_stepper " << time_stepper_name << std::endl
        << "initial_m " << initial_m_name << std::endl
        << "h_app " << h_app_name << std::endl
        << "mesh " << mesh_name << std::endl;
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

    InitialM::InitialMFctPt initial_m_fct_pt;
    HApp::HAppFctPt h_app_fct_pt;
    TimeStepper* time_stepper_pt;
    Mesh* mesh_pt;

    // Strings for input to factory functions
    std::string time_stepper_name;
    std::string initial_m_name;
    std::string h_app_name;
    std::string mesh_name;

  };


}

#endif
