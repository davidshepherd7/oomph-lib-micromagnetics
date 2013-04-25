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


// Basic oomph-lib headers
#include "generic.h"

#include <ostream>
#include <utility>


// Floating point debugging
#include <fenv.h>

// Quicker to use vector functions
#include "./vector_helpers.h"


// All my micromag headers and .cc
#include "./micromagnetics_element.h"
#include "./micromagnetics_element.cc"
#include "./micromagnetics_boundary_element.h"
#include "./micromagnetics_boundary_element.cc"
#include "./micromagnetics_flux_element.h"
#include "./magnetic_materials.h"
#include "./micromagnetics_preconditioners.h"
#include "./magnetics_helpers.h"

// Timesteppers
#include "./midpoint_method.h"
#include "./old_midpoint_method.h"

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


  class MyDocInfo : public DocInfo
  {
  public:
    /// Default constructor
    MyDocInfo() : DocInfo(), output_jacobian("never") {}
    MyDocInfo(const std::string& directory,
              const std::string& output_jacobian="never")
    : DocInfo(directory), output_jacobian(output_jacobian)
      {}

    std::string output_jacobian;
  };


  // Need to template so that we can construct the right mesh
  template<class ELEMENT> class MyCliArgs
  {
    public:

    MyCliArgs(int argc, char *argv[])
      : dt(1e-6), tmax(1.0), tol(0.0), refinement(1),

        outdir("results"),
        output_jacobian("never"),

        initial_m_fct_pt(0),
        h_app_fct_pt(0),
        time_stepper_pt(0),
        mesh_pt(0),

        time_stepper_name("bdf2"),
        initial_m_name("z"),
        h_app_name("minus_z"),
        mesh_name("sq_square")
      {
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

        time_stepper_pt = time_stepper_factory(time_stepper_name, tol);
        initial_m_fct_pt = InitialM::initial_m_factory(initial_m_name);
        h_app_fct_pt = HApp::h_app_factory(h_app_name);
        mesh_pt = mesh_factory(mesh_name, refinement);

        // etc. for precond?
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


  private:

    // Strings for input to factory functions
    std::string time_stepper_name;
    std::string initial_m_name;
    std::string h_app_name;
    std::string mesh_name;


    /// \short Make a timestepper from an input argument. Assumption: this will be
    /// passed into a problem, which will delete the pointer when it's
    /// done.
    TimeStepper* time_stepper_factory(const std::string& ts_name,
                                      double tol)
    {
      if(ts_name == "bdf1")
        {
          return new BDF<1>(adaptive_flag());
        }
      else if(ts_name == "bdf2")
        {
          return new BDF<2>(adaptive_flag());
        }
      else if(ts_name == "midpoint")
        {
          return new MidpointMethod(adaptive_flag());
          //??ds add access to interp points, fudge factor?
        }
      else if(ts_name == "bdf1_midpoint")
        {
          return new BDFMidpointMethod(adaptive_flag());
          //??ds add access to interp points, fudge factor?
        }
      else
        {
          throw OomphLibError("Unrecognised timestepper name " + ts_name,
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }


    /// \short Make a mesh as specified by an input argument. Refined according to
    /// the given refinement leve (in some way appropriate for that mesh
    /// type). Assumption: this will be passed into a problem, which will
    /// delete the pointer when it's done.
    Mesh* mesh_factory(const std::string& mesh_name,
                       int refinement_level)
    {
      // Lower case mesh names!

      if(mesh_name == "sq_square")
        {
          double lx = 1.0;
          unsigned nx = 5 * std::pow(2, refinement_level);
          return new SimpleRectangularQuadMesh<ELEMENT>
            (nx, nx, lx, lx, time_stepper_pt);
        }
      else if(mesh_name == "ut_square")
        {
          return new TriangleMesh<ELEMENT>
            ("../meshes/square." + to_string(refinement_level) + ".node",
             "../meshes/square." + to_string(refinement_level) + ".ele",
             "../meshes/square." + to_string(refinement_level) + ".poly",
             time_stepper_pt);
        }
      else if(mesh_name == "st_cubeoid")
        {
          double lx = 30, ly = lx, lz = 100;
          unsigned nx = std::pow(2, refinement_level);
          return new SimpleCubicTetMesh<ELEMENT>
            (nx, nx, int(lz/lx)*nx, lx, ly, lz, time_stepper_pt);
        }
      else
        {
          throw OomphLibError("Unrecognised mesh name " + mesh_name,
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }

  };

}

#endif
