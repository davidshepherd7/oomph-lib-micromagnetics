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
#include "./midpoint_method.h"

// Basic oomph-lib headers
#include "generic.h"

#include <ostream>
#include <utility>

// Meshes
#include "meshes/simple_rectangular_quadmesh.h"


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






  // Need to template so that we can construct the right mesh
  template<class ELEMENT>
  struct MyCliArgs
  {
    public:

    MyCliArgs(int argc, char *argv[])
      : dt(1e-6), tmax(1.0), tol(0.0), refinement(1),

        outdir("results"),

        initial_m_fct_pt(0),
        h_app_fct_pt(0),
        time_stepper_pt(0),
        mesh_pt(0),

        time_stepper_name("bdf2"),
        initial_m_name("z"),
        h_app_name("minus_z"),
        mesh_name("SQ_square")
      {
        // Store command line args
        CommandLineArgs::setup(argc,argv);
        CommandLineArgs::specify_command_line_flag("-dt", &dt);
        CommandLineArgs::specify_command_line_flag("-tmax", &tmax);
        CommandLineArgs::specify_command_line_flag("-tol", &tol);
        CommandLineArgs::specify_command_line_flag("-ref", &refinement);

        CommandLineArgs::specify_command_line_flag("-outdir", &outdir);

        CommandLineArgs::specify_command_line_flag("-ts", &time_stepper_name);

        CommandLineArgs::specify_command_line_flag("-initm", &initial_m_name);
        CommandLineArgs::specify_command_line_flag("-happ", &h_app_name);
        CommandLineArgs::specify_command_line_flag("-mesh", &mesh_name);

        CommandLineArgs::parse_and_assign();
        CommandLineArgs::output();
        CommandLineArgs::doc_specified_flags();

        // Make sure all strings are lower case
        to_lower(time_stepper_name);
        to_lower(initial_m_name);
        to_lower(h_app_name);
        to_lower(mesh_name);

        time_stepper_pt = time_stepper_factory(time_stepper_name, tol);
        initial_m_fct_pt = InitialM::initial_m_factory(initial_m_name);
        h_app_fct_pt = HApp::h_app_factory(h_app_name);
        mesh_pt = mesh_factory(mesh_name, refinement);

        // etc. for precond? meshes?
      }

    // Adaptive if a tolerance has been set
    bool adaptive_flag() {return tol != 0.0;}

    // Variables
    double dt;
    double tmax;
    double tol;
    int refinement;

    std::string outdir;

    InitialM::InitialMFctPt initial_m_fct_pt;
    HApp::HAppFctPt h_app_fct_pt;
    TimeStepper* time_stepper_pt;
    Mesh* mesh_pt;

  private:

    // Strings for input to factory functions
    std::string time_stepper_name;
    std::string h_app_name;
    std::string initial_m_name;
    std::string mesh_name;

    // Make a timestepper from an input argument. Assumption: this will be
    // passed into a problem, which will delete the pointer when it's done.
    TimeStepper* time_stepper_factory(const std::string& ts_name,
                                      double tol)
    {
      if(ts_name == "bdf1")
        {
          return new BDF<1>(adaptive_flag());
        }
      if(ts_name == "bdf2")
        {
          return new BDF<2>(adaptive_flag());
        }
      if(ts_name == "midpoint")
        {
          return new MidpointMethod(adaptive_flag());
          //??ds add access to interp points, fudge factor?
        }
      else
        {
          throw OomphLibError("Unrecognised timestepper name " + ts_name,
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }

    Mesh* mesh_factory(const std::string& mesh_name,
                       int refinement_level)
    {
      if(mesh_name == "SQ_square")
        {
          double lx = 1.0;
          unsigned nx = 5 * std::pow(2, refinement_level);
          return new SimpleRectangularQuadMesh<ELEMENT>
            (nx, nx, lx, lx, time_stepper_pt);
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
