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

// Basic oomph-lib headers
#include "generic.h"

#include <ostream>
#include <utility>


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






  struct MyCliArgs
  {
    public:

    MyCliArgs(int argc, char *argv[])
      : nx(10), ny(10), dt(1e-6), tmax(1.0), tol(0.0),
        timestepper("bdf2"), outdir("results")
      {
        // Store command line args
        CommandLineArgs::setup(argc,argv);
        CommandLineArgs::specify_command_line_flag("-nx", &nx);
        CommandLineArgs::specify_command_line_flag("-ny", &ny);

        CommandLineArgs::specify_command_line_flag("-dt", &dt);
        CommandLineArgs::specify_command_line_flag("-tmax", &tmax);
        CommandLineArgs::specify_command_line_flag("-tol", &tol);

        CommandLineArgs::specify_command_line_flag("-outdir", &outdir);

        CommandLineArgs::specify_command_line_flag("-timestepper", &timestepper);

        CommandLineArgs::parse_and_assign();
        CommandLineArgs::output();
        CommandLineArgs::doc_specified_flags();
      }

    // Variables
    unsigned nx;
    unsigned ny;

    double dt;
    double tmax;
    double tol;

    std::string timestepper;

    std::string outdir;

    // Adaptive if eps has been set
    bool adaptive_flag() {return tol != 0.0;}
  };

}

#endif
