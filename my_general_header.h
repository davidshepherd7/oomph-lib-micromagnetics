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


  class MyDocInfo
  {

  public:

    /// \short Constructor. Default settings: Current directory, step `0',
    /// label="", full documentation enabled and output directory is not required
    /// to exist when set_directory() is called.
    MyDocInfo() :
      Directory("."), Doc_flag(true), Label(""), Timestep(0),
      Newton_step(0), Directory_must_exist(false)
    {}

    /// Output directory
    std::string directory() const {return Directory;}

    /// \short Set output directory (we try to open a file in it
    /// to see if the directory exists -- if it doesn't we'll
    /// issue a warning -- or, if directory_must_exist()==true,
    /// throw and OomphLibError
    void set_directory(const std::string& directory);

    /// \short Enable documentation
    void enable_doc() {Doc_flag=true;}

    /// \short Disable documentation
    void disable_doc() {Doc_flag=false;}

    /// \short Are we documenting?
    bool is_doc_enabled() const {return Doc_flag;}

    //??ds
    unsigned timestep() const {return Timestep;}
    void next_timestep() {Timestep++; Newton_step =0;}

    unsigned newton_step() const {return Newton_step;}
    void next_newton_step() {Newton_step++;}

    /// String used (e.g.) for labeling output files
    std::string& label() {return Label;}

    /// String used (e.g.) for labeling output files. Const version.
    std::string label() const {return Label;}

    /// \short Call to throw an error if directory does not exist
    void enable_error_if_directory_does_not_exist() {Directory_must_exist=true;}

    /// \short Call to issue a warning if the directory does not exists
    void disable_error_if_directory_does_not_exist() {Directory_must_exist=false;}

  private:

    /// Directory name
    std::string Directory;

    /// Doc or don't?
    bool Doc_flag;

    /// String to label output file, say
    std::string Label;

    //??ds
    unsigned Timestep;
    unsigned Newton_step;

    /// Boolean flag to decide response if an output
    /// directory doesn't exist: If true, we terminate
    /// code execution by throwing an OomphLibError rather than
    /// just issuing a warning.
    bool Directory_must_exist;
  };

  void MyDocInfo::set_directory(const std::string& directory_)
  {
    // Try to open a file in output directory
    std::ostringstream filename;
    filename << directory_ << "/.dummy_check.dat";
    std::ofstream some_file;
    some_file.open(filename.str().c_str());
    if (!some_file.is_open())
      {
        //Construct the error message
        std::string error_message = "Problem opening output file.\n";
        error_message += "I suspect you haven't created the output directory ";
        error_message += directory_;
        error_message += "\n";

        //Issue a warning if the directory does not have to exist
        if(!Directory_must_exist)
          {
            //Create an Oomph Lib warning
            OomphLibWarning(error_message,"set_directory()",
                            OOMPH_EXCEPTION_LOCATION);
          }
        //Otherwise throw an erro
        else
          {
            error_message += "and the Directory_must_exist flag is true.\n";
            throw OomphLibError(error_message,
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
      }
    //Write to the dummy file
    some_file << "Dummy file, opened to check if output directory " << std::endl;
    some_file << "exists. Can be deleted...." << std::endl;
    some_file.close();
    // Set directory
    Directory=directory_;
  }


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
