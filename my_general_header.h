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
using namespace CommandLineArgs;
using namespace StringConversion;


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

      else if((solver_name == "cg-amg") || (solver_name == "poisson"))
        {
#ifdef OOMPH_HAS_HYPRE
          HyprePreconditioner* amg_pt = new HyprePreconditioner;
          amg_pt->hypre_method() = HyprePreconditioner::BoomerAMG;

          // Use good Poisson settings for 3D (2D problems should be ok
          // with the same ones I hope...).
          Hypre_default_settings::set_defaults_for_3D_poisson_problem(amg_pt);

          IterativeLinearSolver* cg_pt = new CG<CRDoubleMatrix>;
          cg_pt->preconditioner_pt() = amg_pt;

          solver_pt = cg_pt;
#else // If no Hypre then give an error
          throw OomphLibError("Don't have Hypre.", OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
#endif
        }

      else if(solver_name == "gmres-amg")
        {
          IterativeLinearSolver* gmres_pt = new GMRES<CRDoubleMatrix>;
          solver_pt = gmres_pt;
#ifdef OOMPH_HAS_HYPRE
          HyprePreconditioner* amg_pt = new HyprePreconditioner;
          amg_pt->hypre_method() = HyprePreconditioner::BoomerAMG;
          gmres_pt->preconditioner_pt() = amg_pt;
#else // If no Hypre then give an error
          OomphLibWarning("Don't have Hypre, using exact preconditioner.",
                              OOMPH_CURRENT_FUNCTION,OOMPH_EXCEPTION_LOCATION);
          gmres_pt->preconditioner_pt() = new SuperLUPreconditioner;
#endif
        }
      else if(solver_name == "fdlu")
        {
          solver_pt = new FD_LU;
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


    Vector<double> doc_times_factory(std::string &label, const double &t_max)
      {
        Vector<double> doc_times;

        if(label == "all")
          {
            // Do nothing: empty vector = output at every step.
          }

        // Otherwise we have a number giving the interval between doc
        // times.
        else
          {
            double doc_interval = std::strtod(label.c_str(), 0);

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
    MyCliArgs() : time_stepper_pt(0) {}

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

        specify_command_line_flag("-outdir", &outdir);
        outdir = "results";

        specify_command_line_flag("-output_jac", &output_jacobian);
        output_jacobian = "never";

        specify_command_line_flag("-ts", &time_stepper_name);
        time_stepper_name = "bdf2";

        specify_command_line_flag("-solver", &solver_name);
        solver_name = "superlu";

        specify_command_line_flag("-doc-interval", &doc_times_interval);
        doc_times_interval = "0.1";
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

      // Build all the pointers to stuff
      time_stepper_pt = Factories::time_stepper_factory(time_stepper_name, tol);
      solver_pt = Factories::linear_solver_factory(solver_name);
      // etc. for precond?

      doc_times = Factories::doc_times_factory(doc_times_interval, tmax);
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
        << "solver_name " << solver_name << std::endl

        << "doc_times_interval " << doc_times_interval << std::endl;
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

    std::string doc_times_interval;
    Vector<double> doc_times;

  };



}

#endif
