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
#include "../../src/generic/general_purpose_preconditioners.h"


#include "./magnetics_helpers.h"
#include "./sum_of_matrices_preconditioner.h"

#include "micromagnetics_element.h" //??ds try to get rid of this?

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

  /// Function type for use as initial condition. Slight overhead of
  /// vectors is worth it even in cases with one dof/space dimensions for
  /// the generality. No overhead for returning a vector due to return
  /// value optimisation.
  typedef Vector<double> (*InitialConditionFctPt)(double t,
                                                  const Vector<double>&x);

  namespace Factories
  {

    /// Check if prefix is a prefix of test_string.
    inline bool has_prefix(const std::string& prefix,
                           const std::string& test_string)
    {
      return test_string.find(prefix) == 0;
    }

    /// Remove the string prefix from the start of test_string. If it isn't
    /// a prefix of it then throw an error.
    inline std::string rest_of_name(const std::string& prefix,
                                    const std::string& test_string)
    {
      if(!has_prefix(prefix, test_string))
        {
          std::string err = "This string doesn't have the prefix to remove!";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }

      return std::string(test_string, prefix.size(), string::npos);
    }

    /// Throw an unrecognised name error from factory function "function".
    inline void unrecognised_name(const std::string& name,
                                  const std::string& function)
    {
      std::string err = "Unrecognised name " + name;
      err += " in factory function " + function + ".";
      throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION, function.c_str());
    }

    /// Make an explicit time stepper
    inline ExplicitTimeStepper* explicit_time_stepper_factory
    (const std::string& ts_name)
    {
      if(ts_name == "rk4")
        {
          return new RungeKutta<4>;
        }
      else if(ts_name == "rk2")
        {
          return new RungeKutta<2>;
        }
      else if(ts_name == "lrk4")
        {
          return new LowStorageRungeKutta<4>;
        }
      else if (ts_name == "ebdf3")
        {
          return new EBDF3;
        }
      else if (ts_name == "euler")
        {
          return new Euler;
        }
      else
        {
          throw OomphLibError("Unrecognised time stepper name "
                              + ts_name,
                              OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
        }
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
          ExplicitTimeStepper* pred_pt = explicit_time_stepper_factory(mp_pred_name);
          mp_pt->set_predictor_pt(pred_pt);
          return mp_pt;
        }
      else if(ts_name == "midpoint-bdf")
        {
          MidpointMethodByBDF* mp_pt = new MidpointMethodByBDF(adaptive_flag);
          ExplicitTimeStepper* pred_pt = explicit_time_stepper_factory(mp_pred_name);
          mp_pt->set_predictor_pt(pred_pt);
          return mp_pt;
        }
      else if(ts_name == "steady")
        {
          // 2 steps so that we have enough space to do reasonable time
          // derivative estimates in e.g. energy derivatives.
          return new Steady<3>;
        }
      else
        {
          return 0;
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
    typedef Mesh* (*MeshFactoryFctPt)(const std::string&, int,
                                     TimeStepper*, double, unsigned);


    /// Create a vector of meshes with the names given in mesh details and
    /// shifted as also specified in mesh_details.
    inline Vector<Mesh*> multimesh_factory(MeshFactoryFctPt underlying_factory,
                                           Vector<ShiftedMeshDetails>& mesh_details,
                                           int refinement_level,
                                           TimeStepper* time_stepper_pt,
                                           double scaling,
                                           unsigned nnode1d)
    {
#ifdef PARANOID
      if(underlying_factory == 0)
        {
          std::string err = "Null underlying mesh factory pointer";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif

      Vector<Mesh*> mesh_pts;

      const unsigned nj = mesh_details.size();
      for(unsigned j=0; j<nj; j++)
        {
          // Build it without scaling
          Mesh* mesh_pt = underlying_factory(mesh_details[j].mesh_name,
                                             refinement_level,
                                             time_stepper_pt,
                                             1.0,
                                             nnode1d);

          // Shift it
          shift_mesh(mesh_details[j].xshift,
                     mesh_details[j].yshift,
                     mesh_details[j].zshift,
                     mesh_pt);

          // Scale it. We do this after shifting so that we are scaling the
          // entire multi-mesh, including gaps between meshes. This is much
          // more intuitive.
          scale_mesh(scaling, mesh_pt);

          // Add to list
          mesh_pts.push_back(mesh_pt);
        }

      return mesh_pts;
    }


    /// Create a pair of meshes near to each other (shifted along x).
    inline Vector<Mesh*> simple_multimesh_factory(MeshFactoryFctPt underlying_factory,
                                                  const std::string& mesh_name,
                                                  int refinement_level,
                                                  TimeStepper* time_stepper_pt,
                                                  double xshift,
                                                  double scaling,
                                                  unsigned nnode1d)
    {
      Vector<ShiftedMeshDetails> inputs(2);
      inputs[0].mesh_name = mesh_name;
      inputs[0].xshift = -xshift;


      inputs[1].mesh_name = mesh_name;
      inputs[1].xshift = +xshift;

      return multimesh_factory(underlying_factory,
                               inputs, refinement_level,
                               time_stepper_pt, scaling, nnode1d);
    }

    // Create lots of meshes near to each other (spaced out along x, y).
    inline Vector<Mesh*> simple_many_multimesh_factory
    (MeshFactoryFctPt underlying_factory, const std::string& mesh_name,
     int refinement_level, TimeStepper* time_stepper_pt,
     double xspacing, double yspacing, double scaling, unsigned nnode1d)
    {
      // Idea is to just make a square (or rect.) grid of meshes, spaced
      // according to xspacing and yspacing.

      // Decided number of meshes
      unsigned n_along_one_direction = 4;
      unsigned n_total = n_along_one_direction*n_along_one_direction;
      Vector<ShiftedMeshDetails> inputs(n_total);

      // Keep track of shift distances
      double basex = 0.0;
      double basey = 0.0;

      // Loop through setting the shift distances for each one
      for(unsigned j=0; j<n_along_one_direction; j++)
        {
          for(unsigned i=0; i<n_along_one_direction; i++)
            {
              // Store the mesh details
              inputs[j*n_along_one_direction+i].mesh_name = mesh_name;
              inputs[j*n_along_one_direction+i].xshift = basex;
              inputs[j*n_along_one_direction+i].yshift = basey;

              // Next column: move along y
              basey += yspacing;
            }

          // Next row: zero the y shift and move along x
          basex += xspacing;
          basey = 0.0;
        }

      // Construct the meshes using the general multimesh factory
      return multimesh_factory(underlying_factory,
                               inputs, refinement_level,
                               time_stepper_pt, scaling, nnode1d);
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

    inline Vector<unsigned> dof_to_block_factory(const std::string& _name)
      {
        const std::string name = to_lower(_name);

        // Make an element to look up indicies from
        TMicromagElement<2,2> dummy_ele;

        const unsigned ndof = dummy_ele.ndof_types(); //??ds unsafe?
        Vector<unsigned> dof_to_block(ndof);

        if(name == "none")
          {
            // identity mapping
            for(unsigned j=0; j<ndof; j++)
              {
                dof_to_block[j] = j;
              }
          }

        // All m values in one block, others left alone.
        // [0, 1, 2, 2, 2, 3, 4]
        else if(name == "group-m")
          {
            unsigned k = 0;

            // Normal phi/phi1
            dof_to_block[dummy_ele.phi_index_micromag()] = k++;
            dof_to_block[dummy_ele.phi_1_index_micromag()] = k++;

            // m all in one block
            for(unsigned j=0; j<3; j++)
              {
                int index = dummy_ele.m_index_micromag(j);
                dof_to_block[index] = k;
              }
            k++;

            // boundary phi/phi1 ??ds assume they are at the end...
            dof_to_block[5] = k++;
            dof_to_block[6] = k++;
          }
        else if(name == "group-m-phi-phi-boundary")
          {
            unsigned k = 0;

            // All phi into one block
            dof_to_block[dummy_ele.phi_index_micromag()] = k;
            dof_to_block[5] = k;
            k++;

            // Simiarly for phi1
            dof_to_block[dummy_ele.phi_1_index_micromag()] = k;
            dof_to_block[6] = k;
            k++;

            // m all in one block
            for(unsigned j=0; j<3; j++)
              {
                int index = dummy_ele.m_index_micromag(j);
                dof_to_block[index] = k;
              }
            k++;
          }
        else
          {
            std::string err = "Unrecognised blocking name";
            err += name;
            throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                OOMPH_CURRENT_FUNCTION);
          }

        return dof_to_block;
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

      else if(prec_name == "ilu0")
        { prec_pt = new ILUZeroPreconditioner<CRDoubleMatrix>; }

      else if(prec_name == "identity")
        { prec_pt = new IdentityPreconditioner; }

      else if(prec_name == "none")
        { prec_pt = 0; }

      else if(prec_name == "exact")
        { prec_pt = new SuperLUPreconditioner; }

      else if(prec_name == "blockexact")
        {
          prec_pt = new ExactBlockPreconditioner<CRDoubleMatrix>;
        }

      else if(prec_name == "blocklt")
        {
          BlockTriangularPreconditioner<CRDoubleMatrix>* bp_pt
            = new BlockTriangularPreconditioner<CRDoubleMatrix>;
          bp_pt->lower_triangular();
          prec_pt = bp_pt;
        }

      else if(prec_name == "blockut")
        {
          BlockTriangularPreconditioner<CRDoubleMatrix>* bp_pt
            = new BlockTriangularPreconditioner<CRDoubleMatrix>;
          bp_pt->upper_triangular();
          prec_pt = bp_pt;
        }


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

      solver_pt = linear_solver_factory(solver_name);


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
              Preconditioner* ul_prec = preconditioner_factory(ul_prec_name);
              prec_pt = new MainMatrixOnlyPreconditioner(ul_prec);

              its_pt->preconditioner_pt() = prec_pt;
            }
          // Otherwise just make a normal preconditioner
          else if(prec_name != "none")
            {
              prec_pt = preconditioner_factory(prec_name);

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

      use_fd_jacobian = command_line_flag_has_been_set("-fd-jac");

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
        << "newton-tol " << newton_tol << std::endl

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
    double newton_tol;
    bool use_fd_jacobian;
    unsigned max_steps;

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
