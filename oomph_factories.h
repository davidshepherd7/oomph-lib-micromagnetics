#ifndef OOMPH_FACTORIES_H
#define OOMPH_FACTORIES_H

#include "micromag_types.h"
#include <string>


// Includes needed for (templated) surface mesh factory.
#include "../../src/generic/mesh.h"
#include "../../src/generic/elements.h"


// Header containing various factory functions (functions which create
// objects).

// The idea is that we often want to be able to create various different
// objects (or different configurations of one object) depending on the
// command line options. We also want this code to be easily reusable in
// various different places (e.g. when creating objects for self tests). So
// we wrap up all the object creation details inside functions.

namespace oomph
{

  // Try to only use forward decls (rather than #includes) in this header
  // because otherwise it will end up huge as it will need to include
  // pretty much everything. Then lots of other files will include this to
  // get access to the factory functions and we end up with everything
  // included in everything -> v. long compile times!
  class Mesh;
  class TimeStepper;
  class ExplicitTimeStepper;
  class LinearSolver;
  class Preconditioner;
  class Integral;



  namespace Factories
  {

    // First some simple helpers
    // ============================================================

    /// Throw an unrecognised name error from factory function "function".
    inline void unrecognised_name(const std::string& name,
                                  const std::string& function)
    {
      if(name == "")
        {
          std::string err = "Got the empty string as a name";
          err += " in factory function " + function + ".";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION, function.c_str());
        }
      else
        {
          std::string err = "Unrecognised name " + name;
          err += " in factory function " + function + ".";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION, function.c_str());
        }
    }

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

      return std::string(test_string, prefix.size(), std::string::npos);
    }


    // General purpose factories
    // ============================================================

    /// \short Make a timestepper from an input argument using
    /// new. Assumption: this will be passed into a problem, which will
    TimeStepper* time_stepper_factory(const std::string& ts_name,
                                      const std::string& mp_pred_name="ebdf3",
                                      const int& mp_update_pinned=-1);

    /// Make an explicit time stepper
    ExplicitTimeStepper* explicit_time_stepper_factory(const std::string& ts_name);


    /// Create a new linear solver as specified by the solver name and
    /// matrix type.
    LinearSolver* linear_solver_factory(const std::string& _solver_name,
                                        const std::string& matrix_type,
                                        const double& krylov_tol,
                                        const unsigned& max_iter,
                                        const bool& throw_on_convergence_fail);

    /// Create a new (general-purpose) preconditioner as specified by the
    /// name.
    Preconditioner* preconditioner_factory(const std::string& _prec_name);

    /// \short Construct a list of times to output the full solution at
    /// based on command line input in label.
    Vector<double> doc_times_factory(const double& doc_interval,
                                     const double &t_max);



    // // Problem factories
    // // ============================================================

    // /// Fill in and build
    // void my_problem_factory_helper(TimeStepper* ts_pt,
    //                                ExplicitTimeStepper* expl_ts_pt,
    //                                LinearSolver* linear_solver_pt,
    //                                const double& newton_tol,
    //                                const double& error_norm_limit,
    //                                bool disable_explicit_solver_optimisations,
    //                                const std::string& outdir,
    //                                const std::string& args_info,
    //                                const std::string& output_jacobian,
    //                                const Vector<double>& doc_times,
    //                                MyProblem& new_problem);


    // LLGProblem* llg_problem_factory();


    // Composite mesh creation
    // ============================================================

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

    /// Create a vector of meshes with the names given in mesh details and
    /// shifted as also specified in mesh_details.
    Vector<Mesh*> multimesh_factory(MeshFactoryFctPt underlying_factory,
                                    Vector<ShiftedMeshDetails>& mesh_details,
                                    int refinement_level,
                                    TimeStepper* time_stepper_pt,
                                    double scaling,
                                    unsigned nnode1d);

    /// Create a pair of meshes near to each other (shifted along x).
    Vector<Mesh*> simple_multimesh_factory(MeshFactoryFctPt underlying_factory,
                                           const std::string& mesh_name,
                                           int refinement_level,
                                           TimeStepper* time_stepper_pt,
                                           double xshift,
                                           double scaling,
                                           unsigned nnode1d);

    /// Create lots of meshes near to each other (spaced out along x, y).
    Vector<Mesh*> simple_many_multimesh_factory
    (MeshFactoryFctPt underlying_factory,
     const std::string& mesh_name,
     int refinement_level,
     TimeStepper* time_stepper_pt,
     double xspacing, double yspacing, double scaling, unsigned nnode1d);


    /// Construct a surface mesh of elements ELEMENT on boundaries listed
    /// in "boundaries" on the mesh given by bulk_mesh_pt. Unfortunately
    /// because it's templated we have to include the entire definition in
    /// the header or limit the template parameters to a pre-built list.
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


    /// Construct a Gaussian quadrature object
    Integral* gauss_integration_factory(const unsigned& dim,
                                        const unsigned& nnode_1d,
                                        const ElementGeometry::ElementGeometry& elgeom);


    /// Construct an element
    template<class ELEMENT> FiniteElement* general_element_factory()
    {
      return new ELEMENT;
    }

  }
}

#endif
