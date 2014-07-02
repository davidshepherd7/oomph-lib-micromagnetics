#ifndef OOMPH_LLG_FACTORIES_H
#define OOMPH_LLG_FACTORIES_H

#include <string>
#include "../../src/generic/Vector.h"

#include "micromag_types.h"


namespace oomph
{

  // Try to only use forward decls (rather than #includes) in this header
  // because otherwise it will end up huge as it will need to include
  // pretty much everything. Then lots of other files will include this to
  // get access to the factory functions and we end up with everything
  // included in everything -> v. long compile times!
  class BoundaryElementHandlerBase;
  class Preconditioner;
  class LLGResidualCalculator;
  class Mesh;
  class Integral;
  class MicromagBEMElementEquations;


  namespace Factories
  {
    /// Make a bem handler object via new. Integers not bools in some
    /// places so that we can use -1 as "unset" and determine a good
    /// default. Only required argument is the meshes.
    BoundaryElementHandlerBase* bem_handler_factory
    (const Vector<Mesh*>& output_mesh_pts,
     const CornerDataInput* input_corner_data_pt=0,
     int hierarchical_bem=-1,
     bool disable_corner_angles=false,
     int numerical_int_bem=-1,
     bool allow_pinned_boundary_values=false);

    /// Create a dof to block mapping for llg block preconditioners based
    /// on _name.
    Vector<unsigned> dof_to_block_factory(const std::string& _name);

    /// \short Make a mesh as specified by an input argument. Refined
    /// according to the given refinement level (in some way appropriate
    /// for that mesh type). Assumption: this will be passed into a
    /// problem, which will delete the pointer when it's done.
    Mesh* llg_mesh_factory(const std::string& _mesh_name,
                           int refinement_level,
                           TimeStepper* time_stepper_pt,
                           double scaling_factor=1.0,
                           unsigned nnode1d = 2);

    LLGResidualCalculator* residual_calculator_factory(const std::string& residual);

    /// \short Create a variable order quadrature object based on the
    /// dimension and shape of the element. Only works for some element
    /// types.
    Integral* variable_order_integrator_factory(const FiniteElement* const el_pt);

    /// \short Create a function to create bem elements based on the
    /// elements used in the bulk mesh.
    BEMElementFactoryFctPt bem_element_factory_factory
    (const FiniteElement* bulk_ele_pt);

    /// \short very simple function: create a new face element of type
    /// ELEMENT.
    template<class ELEMENT>
    MicromagBEMElementEquations* bem_element_factory(FiniteElement* ele,
                                                     const int& face)
    {
      return new ELEMENT(ele, face);
    }

    FluxMeshFactoryFctPt
    mm_flux_mesh_factory_factory(const FiniteElement* bulk_ele_pt);

    /// Construct preconditioner for a complete micromagnetics problem.
    Preconditioner* micromag_preconditioner_factory
    (const std::string& ms_prec, const std::string& llg_prec,
     const std::string& llg_sub_prec);

    /// Construct preconditioner for an llg problem (no magnetostatics).
    Preconditioner* llg_preconditioner_factory
    (const std::string& llg_prec, const std::string& llg_sub_prec);

    /// Construct preconditioner for the 2x2 upper left block of an llg
    /// problem.
    Preconditioner* llg_sub_preconditioner_factory(const std::string& llg_sub_prec);

    /// Pick an applied field function pointer
    HAppFctPt h_app_factory(const std::string& field_name);

    /// Pick an initial magnetisation function pointer. ??ds have to pass
    /// in wave solution parameter here but I don't like it... Need a
    /// **kwargs type construct like python so we can pass through args!
    InitialMFct* initial_m_factory(const std::string& m_name,
                                   const double& wave_solution_c=0);

    /// Create a rescaled reduced integration scheme for this element.
    Integral* rescaled_nodal_quadrature_factory(const FiniteElement* ele_pt,
                                                   const double& mean_size);

    /// Create any other reduced intergration scheme for this element.
    Integral* nodal_quadrature_factory(const FiniteElement* ele_pt,
                                          const double& mean_size);

    /// Factory function to return a null integration scheme.
    inline Integral* no_nodal_quadrature_scheme_factory
    (const FiniteElement* ele_pt, const double& mean_size)
    {
      return 0;
    }

    /// Pick the function which will be used by each element to create an
    /// integration scheme. ugh... factories for factories...
    NodalQuadratureFactoryFctPt
    nodal_quadrature_factory_factory(const std::string& label);

  } // end of Factories

} // End of oomph namespace

#endif
