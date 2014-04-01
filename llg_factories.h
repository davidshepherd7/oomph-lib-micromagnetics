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
  class BoundaryElementHandler;
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
    BoundaryElementHandler* bem_handler_factory
    (const Vector<Mesh*>& output_mesh_pts, const CornerDataInput* input_corner_data_pt=0,
     int hierarchical_bem=-1, bool disable_corner_angles=false,
     int numerical_int_bem=-1);

    /// Fill in and build an existing bem handler object.
    void bem_handler_factory(BoundaryElementHandler& new_bem_handler,
                             const Vector<Mesh*>& output_mesh_pts,
                             const CornerDataInput* input_corner_data_pt=0,
                             int hierarchical_bem=-1,
                             bool disable_corner_angles=false,
                             int numerical_int_bem=-1);


    /// Create an llg block preconditioner based on a string. Possibly
    /// obsolete ??ds
    Preconditioner* block_llg_factory(const std::string& _prec_name);

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
  }

} // End of oomph namespace

#endif