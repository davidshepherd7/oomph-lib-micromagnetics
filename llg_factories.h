#ifndef OOMPH_LLG_FACTORIES_H
#define OOMPH_LLG_FACTORIES_H

#include <string>

namespace oomph
{

  // Try to only use forward decls (rather than #includes) in this header
  // because otherwise it will end up huge as it will need to include
  // pretty much everything. Then lots of other files will include this to
  // get access to the factory functions and we end up with everything
  // included in everything -> v. long compile times!
  class BoundaryElementHandler;
  class Preconditioner;


  namespace Factories
  {

    // LLG specific factories
    // ============================================================

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
  }

} // End of oomph namespace

#endif
