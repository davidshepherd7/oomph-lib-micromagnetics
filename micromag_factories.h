#ifndef OOMPH_MICROMAG_FACTORIES_H
#define OOMPH_MICROMAG_FACTORIES_H

#include "../../src/generic/Vector.h"

namespace oomph
{

  // Move into types header?

  /// Type to hold data on which boundaries of which meshes should be
  /// included in the bem.
  class Mesh;
  typedef Vector<std::pair<unsigned, const Mesh*> > BemBoundaryData;

  /// Type to hold data on the locations of sharp corners in bem meshes
  typedef Vector<std::pair<Vector<double>, double> > CornerDataInput;


  // Try to only use forward decls (rather than #includes) in this header
  // because otherwise it will end up huge as it will need to include
  // pretty much everything. Then lots of other files will include this to
  // get access to the factory functions and we end up with everything
  // included in everything -> v. long compile times!
  class BoundaryElementHandler;

  namespace factories
  {

    void bem_handler_factory(const BemBoundaryData& bem_boundaries,
                             const unsigned& phi_index,
                             const unsigned& phi_1_index,
                             const CornerDataInput& input_corner_data_pt,
                             bool use_hlib,
                             BoundaryElementHandler& new_bem_handler);


    BoundaryElementHandler* bem_handler_factory
    (const BemBoundaryData& bem_boundaries,
     const unsigned& phi_index,
     const unsigned& phi_1_index,
     const CornerDataInput& input_corner_data_pt,
     bool use_hlib);

  }
}

#endif
