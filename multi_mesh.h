#ifndef OOMPH_MULTI_MESH_H
#define OOMPH_MULTI_MESH_H


#include "../../src/generic/mesh.h"
#include "../../src/generic/timesteppers.h"
#include "../../src/generic/Vector.h"


namespace oomph
{

  /// Move all nodes of a mesh by (x, y, z).
  // ??ds Move inside Mesh?
  inline void shift_mesh(const double& x, const double& y, const double& z,
                         Mesh* mesh_pt)
  {
#ifdef PARANOID
    // If given shifts in more dimensions than we have then error
    if(((mesh_pt->node_pt(0)->ndim() == 1) && ((y != 0.0) || (z != 0.0)))
       ||
       ((mesh_pt->node_pt(0)->ndim() == 2) && (z != 0.0)))
      {
        std::string err = "Warning tried to shift a mesh in more dimensions than possible";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif

    // For each node shift any positions that exist
    for(unsigned nd=0, nnd=mesh_pt->nnode(); nd<nnd; nd++)
      {
        Node* nd_pt = mesh_pt->node_pt(nd);
        nd_pt->x(0) += x;
        if(nd_pt->ndim() > 1) {nd_pt->x(1) += y;}
        if(nd_pt->ndim() > 2) {nd_pt->x(2) += z;}
      }
  }

  // /// Create two meshes near each other.
  // inline Vector<Mesh*> generate_my_multi_mesh(int refinement_level,
  //                                             TimeStepper* ts_pt,
  //                                             unsigned nnode1d)
  // {
  //   Vector<Mesh*> meshes(2, 0);

  //   meshes[0] = SemiImplicitFactories::llg_mesh_factory
  //     ("sq_square", refinement_level, ts_pt, nnode1d);

  //   meshes[1] = SemiImplicitFactories::llg_mesh_factory
  //     ("sq_square", refinement_level, ts_pt, nnode1d);

  //   shift_mesh(+3.0, 0, 0, meshes[1]);

  //   return meshes;
  // }

} // End of oomph namespace

#endif
