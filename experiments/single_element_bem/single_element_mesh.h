#ifndef OOMPH_SINGLE_ELEMENT_MESH_H
#define OOMPH_SINGLE_ELEMENT_MESH_H

/*
  description of file goes here
*/

// // Config header generated by autoconfig
// #ifdef HAVE_CONFIG_H
// #include <oomph-lib-config.h>
// #endif

// OOMPH-LIB headers
#include "generic/mesh.h"
#include "generic/matrices.h"
#include "generic/quadtree.h"
#include "generic/quad_mesh.h"


namespace oomph
{

 // ============================================================
 ///
 // ============================================================
 template <class ELEMENT>
 class SingleTetElementMesh : public TetMeshBase
 {
 public:

  /// Default constructor
  SingleTetElementMesh
  (const double &Lx, const double &Ly, const double &Lz,
   TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

  /// Destructor
  ~SingleTetElementMesh() {}

 private:

  /// Inaccessible copy constructor
  SingleTetElementMesh(const SingleTetElementMesh &dummy)
  {BrokenCopy::broken_copy("SingleTetElementMesh");}

  /// Inaccessible assignment operator
  void operator=(const SingleTetElementMesh &dummy)
  {BrokenCopy::broken_assign("SingleTetElementMesh");}
 };

 template <class ELEMENT>
 SingleTetElementMesh<ELEMENT>::SingleTetElementMesh
 (const double &Lx, const double &Ly, const double &Lz,
  TimeStepper* time_stepper_pt)
 {

  // Set the number of boundaries
  set_nboundary(4);

  // Create the element
  Element_pt.resize(1);
  Element_pt[0] = new ELEMENT;

  // Read out the number of linear points in the element
  unsigned n_p = dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();

  if(n_p != 2)
   {
    throw OomphLibError("nnode_1d > 2 is not yet implemented,
OOMPH_CURRENT_FUNCTION,
                        "SingleTetElementMesh::SingleTetElementMesh()",
                        OOMPH_EXCEPTION_LOCATION);
   }

  // Can now allocate the store for the nodes
  Node_pt.resize((1 + (n_p-1))*(1 + (n_p-1)));


  // Start creating the element
  // ============================================================

  // Note: have to use stupid c-arrays here because we don't have
  // initialiser lists for c++ vectors yet (they're in c++11).

  // Make an array of positions
  unsigned nnode = 4;
  double node_positions[][3] = {{0,0,0}, {1,0,0}, {0,0,1}, {0,1,0}};

  // Make an array of boundaries that the nodes are on
  unsigned boundaries[][3] = {{0,1,2}, {0,2,3}, {0,1,3}, {1,2,3}};

  // Build corner nodes
  for(unsigned nd=0; nd<nnode; nd++)
   {
    // Make the node. All nodes are boundary nodes here.
    Node_pt[nd] =
     finite_element_pt(0)->construct_boundary_node(nd,time_stepper_pt);

    // Add pointer to element
    finite_element_pt(0)->node_pt(nd) = Node_pt[nd];

    // Set positions
    Node_pt[nd]->x(0) = node_positions[nd][0];
    Node_pt[nd]->x(1) = node_positions[nd][1];
    Node_pt[nd]->x(2) = node_positions[nd][2];

    for(unsigned b=0; b < 3; b++)
     {
      add_boundary_node(boundaries[nd][b], Node_pt[nd]);
     }
   }

  // Not implementing additional nodes: don't need them and I have nothing
  // to test it on.

  // Finish boundary setup
  setup_boundary_element_info();
 }


} // End of oomph namespace

#endif
