#ifndef OOMPH_SINGLE_ELEMENT_MESH_H
#define OOMPH_SINGLE_ELEMENT_MESH_H


namespace oomph
{

  template <class ELEMENT>
  class SingleElementMesh : public QuadMeshBase
  {
  public:

    SingleElementMesh(TimeStepper* time_stepper_pt=&Mesh::Default_TimeStepper)
    {
      // Create element
      ELEMENT* ele_pt = new ELEMENT;

      // Mesh can only be built with 2D Qelements with 2 nodes per side
      MeshChecker::assert_geometric_element<QElementGeometricBase,ELEMENT>(2, 2);

      // Put into mesh
      this->add_element_pt(ele_pt);

      set_nboundary(4);
      Node_pt.resize(4);

      // Copy nodes into node pt
      for(unsigned nd=0, nnd=ele_pt->nnode(); nd<nnd; nd++)
        {
          Node_pt[nd] = ele_pt->construct_boundary_node(nd, time_stepper_pt);
          ele_pt->node_pt(nd) = Node_pt[nd];
        }

      node_pt(0)->x(0) = 0.0;
      node_pt(0)->x(1) = 0.0;
      add_boundary_node(0, node_pt(0));
      add_boundary_node(3, node_pt(0));

      node_pt(1)->x(0) = 1.0;
      node_pt(1)->x(1) = 0.0;
      add_boundary_node(0, node_pt(1));
      add_boundary_node(1, node_pt(1));

      node_pt(2)->x(0) = 0.0;
      node_pt(2)->x(1) = 1.0;
      add_boundary_node(1, node_pt(2));
      add_boundary_node(2, node_pt(2));

      node_pt(3)->x(0) = 1.0;
      node_pt(3)->x(1) = 1.0;
      add_boundary_node(2, node_pt(3));
      add_boundary_node(3, node_pt(3));

    }
  };
} // End of oomph namespace

#endif
