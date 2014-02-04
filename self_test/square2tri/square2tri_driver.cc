
#include "generic.h"
#include "micromag.h"

#include "../../../../src/meshes/simple_rectangular_tri_mesh.h"
#include "../../../../src/meshes/simple_rectangular_quadmesh.h"

namespace oomph
{

  /// The trivial factory function for an ELEMENT.
  template<class ELEMENT>
  inline  FiniteElement* new_element() { return new ELEMENT;}

  /// Create a 2d triangle mesh from a quad mesh. Nodes are reused so don't
  /// delete them.
  inline void trianglify_mesh(const Mesh& square_mesh,
                              TriangleMeshBase& tri_mesh,
                              ElementFactoryFctPt element_factory_fpt)
  {
#ifdef PARANOID
    if(tri_mesh.nelement() != 0)
      {
        std::string err = "non-empty element list in input mesh";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
    if(tri_mesh.nelement() != 0)
      {
        std::string err = "non-empty node list in input mesh";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
    if(square_mesh.finite_element_pt(0)->dim() != 2)
      {
        std::string err = "Only implemented for 2D quads with nnode1d = 2.";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
    if(square_mesh.finite_element_pt(0)->nnode() != 4)
      {
        std::string err = "Only implemented for 2D quads with nnode1d = 2.";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif

    // Just copy node pointers
    for(unsigned nd=0, nnd=square_mesh.nnode(); nd<nnd; nd++)
      {
        Node* nd_pt = square_mesh.node_pt(nd);
        tri_mesh.add_node_pt(nd_pt);
      }


    // Create new elements
    for(unsigned ele=0, nele=square_mesh.nelement(); ele<nele; ele++)
      {
        const FiniteElement* qele_pt = square_mesh.finite_element_pt(ele);

#ifdef PARANOID
        if(dynamic_cast<const FaceElement*>(qele_pt) != 0)
          {
            throw OomphLibError("Function not yet implemented for face elements",
                                OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
          }
#endif

        // Make tri elements
        FiniteElement* tele1_pt = element_factory_fpt();
        FiniteElement* tele2_pt = element_factory_fpt();

        // Add to mesh
        tri_mesh.add_element_pt(tele1_pt);
        tri_mesh.add_element_pt(tele2_pt);


        // Set element properties
        tele1_pt->set_n_node(3);
        tele1_pt->set_dimension(2);
        tele1_pt->set_nnodal_position_type(qele_pt->nnodal_position_type());

        tele2_pt->set_n_node(3);
        tele2_pt->set_dimension(2);
        tele2_pt->set_nnodal_position_type(qele_pt->nnodal_position_type());


        // Add nodes to elements
        tele1_pt->node_pt(0) = qele_pt->node_pt(1);
        tele1_pt->node_pt(1) = qele_pt->node_pt(2);
        tele1_pt->node_pt(2) = qele_pt->node_pt(0);

        tele2_pt->node_pt(0) = qele_pt->node_pt(1);
        tele2_pt->node_pt(1) = qele_pt->node_pt(2);
        tele2_pt->node_pt(2) = qele_pt->node_pt(3);
      }

    // Use TriangleMeshBase to sort out boundary stuff
    tri_mesh.set_nboundary(square_mesh.nboundary());
    tri_mesh.setup_boundary_element_info();

  }


  /// Templated version which just creates elements using new.
  template<class ELEMENT>
  inline void trianglify_mesh(const Mesh& square_mesh,
                              TriangleMeshBase&  tri_mesh)
  {
    // trianglify_mesh(square_mesh, tri_mesh, &new_a);
    trianglify_mesh(square_mesh, tri_mesh, &new_element<ELEMENT>);
  }

}


using namespace oomph;

// Test it
int main()
{

  // Make a simple square mesh
  SimpleRectangularQuadMesh<QTFPoissonElement<2, 2> > square_mesh(10, 10, 2, 2);

  // Convert
  TriangleMeshBase* tri_mesh_pt = new TriangleMeshBase;
  trianglify_mesh<TTFPoissonElement<2, 2> >(square_mesh, *tri_mesh_pt);

  // Dump
  tri_mesh_pt->dump("Validation/tri_mesh", false);

  // Make a real tri mesh in same geometry to compare with
  SimpleRectangularTriMesh<TTFPoissonElement<2, 2> > real_tri_mesh(10, 10, 2, 2);
  real_tri_mesh.dump("Validation/real_tri_mesh", false);

  return 0;
}
