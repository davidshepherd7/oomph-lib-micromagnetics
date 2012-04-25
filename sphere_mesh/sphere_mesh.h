#ifndef OOMPH_COMPLETE_SPHERE_MESH_H
#define OOMPH_COMPLETE_SPHERE_MESH_H

/*
  description of file goes here
*/

#include "generic.h"

#include "meshes/eighth_sphere_mesh.h"

using namespace oomph;

namespace oomph
{

  template<class ELEMENT>
  class SphereMesh : public virtual BrickMeshBase
  {

  public:

    /// \short Constructor: Pass radius and timestepper; defaults to
    /// static default timestepper
    SphereMesh(const double& radius, TimeStepper* time_stepper_pt=
	       &Mesh::Default_TimeStepper);

    Mesh* sub_mesh_pt(const unsigned& i)
    {
      return Sub_mesh_pt[i];
    }


  protected :

    /// Radius of the sphere
    double Radius;

    /// Pointers to the sub-meshes used in construction
    Vector<Mesh*> Sub_mesh_pt;

  };


  template<class ELEMENT>
  SphereMesh<ELEMENT>::SphereMesh(const double& radius,
				  TimeStepper* time_stepper_pt) :
    Radius(radius), Sub_mesh_pt(8)
  {

    // Set up number of boundarys, number of elements
    this->set_nboundary(1);

    // Create the sub-meshes and add pointer to the list of pointers
    for(unsigned i=0; i<8; i++)
      Sub_mesh_pt[i] = new EighthSphereMesh<ELEMENT>(radius, time_stepper_pt);

    // Move nodes so that each submesh is a different sector of a sphere.
    // =================================================================

    // Do nothing to submesh 0

    // Flip submesh 1 over the x axis:
    for(unsigned nd=0; nd < Sub_mesh_pt[1]->nnode(); nd++)
      {
	Sub_mesh_pt[1]->node_pt(nd)->x(0) *= -1;
      }

    // Flip submesh 2 over the x axis and y axis:
    for(unsigned nd=0; nd < Sub_mesh_pt[2]->nnode(); nd++)
      {
	Sub_mesh_pt[2]->node_pt(nd)->x(0) *= -1;
	Sub_mesh_pt[2]->node_pt(nd)->x(1) *= -1;
      }

    // Flip submesh 3 over the z axis:
    for(unsigned nd=0; nd < Sub_mesh_pt[3]->nnode(); nd++)
      {
	Sub_mesh_pt[3]->node_pt(nd)->x(2) *= -1;
      }

    // Flip submesh 4 over the y axis:
    for(unsigned nd=0; nd < Sub_mesh_pt[4]->nnode(); nd++)
      {
	Sub_mesh_pt[4]->node_pt(nd)->x(1) *= -1;
      }

    // Flip submesh 5 over the z axis and x axis:
    for(unsigned nd=0; nd < Sub_mesh_pt[5]->nnode(); nd++)
      {
	Sub_mesh_pt[5]->node_pt(nd)->x(2) *= -1;
	Sub_mesh_pt[5]->node_pt(nd)->x(0) *= -1;
      }

    // Flip submesh 6 over the z axis, x axis and y axis:
    for(unsigned nd=0; nd < Sub_mesh_pt[6]->nnode(); nd++)
      {
	Sub_mesh_pt[6]->node_pt(nd)->x(0) *= -1;
	Sub_mesh_pt[6]->node_pt(nd)->x(1) *= -1;
	Sub_mesh_pt[6]->node_pt(nd)->x(2) *= -1;
      }

    // Flip submesh 7 over the z axis and y axis
    for(unsigned nd=0; nd < Sub_mesh_pt[7]->nnode(); nd++)
      {
	Sub_mesh_pt[7]->node_pt(nd)->x(1) *= -1;
	Sub_mesh_pt[7]->node_pt(nd)->x(2) *= -1;
      }


    // Now combine the sub-meshes into a single spherical mesh with a single
    // boundary over the surface of the sphere.
    for(unsigned i_mesh=0; i_mesh<8; i_mesh++)
      {
	// Add nodes to the overall mesh
	for(unsigned nd=0; nd < Sub_mesh_pt[i_mesh]->nnode(); nd++)
	  {
	    this->add_node_pt(Sub_mesh_pt[i_mesh]->node_pt(nd));
	  }

	// Add elements to the overall mesh
	for(unsigned e=0; e < Sub_mesh_pt[i_mesh]->nelement(); e++)
	  {
	    // std::cout << Sub_mesh_pt[i_mesh]->element_pt(e)->ndof() << std::endl;
	    this->add_element_pt(Sub_mesh_pt[i_mesh]->element_pt(e));
	  }

	// Add nodes on curved boundary (b = 0) to overall boundary
	for(unsigned nd=0; nd < Sub_mesh_pt[i_mesh]->nboundary_node(0); nd++)
	  {
	    this->add_boundary_node(0,Sub_mesh_pt[i_mesh]->node_pt(nd));
	  }

	// // Add elements on curved boundary  (b = 0) to overall boundary
	// for(unsigned e=0; e < Sub_mesh_pt[i_mesh]->nboundary_element(0); e++)
	//   {
	//     this->add_boundary_element(0,
	//   }
      }

    for(unsigned e=0; e < this->nelement(); e++)
      {
	std::cout << this->element_pt(e)->ndof() << std::endl;
      }

    // Setup boundary element lookup schemes
    setup_boundary_element_info();
  }


} // End of oomph namespace

#endif
