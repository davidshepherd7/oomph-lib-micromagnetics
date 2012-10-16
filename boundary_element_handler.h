#ifndef OOMPH_BOUNDARY_ELEMENT_HANDLER_H
#define OOMPH_BOUNDARY_ELEMENT_HANDLER_H

/*
  TODO:

  * Store integration scheme?

  * Parallelisation?

  * some way to check if the mesh has changed?

  * Fix for non-rectangle boundaries!!

  */

#include "generic.h"

using namespace oomph;

namespace oomph
{

 // =================================================================
 /// Simple class to store a list of angles associated with nodes of the
 /// boundary element mesh for assembly of the matrix.
 // =================================================================
 class CornerAngleList
 {
 public:
  /// Default constructor
  CornerAngleList() {}

  /// Destructor
  ~CornerAngleList() {}

  /// Set up corners for a rectangular/cubeoid shape
  void set_up_rectangular_corners(const Mesh* const mesh_pt)
  {
   unsigned nnode = mesh_pt->nnode();
   Corners.assign(nnode,0.0);

#ifdef PARANOID
   if(nnode == 0)
    {
     std::ostringstream error_msg;
     error_msg << "No nodes in mesh.";
     throw OomphLibError(error_msg.str(),
                         "CornerAngleList::set_up_rectangular_corners",
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif

   // Get dimension (from first node)
   unsigned dim = mesh_pt->node_pt(0)->ndim();

   // The angle is assumed to be 1/(2^dim), since it is rectangular or
   // cubeoid.
   double angle = 1 / std::pow(2,dim);

   // Loop over nodes
   for(unsigned nd=0; nd<nnode; nd++)
    {
     Node* nd_pt = mesh_pt->node_pt(nd);

     // Get list of boundaries that this node is on
     std::set<unsigned>* boundaries_pt;
     nd_pt->get_boundaries_pt(boundaries_pt);

#ifdef PARANOID
     if(boundaries_pt == 0)
      {
       std::ostringstream error_msg;
       error_msg << "Node is not on any boundaries, this probably "
                 << "means something has gone wrong, maybe you passed "
                 << "in the bulk mesh?";
       throw OomphLibError(error_msg.str(),
                           "CornerAngleList::set_up_rectangular_corners",
                           OOMPH_EXCEPTION_LOCATION);
      }
#endif

     // If it is on dim many boundaries then this is a corner, otherwise
     // assume it is a smooth point and so the angle is 0.5
     if(boundaries_pt->size() == dim)
      {
       Corners[nd] = angle;
      }
     else
      {
       // Angle = pi/2pi or 2pi/4pi in 2 or 3 dimensions respectively.
       Corners[nd] = 0.5;
      }
    }
  }

  /// Set up corners for a smooth mesh (i.e. no sharp corners).
  void set_up_smooth_mesh(const Mesh* const mesh_pt)
  {
   Corners.assign(mesh_pt->nnode(),0.5);
  }

  /// Add the contribution due to corners to the diagonal of the boundary
  /// matrix.
  void add_corner_contributions(DenseDoubleMatrix& bem_matrix) const
  {
   // Assume that the bem matrix is a densedoublematrix so that we can write
   // to it with operator().

#ifdef PARANOID
   // Check that the list has been set up
   if(!(is_set_up()))
    {
     std::ostringstream error_msg;
     error_msg << "Corner list has not been set up.";
     throw OomphLibError(error_msg.str(),
                         "CornerAngleList::add_corner_contributions",
                         OOMPH_EXCEPTION_LOCATION);
    }

   // Check that it is the correct size
   if(bem_matrix.nrow() != Corners.size())
    {
     std::ostringstream error_msg;
     error_msg << "Corners list is the wrong size for the matrix.";
     throw OomphLibError(error_msg.str(),
                         "",
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif

   // Add the fractional angles
   for(unsigned nd=0, s=Corners.size(); nd<s; nd++)
    {
     bem_matrix(nd,nd) += Corners[nd];
    }
  }

  /// Check if the list has been set up.
  bool is_set_up() const
  {
   // If there is a non-zero length vector something has been set up.
   return (Corners.size() != 0);
  }

 private:

  /// Storage for the location and angle of the corners.
  Vector<double> Corners;

  /// Inaccessible copy constructor
  CornerAngleList(const CornerAngleList& dummy)
  {BrokenCopy::broken_copy("CornerAngleList");}

  /// Inaccessible assignment operator
  void operator=(const CornerAngleList& dummy)
  {BrokenCopy::broken_assign("CornerAngleList");}
 };

 // =================================================================
 /// A class implementing all the boundary element methods stuff needed for
 /// the hybrid method in micromagnetics. Problems can then simply contain
 /// an instance of this class to gain access to hybrid method functions
 /// (see the "Prefer composition over inheritance" principle).
 ///
 /// In particular this allows the semi-implict "problem" class to not
 /// inherit from problem. This is good because many functions defined in
 /// Problem (for fully implicit methods) do not make sense for the
 /// semi-implicit "problem".
 // =================================================================
 template<class BEM_ELEMENT>
 class BoundaryElementHandler
 {

  // Note: in BEM input (output) index == jacobian col (row) == phi_1 (phi)
  // global equation numbers respectively.

 public:

  /// Default constructor
  BoundaryElementHandler() {}

  /// Destructor
  ~BoundaryElementHandler() {}

  /// Put the (output) values of the bem into a vector.
  void get_bem_values(DoubleVector& bem_values) const;

  /// Build the mesh, lookup schemes and matrix in that order.
  void build()
  {
   // Construct the mesh using on boundaries specified in Bem_boundaries
   build_bem_mesh();

   // Construct the lookup schemes
   Input_lookup.build(bem_mesh_pt(), input_index());
   Output_lookup.build(bem_mesh_pt(), output_index());

   // Set up the corners ??ds this is stupid!
   corner_list_pt()->set_up_rectangular_corners(bem_mesh_pt());

   // Construct the (dense) matrix
   build_boundary_matrix();
  }

  /// Use BEM on boundary b of the mesh.
  void set_bem_boundary(const unsigned& b, const Mesh* const mesh_pt)
  {
   std::pair<unsigned, const Mesh*> bound = std::make_pair(b,mesh_pt);
   Bem_boundaries.push_back(bound);
  }

  /// Use BEM on all boundaries in the mesh
  void set_bem_all_boundaries(const Mesh* mesh_pt)
  {
   for(unsigned b = 0; b < mesh_pt->nboundary(); b++)
    {
     set_bem_boundary(b, mesh_pt);
    }
  }

  // Access functions:
  // ============================================================

  /// Access to the pointer to the boundary element method mesh. Use
  /// pointer for consistency with everything else ever even though we
  /// store the actual value in the class.
  const Mesh* bem_mesh_pt() const {return &Bem_mesh;}

  /// Const access to the boundary matrix
  const DenseDoubleMatrix* boundary_matrix_pt() const
  {return &Boundary_matrix;}

  /// Non-const access to the boundary matrix
  DenseDoubleMatrix* boundary_matrix_pt() {return &Boundary_matrix;}

  /// \short Non-const access function for Input_index.
  unsigned& input_index() {return Input_index;}

  /// \short Const access function for Input_index.
  unsigned input_index() const {return Input_index;}

  /// \short Non-const access function for Output_index.
  unsigned& output_index() {return Output_index;}

  /// \short Const access function for Output_index.
  unsigned output_index() const {return Output_index;}

  /// \short Const access function for a pointer to the map of sharp corner
  /// angles at nodes.
  const CornerAngleList* corner_list_pt() const
  {return &Corner_list;}

  /// \short Non-const access function for a pointer to the map of sharp corner
  /// angles at nodes.
  CornerAngleList* corner_list_pt() {return &Corner_list;}

  /// \short Pointer to a lookup between output value global equation
  /// numbers and node numbers within mesh.
  const NodeGlobalNumbersLookup* output_lookup_pt() const
  {return &Output_lookup;}

  /// \short Alias of output_lookup_pt for when we are working with Jacobians
  /// (they are the same lookup).
  const NodeGlobalNumbersLookup* row_lookup_pt() const
  {return output_lookup_pt();}

  /// \short Pointer to a lookup between input value global equation
  /// numbers and node numbers within mesh.
  const NodeGlobalNumbersLookup* input_lookup_pt() const
  {return &Input_lookup;}

  /// \short Alias of input_lookup_pt for when we are working with Jacobians
  /// (they are the same lookup).
  const NodeGlobalNumbersLookup* col_lookup_pt() const
  {return input_lookup_pt();}

  /// \short Get an appropriate linear algebra distribution for working
  /// with the boundary matrix.
  void get_bm_distribution(LinearAlgebraDistribution& dist) const;

 private:

  /// \short Lookup between output value global equation numbers and node
  /// numbers within mesh.
  NodeGlobalNumbersLookup Output_lookup;

  /// \short Lookup between input value global equation numbers and node
  /// numbers within mesh.
  NodeGlobalNumbersLookup Input_lookup;

  /// The pointer to the "boundary element" mesh (as in boundary element method
  /// not finite elements on the boundary).
  Mesh Bem_mesh;

  /// A list of the boundaries (on various meshes) to which the boundary
  /// element method should be applied.
  Vector<std::pair<unsigned, const Mesh*> > Bem_boundaries;

  /// Pointer to storage for the list of nodal angles/solid angles.
  CornerAngleList Corner_list;

  /// The (local/elemental) index of the dof we take input values from.
  unsigned Input_index;

  /// The (local/elemental) index of the dof we are determining the
  /// boundary conditions for.
  unsigned Output_index;

  /// Matrix to store the relationship between phi_1 and phi on the boundary
  DenseDoubleMatrix Boundary_matrix;

  /// Construct BEM elements on boundaries listed in Bem_boundaries and add
  /// to the Bem_mesh.
  void build_bem_mesh();

  /// Construct the boundary matrix using the Bem_mesh.
  void build_boundary_matrix();

  /// \short Get the mapping between the global equation numbering and
  /// the boundary equation numbering.
  void create_global_boundary_equation_number_maps();

  /// Inaccessible copy constructor
  BoundaryElementHandler(const BoundaryElementHandler& dummy)
  {BrokenCopy::broken_copy("BoundaryElementHandler");}

  /// Inaccessible assignment operator
  void operator=(const BoundaryElementHandler& dummy)
  {BrokenCopy::broken_assign("BoundaryElementHandler");}

 };


 //==========================================================================
 /// Get the fully assembled boundary matrix in dense storage.
 //==========================================================================
 template<class BEM_ELEMENT>
 void BoundaryElementHandler<BEM_ELEMENT>::
 build_boundary_matrix()
 {

#ifdef PARANOID
  if(0)
   {
    std::ostringstream error_msg;
    error_msg << "";
    throw OomphLibError(error_msg.str(),
			"",
			OOMPH_EXCEPTION_LOCATION);
   }

  // Check the corner list has been set up
  if((corner_list_pt() == 0) || !(corner_list_pt()->is_set_up()))
   {
    std::ostringstream error_msg;
    error_msg << "List of the sharp corners of the mesh has not been set up.";
    throw OomphLibError(error_msg.str(),
			"BoundaryElementHandler::build_bem_mesh",
			OOMPH_EXCEPTION_LOCATION);
   }

#endif

  // Get the number of nodes in the boundary problem
  unsigned long n_node = bem_mesh_pt()->nnode();

  // Initialise and resize the boundary matrix
  Boundary_matrix.resize(n_node,n_node);
  Boundary_matrix.initialise(0.0);

  // Loop over all elements in the BEM mesh
  unsigned long n_bem_element = bem_mesh_pt()->nelement();
  for(unsigned long e=0;e<n_bem_element;e++)
   {
    // Get the pointer to the element (and cast to FiniteElement)
    BEM_ELEMENT* elem_pt = dynamic_cast<BEM_ELEMENT* >
     (bem_mesh_pt()->element_pt(e));

    // Find number of nodes in the element
    unsigned long n_element_node = elem_pt->nnode();

    // Set up and initialise matrix
    DenseMatrix<double> element_boundary_matrix(n_element_node,n_node,0.0);

    // Fill the matrix
    elem_pt->fill_in_contribution_to_boundary_matrix(element_boundary_matrix);

    // Loop over the nodes in this element (to copy results into final matrix)
    for(unsigned l=0;l<n_element_node;l++)
     {
      // Get the node number (in the bem mesh) from the global equation number.
      unsigned global_l_number = elem_pt->node_pt(l)->eqn_number(input_index());
      unsigned l_number = input_lookup_pt()->global_to_node(global_l_number);

      // Loop over all nodes in the mesh and add contributions from this element
      for(unsigned long s_nd=0; s_nd<n_node; s_nd++)
       {
	Boundary_matrix(l_number,s_nd) -= element_boundary_matrix(l,s_nd);
	// I think the sign here is negative because the lindholm formula
	// is for +ve but our final equation has negative kernel...
       }
     }
   }

  // Lindholm formula does not contain the solid angle contribution so add
  // it.
  corner_list_pt()->add_corner_contributions(Boundary_matrix);

 }


 //======================================================================
 /// Build the mesh of bem elements.
 //======================================================================
 template<class BEM_ELEMENT>
 void BoundaryElementHandler<BEM_ELEMENT>::
 build_bem_mesh()
 {
#ifdef PARANOID
  // Check list of BEM boundaries is not empty
  if(Bem_boundaries.size() == 0)
   {
    std::ostringstream error_msg;
    error_msg << "No BEM boundaries are set so there is no need"
	      << " to call build_bem_mesh().";
    throw OomphLibWarning(error_msg.str(),
			  "BoundaryElementHandler::build_bem_mesh",
			  OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // Create a set to temporarily store the list of boundary nodes (we use a
  // set because they automatically detect duplicates).
  std::set<Node*> node_set;

  // Loop over entries in Bem_boundaries vector.
  for(unsigned i=0; i < Bem_boundaries.size(); i++)
   {
    // Get mesh pointer and boundary number from vector.
    const unsigned b = Bem_boundaries[i].first;
    const Mesh* mesh_pt = Bem_boundaries[i].second;

    // Loop over the nodes on boundary b adding to the set of nodes.
    for(unsigned n=0, nnd=mesh_pt->nboundary_node(b); n<nnd;n++)
     {
      node_set.insert(mesh_pt->boundary_node_pt(b,n));
     }

    // Loop over the elements on boundary b creating bem elements
    for(unsigned e=0, ne=mesh_pt->nboundary_element(b); e<ne;e++)
     {
      // Create the corresponding BEM Element
      BEM_ELEMENT* bem_element_pt =
       new BEM_ELEMENT (mesh_pt->boundary_element_pt(b,e),
			mesh_pt->face_index_at_boundary(b,e));

      // Add the new BEM element to the BEM mesh
      Bem_mesh.add_element_pt(bem_element_pt);

      //??ds
      // // Set integration pointer
      // bem_element_pt->set_integration_scheme(bem_integration_scheme_pt());

      // Set the mesh pointer
      bem_element_pt->set_boundary_mesh_pt(&Bem_mesh);
     }
   }

  // Iterate over all nodes in the set and add them to the BEM mesh
  std::set<Node*>::iterator it;
  for(it=node_set.begin(); it!=node_set.end(); it++)
   {
    Bem_mesh.add_node_pt(*it);
   }

 }

 // =================================================================
 /// If the boundary matrix is distributed then get its
 /// distribution. Otherwise return a "non-distributed distribution" with
 /// the correct number of rows and a null comm pointer.
 // =================================================================
 template<class BEM_ELEMENT>
 void BoundaryElementHandler<BEM_ELEMENT>::
 get_bm_distribution(LinearAlgebraDistribution& dist) const
 {
  // Try to cast to a distributed object.
  const DistributableLinearAlgebraObject* dist_bm_pt =
   dynamic_cast<const DistributableLinearAlgebraObject* >
   (boundary_matrix_pt());

  // If it's not distributable (i.e. if the cast failed) then make a dummy
  // one, otherwise copy the distribution.
  if(dist_bm_pt == 0)
   {
    dist.build(0, boundary_matrix_pt()->nrow(), false);
   }
  else
   {
    dist.build(dist_bm_pt->distribution_pt());
   }
 }


 // =================================================================
 /// Put the output values from the boundary element method into a
 /// vector.
 // =================================================================
 template<class BEM_ELEMENT>
 void BoundaryElementHandler<BEM_ELEMENT>::
 get_bem_values(DoubleVector& bem_output_values) const
 {
  // Get the boundary matrix linear algebra distribution (if there is one).
  LinearAlgebraDistribution dist;
  get_bm_distribution(dist);

  // Set up double vectors
  DoubleVector input_values(dist);
  bem_output_values.build(dist);

  // Get input values
  for(unsigned nd=0, nnode=bem_mesh_pt()->nnode(); nd<nnode; nd++)
   {
    input_values[nd] = bem_mesh_pt()->node_pt(nd)->value(input_index());
   }

  // Matrix multiply to get output values
  boundary_matrix_pt()->multiply(input_values, bem_output_values);
 }


} // End of oomph namespace

#endif
