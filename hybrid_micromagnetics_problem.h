#ifndef OOMPH_MICROMAGNETICS_PROBLEM_H
#define OOMPH_MICROMAGNETICS_PROBLEM_H


/*
description of file goes here
*/

#include "generic.h"
#include "./magnetic_materials.h"

using namespace oomph;

namespace oomph
{

  template<class BULK_ELEMENT,
	   template<class BULK_ELEMENT,unsigned DIM> class BEM_ELEMENT,
	   unsigned DIM>
  class HybridMicromagneticsProblem : public Problem
  {

  public:

    HybridMicromagneticsProblem()
      : Flux_mesh_pt(0), Bem_mesh_pt(0),
	Bulk_mesh_pt(0), Bulk_mesh_magnetic_parameters_pt(0)
    {};

    virtual void doc_solution(DocInfo& doc_info) = 0;

    virtual void set_initial_condition() = 0;

    /// Finish off construction of problem (once mesh and magnetic properties
    /// have been set up).
    void finish_building_hybrid_problem();

    void build_bem_mesh(Mesh* bem_mesh_pt) const;

    void build_boundary_matrix();

    void create_flux_elements(const unsigned& b, Mesh* const &bulk_mesh_pt,
			      Mesh* const &surface_mesh_pt);

    /// \short Get the mapping between the global equation numbering and
    /// the boundary equation numbering.
    void create_global_boundary_equation_number_maps();

    /// Get the boundary equation number from the global equation number
    unsigned get_boundary_equation_number(const Node* const boundary_node) const;

    /// Get a doublevector of the values of phi_1 on the boundary, ordered by
    /// boundary equation number.
    void get_boundary_phi_1(DoubleVector& boundary_phi_1) const;

    /// Update the values of phi on the boundary
    void get_bem_phi_values(DoubleVector& bem_phi_values) const;

    /// Access to the pointer to the boundary element method mesh
    Mesh* bem_mesh_pt() const {return Bem_mesh_pt;}

    Mesh* bulk_mesh_pt() const {return Bulk_mesh_pt;}
    Mesh*& bulk_mesh_pt() {return Bulk_mesh_pt;}

    MagneticParameters* magnetic_parameters_pt() const
    {return Bulk_mesh_magnetic_parameters_pt;}

    MagneticParameters*& magnetic_parameters_pt()
    {return Bulk_mesh_magnetic_parameters_pt;}

    /// Get the index of phi for use in BEM mapping
    unsigned phi_index() const {return Phi_index;}

    /// Get the index of phi_1 for use in BEM mapping
    unsigned phi_1_index() const {return Phi_1_index;}

    /// Get the flux mesh pointer.
    Mesh* flux_mesh_pt() const {return Flux_mesh_pt;}

    /// Access to the boundary matrix
    DenseDoubleMatrix* boundary_matrix_pt()
    {return &Boundary_matrix;}

    void get_residuals(DoubleVector& residuals)
    {
      //std::cout << "Calling your get_residuals" << std::endl;
      Problem::get_residuals(residuals);
      insert_bem_phi_residual_contribution(residuals);
    }

    void insert_bem_phi_residual_contribution(DoubleVector& residuals) const
    {
      // calculate values for phi from bem
      DoubleVector bem_phi_values;
      get_bem_phi_values(bem_phi_values);

      // for each node in bem mesh
      for(unsigned nd=0; nd< bem_mesh_pt()->nnode(); nd++)
	{
	  // Get a pointer
	  Node* nd_pt = bem_mesh_pt()->node_pt(nd);

	  // get the the bem equation number
	  unsigned bem_eqn_num = get_boundary_equation_number(nd_pt);
	  unsigned global_eqn_num = nd_pt->eqn_number(phi_index());

	  // insert appropriate value into residuals
	  double r = bem_phi_values[bem_eqn_num]
	    - nd_pt->value(phi_index());

	  residuals[global_eqn_num] = r;
	}
    }

    /// Overload get_jacobian to include the boundary matrix in a sparse form.
    void get_jacobian(DoubleVector& residuals, CRDoubleMatrix& jacobian)
    {
      std::cout << "Calling your get_jacobian function with bem added directly." << std::endl;

      // Get the fem jacobian (in the same distribution pattern as original).
      LinearAlgebraDistribution* dist_pt = jacobian.distribution_pt();
      CRDoubleMatrix sparse_jacobian(dist_pt);
      Problem::get_jacobian(residuals,sparse_jacobian);

      // Finish off the residual calculation
      insert_bem_phi_residual_contribution(residuals);

       // Create a sum of matrices holding the total jacobian
      SumOfMatrices lazy_sum;
      lazy_sum.main_matrix_pt() = &sparse_jacobian;

      lazy_sum.add_matrix(boundary_matrix_pt(),
			  &Global_boundary_equation_num_map,
			  &Global_phi_1_num_map,0);

      // Add an identity matrix * -1: phi dependence on itself,
      // delete when done with Jacobain.
      unsigned n_bem_nodes = bem_mesh_pt()->nnode();
      Vector<int> id_row_indices(n_bem_nodes), id_col_indices(n_bem_nodes);
      Vector<double> id_values(n_bem_nodes,-1.0);
      for(int i=0; i<int(n_bem_nodes); i++) {id_row_indices[i] = i; id_col_indices[i] = i;}

      CRDoubleMatrix* neg_id_matrix_pt =
	new CRDoubleMatrix(id_row_indices,id_col_indices,id_values,n_bem_nodes,n_bem_nodes);

      neg_id_matrix_pt->Matrix::sparse_indexed_output("id_matrix");

      lazy_sum.add_matrix(neg_id_matrix_pt,
			  &Global_boundary_equation_num_map,
			  &Global_boundary_equation_num_map,
			  1);

      // Build a new sparse jacobian from data outputted from the sum
      Vector<int> rows, cols;
      Vector<double> values;
      lazy_sum.get_as_indicies(rows,cols,values);
      jacobian.build(dist_pt,rows,cols,values,
    		     lazy_sum.nrow(),lazy_sum.ncol());

      // dump jacobian for tests
      //why are the Matrix:: needed? - because inheritence of crdouble matrix is weird!
      sparse_jacobian.Matrix::sparse_indexed_output(std::string("matrices/sparse_jac_part"));
      jacobian.Matrix::sparse_indexed_output("matrices/cr_jacobian");
      residuals.output("matrices/residual");
    }

  private:

    /// The map between the global equation numbers for phi and the boundary
    /// equation/matrix numbering.
    std::map<long unsigned,long unsigned> Global_boundary_equation_num_map;

    /// The map between the global equation numbers for phi_1 and the boundary
    /// equation/matrix numbering.
    std::map<long unsigned,long unsigned> Global_phi_1_num_map;

    Mesh* Flux_mesh_pt;

    /// The pointer to the boundary element method mesh
    Mesh* Bem_mesh_pt;

    ///
    Mesh* Bulk_mesh_pt;

    /// Store the magnetic parameters of the bulk mesh region.
    //??ds eventually may need to be able to store multiple meshes with
    // different parameters
    MagneticParameters* Bulk_mesh_magnetic_parameters_pt;

    /// Store the index of phi for use in BEM mapping
    unsigned Phi_index;

    /// Store the index of phi_1 for use in BEM mapping
    unsigned Phi_1_index;

    /// Matrix to store the relationship between phi_1 and phi on the boundary
    DenseDoubleMatrix Boundary_matrix;

  };

  // ============================================================
  /// Finish off construction of problem (once mesh has been set up)
  // ============================================================
  template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
  void HybridMicromagneticsProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
  finish_building_hybrid_problem()
  {
#ifdef PARANOID
    if(Bulk_mesh_magnetic_parameters_pt == 0)
      {
	std::ostringstream error_msg;
	error_msg << "Magnetic parameters of bulk mesh not set up.";
	throw OomphLibError(error_msg.str(),
			    "HybridMicromagneticsProblem::build()",
			    OOMPH_EXCEPTION_LOCATION);
      }

    if(Bulk_mesh_pt == 0)
      {
	std::ostringstream error_msg;
	error_msg << "Bulk mesh pointer not set up.";
	throw OomphLibError(error_msg.str(),
			    "HybridMicromagneticsProblem::build()",
			    OOMPH_EXCEPTION_LOCATION);
      }
#endif

    // Finish bulk elements
    // ============================================================
    // Loop over elements in bulk mesh to set function pointers
    for(unsigned i=0; i< bulk_mesh_pt()->nelement(); i++)
      {
	// Upcast from GeneralisedElement to the present element
	BULK_ELEMENT* elem_pt = dynamic_cast<BULK_ELEMENT*>(bulk_mesh_pt()->element_pt(i));

	// Set pointer to continuous time
	elem_pt->time_pt() = time_pt();

	// Set the function pointers and values for magnetic parameters
	elem_pt->magnetic_parameters_pt() = magnetic_parameters_pt();
      }

    // Get the indicies for phi and phi_1 via casting a pointer
    BULK_ELEMENT* elem_pt = dynamic_cast< BULK_ELEMENT* >(bulk_mesh_pt()->element_pt(0));
    Phi_index = elem_pt->phi_index_micromag();
    Phi_1_index = elem_pt->phi_1_index_micromag();


    // Pin the value of phi_1 to zero at a point close to the origin.
    // ??ds this is likely to go wrong if the origin is ont he boundary, fix?!
    // ============================================================
    bool found = false;
    for(unsigned nd=0; nd< bulk_mesh_pt()->nnode(); nd++)
      {
	Vector<double> nd_x(DIM,0.0);
	for(unsigned j=0; j<DIM; j++)
	  nd_x[j] = bulk_mesh_pt()->node_pt(nd)->x(j);
	if ( small(nd_x[0]) && small(nd_x[1]) && small(nd_x[2]) )
	  {
	    bulk_mesh_pt()->node_pt(nd)->pin(Phi_1_index);
	    bulk_mesh_pt()->node_pt(nd)->set_value(Phi_1_index,0.0);
	    found = true;
	    break;
	  }
      }
    if (!(found))
      throw OomphLibError("No node near middle","",OOMPH_EXCEPTION_LOCATION);

    // Create flux elements
    // ============================================================

    // We want Neumann (flux) boundary condtions on phi_1 on all boundaries so
    // create the face elements needed.
    Flux_mesh_pt = new Mesh;
    for(unsigned b=0; b < bulk_mesh_pt()->nboundary(); b++)
      {
	create_flux_elements(b,bulk_mesh_pt(),Flux_mesh_pt);
      }

    // Create BEM elements
    // ============================================================

    // Create integration scheme in case it is needed.
    QVariableOrderGaussLegendre<DIM-1>* variable_scheme_pt
      = new QVariableOrderGaussLegendre<DIM-1>; //??ds generalise to triangles??

    // Create BEM elements on all boundaries and add to BEM mesh
    Bem_mesh_pt = new Mesh;
    build_bem_mesh(bem_mesh_pt());

    // Set boundary mesh pointer in element (static memeber so only do once)
    BEM_ELEMENT<BULK_ELEMENT,DIM>::set_boundary_mesh_pt(bem_mesh_pt());

    // Set pointers in elements
    for(unsigned i_ele = 0; i_ele < bem_mesh_pt()->nelement(); i_ele++)
      {
	// Cast pointer
	BEM_ELEMENT<BULK_ELEMENT,DIM>* ele_pt
	  = dynamic_cast< BEM_ELEMENT< BULK_ELEMENT,DIM>* >
	  (bem_mesh_pt()->element_pt(i_ele));

	// Set integration scheme
	ele_pt->set_integration_scheme(variable_scheme_pt);
      }

    // Build global finite element mesh
    add_sub_mesh(bulk_mesh_pt());
    add_sub_mesh(Flux_mesh_pt);
    build_global_mesh();

    // Setup equation numbering scheme for all the finite elements
    std::cout << "FEM number of equations: " << assign_eqn_numbers() << std::endl;

    // Make the boundary matrix (including setting up the numbering scheme).  Note
    // that this requires the FEM numbering scheme to be already set up.
    build_boundary_matrix();

  }

  //=============================================================================
  /// Get the fully assembled boundary matrix in dense storage.
  //=============================================================================
  template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
  void HybridMicromagneticsProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
  build_boundary_matrix()
  {
    // Create the mapping from global to boundary equations
    create_global_boundary_equation_number_maps();

    //std::cout << Global_boundary_equation_num_map << std::endl;

    // get the number of nodes in the boundary problem
    unsigned long n_node = bem_mesh_pt()->nnode();

    // Initialise and resize the boundary matrix
    Boundary_matrix.resize(n_node,n_node);
    Boundary_matrix.initialise(0.0);

    // Loop over all elements in the BEM mesh
    unsigned long n_bem_element = bem_mesh_pt()->nelement();
    for(unsigned long e=0;e<n_bem_element;e++)
      {
	// Get the pointer to the element (and cast to FiniteElement)
	BEM_ELEMENT<BULK_ELEMENT,DIM>* elem_pt =
	  dynamic_cast < BEM_ELEMENT<BULK_ELEMENT,DIM>* > (bem_mesh_pt()->element_pt(e));

	// Find number of nodes in the element
	unsigned long n_element_node = elem_pt->nnode();

	// Set up and initialise matrix
	DenseMatrix<double> element_boundary_matrix(n_element_node,n_node,0.0);

	// Fill the matrix
	elem_pt->fill_in_contribution_to_boundary_matrix(element_boundary_matrix);

	// Loop over the nodes in this element (to copy results into final matrix)
	for(unsigned l=0;l<n_element_node;l++)
	  {
	    // Get the boundary equation (=node) number from the global one
	    unsigned l_number =
	      this->get_boundary_equation_number
	      (elem_pt->node_pt(l));

	    // Loop over all nodes in the mesh and add contributions from this element
	    for(unsigned long source_node=0; source_node<n_node; source_node++)
	      {
		unsigned source_number =
		  this->get_boundary_equation_number
		  (bem_mesh_pt()->node_pt(source_node));

		Boundary_matrix(l_number,source_number)
		  -= element_boundary_matrix(l,source_node);
		//??Ds I think the sign here is negative because the lindholm
		//formula is for +ve but our final equation has negative
		//kernel...
	      }
	  }
      }

// #ifdef PARANOID
//     // If the mesh is not one I've dealt with here:
//     if((dynamic_cast<RectangularQuadMesh<BULK_ELEMENT>*>(bulk_mesh_pt()) == 0)
//        && (dynamic_cast<SimpleCubicMesh<BULK_ELEMENT>*>(bulk_mesh_pt()) == 0)
//        && (dynamic_cast<TetgenMesh<BULK_ELEMENT>*>(bulk_mesh_pt()) == 0))
//       {
// 	std::ostringstream error_msg;
// 	error_msg << "No corner data for this mesh.";
// 	throw OomphLibError(error_msg.str(),
// 			    "HybridMicromagneticsProblem::build_boundary_matrix()",
// 			    OOMPH_EXCEPTION_LOCATION);
//       }
// #endif

    // Lindholm formula does not contain the solid angle contribution so add it
    // here. Loop over the matrix diagonals adding the angle factor.
    for(unsigned long nd = 0; nd < bem_mesh_pt()->nnode(); nd++)
      {
	// // We know that for rectangles/cubeoids corner nodes are on DIM many
	// // boundaries so check this to find the corners.
	// unsigned n_bound = 0;
	// if(bem_mesh_pt()->node_pt(nd)->is_on_boundary())
	//   {
	//     n_bound = bem_mesh_pt()->node_pt(nd)->get_boundaries_pt()->size();
	//   }
	// //??ds fix this sometime...
	// if(n_bound == DIM)
	//   {
	//     // Get appropriate angle for rectangle/cubeoid
	//     double angle = ((DIM == 2) ? 0.25 : 0.125);
	//     Boundary_matrix(nd,nd) += angle;
	//   }
	// else if(n_bound > DIM)
	//   throw OomphLibError("Something has gone wrong here...",
	// 		      "HybridMicromagneticsProblem::build_boundary_matrix()",
	// 		      OOMPH_EXCEPTION_LOCATION);
	// else
	{
	  // This only accounts for points which are smooth (at least in the limit
	  // of infinite refinement) so the solid angle contribution is
	  // (2*pi)/(4*pi) = 0.5.
	  Boundary_matrix(nd,nd) += 0.5;
	}
      }
  }

  //======================================================================
  /// Create the map between the global equation numbering system and the
  /// boundary equation numbering system.
  //======================================================================
  template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
  void HybridMicromagneticsProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
  create_global_boundary_equation_number_maps()
  {
    // Initialise the map
    Global_boundary_equation_num_map.clear();
    Global_phi_1_num_map.clear();

    // Initialise counters for number of unique boundary nodes included so far
    unsigned k_phi = 0, k_phi_1 = 0;

    // Loop over boundary nodes assigning a boundary equation number to each.
    unsigned n_boundary_node = this->bem_mesh_pt()->nnode();
    for(unsigned i_node=0; i_node<n_boundary_node; i_node++)
      {
	// Get global equation number for phi
	unsigned global_phi_number = this->bem_mesh_pt()->
	  node_pt(i_node)->eqn_number(phi_index());

	// Get global equation number for phi_1
	unsigned global_phi_1_number = this->bem_mesh_pt()->
	  node_pt(i_node)->eqn_number(phi_1_index());

	// Set up the pair ready to input with key="global equation number" and
	// value ="boundary equation number"=k.
	std::pair<unsigned,unsigned> input_pair_phi
	  = std::make_pair(global_phi_number,k_phi);
	std::pair<unsigned,unsigned> input_pair_phi_1
	  = std::make_pair(global_phi_1_number,k_phi_1);

	// Add entry to map and store whether this was a new addition
	bool new_addition_phi = (Global_boundary_equation_num_map.insert(input_pair_phi)
				 ).second;
	bool new_addition_phi_1 = (Global_phi_1_num_map.insert(input_pair_phi_1)
				   ).second;

	// Increment k if this was a new addition to the map
	if(new_addition_phi) k_phi++;
	if(new_addition_phi_1) k_phi_1++;
      }

    //std::cout << Global_boundary_equation_num_map << std::endl;
    // std::cout << Global_phi_1_num_map << std::endl;
  }

  //======================================================================
  /// Given the global equation number return the boundary equation number. Most
  /// of the function is an error check.
  //======================================================================
  template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
  unsigned HybridMicromagneticsProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
  get_boundary_equation_number(const Node* const boundary_node_pt) const
  {
    // Get the global equation number for phi for this node (we could use any
    // equation at this node).
    long global_equation_number = boundary_node_pt->eqn_number(phi_index());

#ifdef PARANOID
    // If the iterator is placed at the end the given global equation number is
    // not in the map, so return an error.
    if(Global_boundary_equation_num_map.find(global_equation_number)
       == Global_boundary_equation_num_map.end())
      {
	std::ostringstream error_stream;
	error_stream << "Global equation number " << global_equation_number
		     << " is not in the global to boundary map.";
	throw OomphLibError(error_stream.str(),
			    "HybridMicromagneticsProblem::get_boundary_equation_number",
			    OOMPH_EXCEPTION_LOCATION);
      }

    //Problems could occur if the index we are using ever has pinned values...
    if(global_equation_number < 0)
      {
	//std::cout << Global_boundary_equation_num_map << std::endl;
	throw OomphLibError("Pinned equation, use a different eq num?",
			    "HybridMicromagneticsProblem::get_boundary_equation_number",
			    OOMPH_EXCEPTION_LOCATION);
      }
#endif

    return Global_boundary_equation_num_map.find(global_equation_number)->second;
  }


  //======================================================================
  /// Build the mesh of bem elements.
  //======================================================================
  template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
  void HybridMicromagneticsProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
  build_bem_mesh(Mesh* bem_mesh_pt) const
  {
    // Create a set to temporarily store the list of boundary nodes
    // (we do it via a set because sets automatically detect duplicates)
    std::set<Node*> node_set;
    std::set<Node*>::iterator it, c_it;

    // Loop over the boundaries
    unsigned n_boundary = bulk_mesh_pt()->nboundary();
    for(unsigned b=0; b<n_boundary; b++)
      {
	// Loop over the boundary nodes on boundary b making a set of nodes
	unsigned n_bound_node = bulk_mesh_pt()->nboundary_node(b);
	for(unsigned n=0;n<n_bound_node;n++)
	  node_set.insert(bulk_mesh_pt()->boundary_node_pt(b,n));

	// Loop over the elements on boundary b creating bem elements
	unsigned n_bound_element = bulk_mesh_pt()->nboundary_element(b);
	for(unsigned e=0;e<n_bound_element;e++)
	  {
	    // Create the corresponding BEM Element
	    BEM_ELEMENT<BULK_ELEMENT,DIM>* bem_element_pt = new BEM_ELEMENT<BULK_ELEMENT,DIM>
	      (bulk_mesh_pt()->boundary_element_pt(b,e),
	       bulk_mesh_pt()->face_index_at_boundary(b,e));

	    // Add the BEM element to the BEM mesh
	    bem_mesh_pt->add_element_pt(bem_element_pt);
	  }
      }

    // Iterate over all nodes in the set and add them to the BEM mesh
    for(it=node_set.begin(); it!=node_set.end(); it++)
      bem_mesh_pt->add_node_pt(*it);
  }

  //======================================================================
  /// Create potential flux boundary condition elements on boundary b.
  //======================================================================
  template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
  void HybridMicromagneticsProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
  create_flux_elements(const unsigned& b, Mesh* const &bulk_mesh_pt,
		       Mesh* const &surface_mesh_pt)
  {
    // How many bulk elements are adjacent to boundary b?
    unsigned n_element = bulk_mesh_pt->nboundary_element(b);

    // Loop over the bulk elements adjacent to boundary b
    for(unsigned e=0;e<n_element;e++)
      {
	// Get pointer to the bulk element that is adjacent to boundary b
	BULK_ELEMENT* bulk_elem_pt = dynamic_cast<BULK_ELEMENT*>
	  (bulk_mesh_pt->boundary_element_pt(b,e));

	// What is the index of the face of the bulk element at the boundary
	int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

	// Build the corresponding prescribed-flux element
	MicromagFluxElement<BULK_ELEMENT>* flux_element_pt =
	  new MicromagFluxElement<BULK_ELEMENT>(bulk_elem_pt,face_index);

	// Pass a pointer to the flux element to the bulk element
	bulk_elem_pt->add_face_element_pt(flux_element_pt);

	// Add the prescribed-flux element to the mesh
	surface_mesh_pt->add_element_pt(flux_element_pt);

      } // End of loop over bulk elements adjacent to boundary b
  }

  template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
  void HybridMicromagneticsProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
  get_boundary_phi_1(DoubleVector& boundary_phi_1) const
  {
    // Initialise vector
    LinearAlgebraDistribution dummy(0,bem_mesh_pt()->nnode(),false);
    boundary_phi_1.build(dummy,0.0);

    for(unsigned i_nd=0; i_nd< bem_mesh_pt()->nnode(); i_nd++)
      {
	Node* node_pt = bem_mesh_pt()->node_pt(i_nd);

	// Get the boundary equation number
	unsigned bdry_eqn_num = get_boundary_equation_number(node_pt);

	// Fill in the value
	boundary_phi_1[bdry_eqn_num] = node_pt->value(phi_1_index());
      }
  }

  //======================================================================
  /// Calculate new phi values using the BEM from multiplying boundary element
  /// matrix and phi_1.
  //======================================================================
  template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
  void HybridMicromagneticsProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
  get_bem_phi_values(DoubleVector& bem_phi_values) const
  {
    // In order to use DoubleVector (for matrix multiplications) we need to have
    // this thingy. In our serial case it just gives a number of rows.
    LinearAlgebraDistribution dummy(0,bem_mesh_pt()->nnode(),false);

    // Assemble a vector of phi_1 values on boundary nodes
    DoubleVector boundary_phi_1;
    get_boundary_phi_1(boundary_phi_1);

    // Dense matrix multiplication to calculate phi (result goes in Bem_phi_values)
    Boundary_matrix.multiply(boundary_phi_1, bem_phi_values);
  }


} // End of oomph namespace

#endif
