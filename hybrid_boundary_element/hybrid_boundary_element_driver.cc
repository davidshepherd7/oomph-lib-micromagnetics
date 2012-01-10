#ifndef OOMPH_HYBRID_BOUNDARY_ELEMENT_DRIVER_H
#define OOMPH_HYBRID_BOUNDARY_ELEMENT_DRIVER_H

/*
  ??ds  description of file goes here
*/

#include "generic.h"
#include "../micromagnetics_boundary_element.h"
#include "meshes/rectangular_quadmesh.h"

using namespace oomph;
using namespace MathematicalConstants;

namespace oomph
{

  template<class ELEMENT>
  class TwoDMicromagProblem : public Problem
  {

  public:

    /// Constructor
    TwoDMicromagProblem(const unsigned& n_x, const unsigned& n_y);

    /// Destructor (empty -- all the cleanup is done in the base class)
    ~TwoDMicromagProblem(){};

    /// Doc the solution
    void doc_solution(DocInfo& doc_info);

    /// Get number of bulk elements
    unsigned get_n_bulk_element() {return N_bulk_element;}

  private:

    /// Pointer to control node at which the solution is documented
    Node* Control_node_pt;

    //  /// Doc info object
    //  DocInfo Doc_info;

    //  /// Trace file
    //  std::ofstream Trace_file;

    /// Update the problem specs before solve
    // Nothing to do here since no dirichlet boundaries.
    void actions_before_newton_solve(){};

    /// Update the problem specs after solve
    void actions_after_newton_solve(){};

    /// Update the problem specs after solve (empty)
    void actions_after_implicit_timestep() {}

    /// Update the problem specs before next timestep
    void actions_before_implicit_timestep(){};

    /// Set initial condition (incl previous timesteps) according to specified function.
    void set_initial_condition(){};

    /// Number of bulk elements
    unsigned N_bulk_element;

  }; // end of problem class

  //=====start_of_constructor===============================================
  ///
  //========================================================================
  template<class ELEMENT>
  TwoDMicromagProblem<ELEMENT>::
  TwoDMicromagProblem(const unsigned& n_x, const unsigned& n_y)
  {
    // Set domain size
    double l_x = 3.0;
    double l_y = 1.0;

    // Allocate steady state timestepper
    add_time_stepper_pt(new Steady<2>);

    // Build mesh
    Problem::mesh_pt() = new RectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y,time_stepper_pt());

    // Store number of bulk elements
    N_bulk_element = mesh_pt()->nelement();

    // Choose a control node at which the solution is documented
    unsigned control_el = unsigned(N_bulk_element/2); // Pick a control element in the middle
    Control_node_pt=mesh_pt()->finite_element_pt(control_el)->node_pt(0);  // Choose its first node as the control node
    std::cout << "Recording trace of the solution at: " << Control_node_pt->x(0) << std::endl;

    // Boundary conditions?
    // No Dirichlet conditions so nothing to do here?
    // Apply conditions on phi_1, phi_2 here?

    // Loop over elements to set pointers
    // ??ds add this when we are actually solving things

    // Setup equation numbering scheme
    std::cout << "Number of equations: " << assign_eqn_numbers() << std::endl;

  } // end of constructor

  //=====================start_of_doc=======================================
  /// ??ds write this!
  //========================================================================
  template<class ELEMENT>
  void TwoDMicromagProblem<ELEMENT>::
  doc_solution(DocInfo& doc_info)
  {

    std::ofstream some_file;
    char filename[100];

    // Number of plot points
    unsigned npts;
    npts=5;

    // Output solution
    sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
	    doc_info.number());
    some_file.open(filename);
    mesh_pt()->output(some_file,npts);
    some_file.close();

    // // Output exact solution
    // //----------------------
    // sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
    // 	    doc_info.number());
    // some_file.open(filename);
    // mesh_pt()->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u);
    // some_file.close();

  } // end of doc

  template<class BULK_ELEMENT, template<class> class FACE_ELEMENT>
  class TwoDBoundaryProblem : public Problem
  {
  public:

    /// Constructor
    TwoDBoundaryProblem(Mesh* bulk_mesh_pt);

    /// Destructor
    ~TwoDBoundaryProblem(){};

    /// Bulk mesh pointer access function
    Mesh* bulk_mesh_pt() const {return Bulk_mesh_pt;}

    /// Build the mesh for this problem
    void build_face_mesh(Mesh* const &bulk_mesh_pt);

    /// Get the boundary element matrix (similar to problem::get_jacobian)
    void get_boundary_matrix(DenseDoubleMatrix& boundary_matrix);

  private:

    /// Create face elements on the i'th boundary and add to boundary mesh
    void create_face_elements(const unsigned &i_bulk_boundary);

    /// Pointer to the bulk mesh
    Mesh* Bulk_mesh_pt;
  };

  //=====start_of_constructor===============================================
  ///
  //========================================================================
  template<class BULK_ELEMENT, template<class> class FACE_ELEMENT>
  TwoDBoundaryProblem<BULK_ELEMENT,FACE_ELEMENT>::
  TwoDBoundaryProblem(Mesh* input_bulk_mesh_pt)
  {
    // Set up bulk mesh pointer
    Bulk_mesh_pt = input_bulk_mesh_pt;

    // Create (empty) boundary mesh
    mesh_pt() = new Mesh();

    // control node?

    // Create face elements on all boundaries and add them to the face mesh
    build_face_mesh(bulk_mesh_pt());

    // Pin the values of phi_1 on all nodes in the boundary mesh
    unsigned n_node = mesh_pt()->nnode();
    for(unsigned n=0; n<n_node; n++)
      {
    	mesh_pt()->node_pt(n)->pin(0); // assuming phi_1 is in storage 0
      }

    // Loop over elements to set mesh pointer
    unsigned n_element = mesh_pt()->nelement();
    for(unsigned i=0;i<n_element;i++)
      {
	// Upcast from GeneralisedElement to the face element
	FACE_ELEMENT<BULK_ELEMENT>* elem_pt =
	  dynamic_cast<FACE_ELEMENT<BULK_ELEMENT> *>(mesh_pt()->element_pt(i));

	// Set boundary mesh pointer in element
	elem_pt->set_mesh_pt(Problem::mesh_pt());
      }

    // Set up numbering scheme
    std::cout << "Number of equations: " << assign_eqn_numbers() << std::endl;

  } // End of constructor

  /// \short Constuct a Mesh of FACE_ELEMENTs along the b-th boundary
  /// of the mesh (which contains elements of type BULK_ELEMENT)
  template<class BULK_ELEMENT, template<class> class FACE_ELEMENT>
  void TwoDBoundaryProblem<BULK_ELEMENT,FACE_ELEMENT>::
  build_face_mesh(Mesh* const &bulk_mesh_pt)
  {
    // Create a set to temporarily store the list of boundary nodes
    // (we do it via a set because sets automatically detect duplicates)
    std::set<Node*> node_set;
    std::set<Node*>::iterator it;

    // Loop over the boundaries
    unsigned n_bulk_boundary = bulk_mesh_pt->nboundary();
    for(unsigned b=0; b<n_bulk_boundary; b++)
      {
	//Find the number of nodes on the boundary
	unsigned n_bound_node = bulk_mesh_pt->nboundary_node(b);

	//Loop over the boundary nodes and add them to the set
	for(unsigned n=0;n<n_bound_node;n++)
	  node_set.insert(bulk_mesh_pt->boundary_node_pt(b,n));


	//Find the number of elements next to the boundary
	unsigned n_bound_element = bulk_mesh_pt->nboundary_element(b);

	//Loop over the elements adjacent to boundary b
	for(unsigned e=0;e<n_bound_element;e++)
	  {
	    //Create the FaceElement
	    FACE_ELEMENT<BULK_ELEMENT>* face_element_pt =
	      new FACE_ELEMENT<BULK_ELEMENT>
	      (bulk_mesh_pt->boundary_element_pt(b,e),
	       bulk_mesh_pt->face_index_at_boundary(b,e));

	    //Add the face element to the face mesh
	    mesh_pt()->add_element_pt(face_element_pt);
	  }
      }

    // Iterate over all elements of the set and add to the mesh
    for(it=node_set.begin(); it!=node_set.end(); it++)
      mesh_pt()->add_node_pt(*it);

    //??ds taken from mesh.h - no idea what this does but maybe useful later...
#ifdef OOMPH_HAS_MPI
    // If the bulk mesh has been distributed then the face mesh is too
    if (this->is_mesh_distributed())
      {
	face_mesh_pt->set_mesh_distributed();
      }
#endif
  }


  //=============================================================================
  /// Get the fully assembled boundary matrix in dense storage.
  //=============================================================================
  template<class BULK_ELEMENT, template<class> class FACE_ELEMENT>
  void TwoDBoundaryProblem<BULK_ELEMENT,FACE_ELEMENT>::
  get_boundary_matrix(DenseDoubleMatrix& boundary_matrix)
  {

    // get the number of nodes in the boundary problem
    unsigned long n_node=mesh_pt()->nnode();

    // resize the boundary matrix
    boundary_matrix.resize(n_node,n_node);
    boundary_matrix.initialise(0.0);

    // Loop over all the elements
    unsigned long n_element = mesh_pt()->nelement();
    for(unsigned long e=0;e<n_element;e++)
      {
	// Get the pointer to the element (and cast to FiniteElement)
	FiniteElement* elem_pt =
	  dynamic_cast<FiniteElement*>(mesh_pt()->element_pt(e));

	// Find number of nodes in the element
	unsigned long n_element_node = elem_pt->nnode();

	// Set up a matrix and dummy residual vector
	DenseMatrix<double> element_boundary_matrix(2,n_node);
	Vector<double> dummy(0);

	// Fill the matrix
	assembly_handler_pt()
	  ->get_jacobian(elem_pt,dummy,element_boundary_matrix);

	// Loop over the nodes in this element
	for(unsigned l=0;l<n_element_node;l++)
	  {
	    //??ds how do I get the global node number?
	    //??ds is there a global node number - only via mesh, no good here :(
	    unsigned source_node = elem_pt->eqn_number(l);

	    std::cout << "source node = " << source_node << std::endl;

	    // Loop over all nodes in the mesh and add contributions from this element
	    for(unsigned long target_node=0;target_node<n_node;target_node++)
	      {
		boundary_matrix(source_node,target_node)
		  += element_boundary_matrix(l,target_node);
	      }
	  }
      }
  }

} // End of oomph namespace

//==========start_of_main=================================================
/// ??ds
//========================================================================
int main()
{
  // Set number of elements in each direction
  unsigned n_x = 2;
  unsigned n_y = 1;

  // Set dimension of problem and number of nodes along each edge
  const unsigned dim = 2;
  const unsigned nnode_1d = 2;

  // Set up the bulk problem
  TwoDMicromagProblem<QMicromagElement<dim,nnode_1d> >
    problem(n_x,n_y);

  // Set up the boundary problem
  TwoDBoundaryProblem<QMicromagElement<dim,nnode_1d>,MicromagFaceElement>
    boundary_problem(problem.mesh_pt());

  // Setup doc info
  DocInfo doc_info;
  doc_info.set_directory("RESLT");

  // Run self tests on problem and boundary problem:
  std::cout << "\n\n\nProblem self-test ";
  if (problem.self_test()==0)
    std::cout << "passed: Problem can be solved." << std::endl;
  else
    throw OomphLibError("Self test failed","main()",OOMPH_EXCEPTION_LOCATION);
  std::cout << "\n\n\nBoundary problem self-test ";
  if (boundary_problem.self_test()==0)
    std::cout << "passed: Problem can be solved." << std::endl;
  else
    throw OomphLibError("Self test failed","main()",OOMPH_EXCEPTION_LOCATION);

  // Dump boundary mesh positions
  unsigned n_node = boundary_problem.mesh_pt()->nnode();
  std::cout << n_node << std::endl;
  for(unsigned i_node=0; i_node < n_node; i_node++)
    {
      std::cout << boundary_problem.mesh_pt()->node_pt(i_node)->position(0)
  		<< ", " << boundary_problem.mesh_pt()->node_pt(i_node)->position(1)
  		<< std::endl;
    }

  std::cout << std::endl;
  std::cout << "number of degrees of freedom = " << boundary_problem.ndof() << std::endl;

  // unsigned i_ele = 3;
  // {
  //   // Get element pointer to ith element
  //   MicromagFaceElement<QMicromagElement<2,2> >* elem_pt =
  //     dynamic_cast<MicromagFaceElement<QMicromagElement<2,2> >*>
  //     (boundary_problem.mesh_pt()->element_pt(i_ele));

  //   // Create matrix to store boundary data for this element
  //   DenseMatrix<double> boundary_matrix(2,n_node,0.0);
  //   //??ds not actually 2 - should be number of nodes in element

  //   // Get boundary element matrix for this element
  //   Vector<double> dummy(0,0.0);
  //   elem_pt->fill_in_contribution_to_jacobian(dummy,boundary_matrix);

  //   // dump for testing:
  //   std::ofstream matrix_file;
  //   char filename[100] = "matrix";
  //   matrix_file.open(filename);
  //   boundary_matrix.output(matrix_file);
  //   matrix_file.close();
  // }

  DoubleVector dummy_vec;
  DenseDoubleMatrix jacobian;
  boundary_problem.get_boundary_matrix(jacobian);

  // // dump for testing:
  // std::ofstream matrix_file;
  // char filename[100] = "fake_jacobian";
  // matrix_file.open(filename);
  // jacobian.output(matrix_file);
  // matrix_file.close();

  return 0;
} //end of main

#endif
