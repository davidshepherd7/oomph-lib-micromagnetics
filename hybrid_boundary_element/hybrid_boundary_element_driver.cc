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
  class TwoDBoundaryProblem : public Problem
  {

  public:

    /// Constructor: Pass number of elements and pointer to source function
    TwoDBoundaryProblem(const unsigned& n_x, const unsigned& n_y);

    /// Destructor (empty -- all the cleanup is done in the base class)
    ~TwoDBoundaryProblem(){};

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

    /// Create face elements on the b-th boundary
    void create_face_elements(const unsigned& b);

    /// Number of bulk elements
    unsigned N_bulk_element;

  }; // end of problem class

  // /// Set number of values stored at each node
  // template<unsigned DIM, unsigned NNODE_1D>
  //   const unsigned MicromagFaceElement<NNODE_1D>::Initial_Nvalue = 2;

  //=====start_of_constructor===============================================
  ///
  //========================================================================
  template<class ELEMENT>
  TwoDBoundaryProblem<ELEMENT>::
  TwoDBoundaryProblem(const unsigned& n_x, const unsigned& n_y)
  {
    // Set domain size
    double l_x = 4.0;
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

    // Create face elements on all elements adjacent to all boundaries
    unsigned n_boundary = mesh_pt()->nboundary();
    for(unsigned i_boundary=0; i_boundary<n_boundary; i_boundary++)
      create_face_elements(i_boundary);

    //??ds boundary conditions?
    // No Dirichlet conditions so nothing to do here

    // Loop over elements to set pointers
    // No pointers to set so nothing to do here

    // Setup equation numbering scheme
    std::cout << "Number of equations: " << assign_eqn_numbers() << std::endl;

  } // end of constructor

  //============start_of_create_face_elements==============================
  /// Create Face Elements on the b-th boundary of the Mesh.
  //=======================================================================
  template<class ELEMENT>
  void TwoDBoundaryProblem<ELEMENT>::
  create_face_elements(const unsigned &b)
  {
    // How many bulk elements are on the boundary?
    unsigned n_boundary_element = mesh_pt()->nboundary_element(b);

    // Loop over these elements
    for(unsigned e=0; e<n_boundary_element; e++)
      {
	// Get pointer to bulk element
	ELEMENT* bulk_element_pt = dynamic_cast<ELEMENT*>
	  (mesh_pt()->boundary_element_pt(b,e));

	// What is the index of the face of the bulk element at the boundary
	int face_index = mesh_pt()->face_index_at_boundary(b,e);

	// Build the corresponding face element
	MicromagFaceElement<ELEMENT>* face_element_pt =
	  new MicromagFaceElement<ELEMENT>(bulk_element_pt,face_index,mesh_pt());

	//Add the prescribed-flux element to the mesh
	mesh_pt()->add_element_pt(face_element_pt);
      }

  }

  //=====================start_of_doc=======================================
  /// ??ds write this!
  //========================================================================
  template<class ELEMENT>
  void TwoDBoundaryProblem<ELEMENT>::
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

} // End of oomph namespace

//==========start_of_main=================================================
/// ??ds
//========================================================================
int main()
{
  // Set number of elements in each direction
  unsigned n_x = 2;
  unsigned n_y = 1;


  //Set up the problem
  TwoDBoundaryProblem<QMicromagElement<2,2> >
    problem(n_x,n_y);

  // Setup doc info
  DocInfo doc_info;
  doc_info.set_directory("RESLT");

  // Check that we're ready to go
  std::cout << "\n\n\nProblem self-test ";
  if (problem.self_test()==0)
    {
      std::cout << "passed: Problem can be solved." << std::endl;
    }
  else
    {
      throw OomphLibError("Self test failed",
			  "main()",
			  OOMPH_EXCEPTION_LOCATION);
    }

  // // Solve the problem
  // problem.newton_solve();

  // //Output solution
  // problem.doc_solution(doc_info);

  // get boundary matrix contribution from first element

  // Eventually loop over all boundary elements?
  unsigned i_ele = problem.get_n_bulk_element();
  {
    // Get element pointer to ith element
    MicromagFaceElement<QMicromagElement<2,2> >* elem_pt =
      dynamic_cast<MicromagFaceElement<QMicromagElement<2,2> >*>
      (problem.mesh_pt()->element_pt(i_ele));

    // Get total number of boundary nodes
    unsigned n_boundary = problem.mesh_pt()->nboundary();
    unsigned n_boundary_nodes(0.0);
    for(unsigned i_boundary=0; i_boundary<n_boundary; i_boundary++)
      n_boundary_nodes += problem.mesh_pt()->nboundary_node(i_boundary);
    std::cout << "There are " << n_boundary_nodes << " nodes on the boundary" << std::endl;

    // Create matrix to store boundary data for this element
    DenseMatrix<double> boundary_matrix(2,n_boundary_nodes,0.0); //??ds not actually 2, jsut for now

    // Get boundary element matrix for this element
    elem_pt->fill_in_elemental_contribution_to_boundary_element_matrix(boundary_matrix);

    // dump for testing:
    std::ofstream matrix_file;
    char filename[100] = "matrix";
    matrix_file.open(filename);
    boundary_matrix.output(matrix_file);
    matrix_file.close();
  }

} //end of main

#endif
