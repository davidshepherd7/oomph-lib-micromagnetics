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
    Problem::mesh_pt() = new Mesh();

    // control node?

    // Create face elements on all boundaries and add them to the face mesh
    unsigned n_bulk_boundary = bulk_mesh_pt()->nboundary();
    for(unsigned i_bulk_boundary=0; i_bulk_boundary < n_bulk_boundary;
	i_bulk_boundary++)
      {
	Bulk_mesh_pt->build_face_mesh<BULK_ELEMENT,FACE_ELEMENT>
	  (i_bulk_boundary,Problem::mesh_pt());
      }

    //??ds remove the duplicate nodes from the face mesh
    //??ds remove + re-add via a set?

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
  // std::cout << "\n\n\nBoundary problem self-test ";
  // if (boundary_problem.self_test()==0)
  //   std::cout << "passed: Problem can be solved." << std::endl;
  // else
  //   throw OomphLibError("Self test failed","main()",OOMPH_EXCEPTION_LOCATION);

  // Dump boundary mesh positions
  unsigned n_node = boundary_problem.mesh_pt()->nnode();
  std::cout << n_node << std::endl;
  for(unsigned i_node=0; i_node < n_node; i_node++)
    {
      std::cout << boundary_problem.mesh_pt()->node_pt(i_node)->position(0)
  		<< ", " << boundary_problem.mesh_pt()->node_pt(i_node)->position(1)
  		<< std::endl;
    }

  // // solve the problem
  // problem.newton_solve();

  // //Output solution
  // problem.doc_solution(doc_info);

  // // Get "jacobian" (actually boundary element matrix) for boundary problem
  // DenseDoubleMatrix* matrix_pt=new DenseDoubleMatrix;
  // DoubleVector dummy;
  // boundary_problem.get_jacobian(dummy,*matrix_pt);

  // Eventually loop over all boundary elements?
  unsigned n_boundary_node = boundary_problem.mesh_pt()->nnode();

  unsigned i_ele = 0;
  {
    // Get element pointer to ith element
    MicromagFaceElement<QMicromagElement<2,2> >* elem_pt =
      dynamic_cast<MicromagFaceElement<QMicromagElement<2,2> >*>
      (boundary_problem.mesh_pt()->element_pt(i_ele));

    // Create matrix to store boundary data for this element
    DenseMatrix<double> boundary_matrix(2,n_boundary_node,0.0);
    //??ds not actually 2 - should be number of nodes in element

    // Get boundary element matrix for this element
    Vector<double> dummy(0,0.0);
    elem_pt->fill_in_contribution_to_jacobian(dummy,boundary_matrix);

    // dump for testing:
    std::ofstream matrix_file;
    char filename[100] = "matrix";
    matrix_file.open(filename);
    boundary_matrix.output(matrix_file);
    matrix_file.close();
  }


  return 0;
} //end of main

#endif
