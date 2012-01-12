#ifndef OOMPH_HYBRID_BOUNDARY_ELEMENT_DRIVER_H
#define OOMPH_HYBRID_BOUNDARY_ELEMENT_DRIVER_H

/*
  ??ds  description of file goes here
*/

#include <map>
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
    std::cout << "FEM number of equations: " << assign_eqn_numbers() << std::endl;

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

    /// \short Get the mapping between the global equation numbering and
    /// the boundary equation numbering.
    void create_global_boundary_equation_number_map();

    /// Get the boundary equation number from the global equation number
    unsigned convert_global_to_boundary_equation_number(const unsigned &global_num)
    {
#ifdef PARANOID
      // Get the location of the global_num key in an iterator
      std::map<unsigned,unsigned>::iterator it
	= Global_boundary_equation_num_map.find(global_num);

      // If the iterator is placed at the end the given global equation number is
      // not in the map, so return an error.
      if(it == Global_boundary_equation_num_map.end())
	{
	  std::ostringstream error_stream;
	  error_stream << "Global equation number " << global_num
			<< " is not in the global to boundary map.";
	  throw OomphLibError(error_stream.str(),
			      "TwoDBoundaryProblem::convert_global_to_boundary_equation_number",
			      OOMPH_EXCEPTION_LOCATION);
	}
#endif
      return ((*Global_boundary_equation_num_map.find(global_num)).second);
    }

  private:

    /// Create face elements on the i'th boundary and add to boundary mesh
    void create_face_elements(const unsigned &i_bulk_boundary);

    /// The map the global equation numbering and the boundary equation numbering.
    std::map<unsigned,unsigned> Global_boundary_equation_num_map;

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
    std::cout << "Boundary element number of equations: " << assign_eqn_numbers() << std::endl;

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

  // Create a map from the global equation numbers to the boundary node.
  // Note since we use the global equation number as a key the map will
  // automatically reject duplicate entries.
  template<class BULK_ELEMENT, template<class> class FACE_ELEMENT>
  void TwoDBoundaryProblem<BULK_ELEMENT,FACE_ELEMENT>::
  create_global_boundary_equation_number_map()
  {
    // Initialise the map
    Global_boundary_equation_num_map.clear();

    // Get index of phi_2 (from the first element of the mesh)
    BULK_ELEMENT* ele_pt = dynamic_cast < BULK_ELEMENT* > (this->mesh_pt()->element_pt(0));
    unsigned phi_2_index = ele_pt->phi_2_index_micromag();

    // Loop over boundary nodes assigning a boundary equation number to each.
    unsigned n_boundary_node = this->mesh_pt()->nnode(), k=0;
    for(unsigned i_node=0; i_node<n_boundary_node; i_node++)
      {
	// Get global equation number for phi_2
	unsigned global_eqn_number = this->mesh_pt()->
	  node_pt(i_node)->eqn_number(phi_2_index);

	// Set up the pair ready to input with key="global equation number" and
	// value ="boundary equation number"=k.
	std::pair<unsigned,unsigned> input_pair = std::make_pair(global_eqn_number,k);

	// Add entry to map and store whether this was a new addition
	bool new_addition = (Global_boundary_equation_num_map.insert(input_pair)
			).second;

	// Increment k if this was a new addition to the map
	if(new_addition) k++;
      }


    // std::cout << "Just messing: " << std::endl;

    // for(unsigned i=0; i<50; i++)
    //   {
    // 	std::cout << convert_global_to_boundary_equation_number(i) << " "
    // 		  << std::endl;
    //   }
  }


  //=============================================================================
  /// Get the fully assembled boundary matrix in dense storage.
  //=============================================================================
  template<class BULK_ELEMENT, template<class> class FACE_ELEMENT>
  void TwoDBoundaryProblem<BULK_ELEMENT,FACE_ELEMENT>::
  get_boundary_matrix(DenseDoubleMatrix& boundary_matrix)
  {

    // get the number of nodes in the boundary problem
    unsigned long n_node = mesh_pt()->nnode();

    // resize the boundary matrix
    boundary_matrix.resize(n_node,n_node);
    boundary_matrix.initialise(0.0);

    // Loop over all the elements
    unsigned long n_element = mesh_pt()->nelement();
    for(unsigned long e=0;e<n_element;e++)
      {
	// Get the pointer to the element (and cast to FiniteElement)
	//??ds might need only finite ele in the end? not sure how to get phi_2_index yet
	FACE_ELEMENT<BULK_ELEMENT>* elem_pt =
	  dynamic_cast<FACE_ELEMENT<BULK_ELEMENT>* >(mesh_pt()->element_pt(e));

	// Find number of nodes in the element
	unsigned long n_element_node = elem_pt->nnode();

	// Set up a matrix and dummy residual vector
	DenseMatrix<double> element_boundary_matrix(n_element_node,n_node);
	Vector<double> dummy(0);

	// Fill the matrix
	assembly_handler_pt()
	  ->get_jacobian(elem_pt,dummy,element_boundary_matrix);

	// Loop over the nodes in this element
	for(unsigned l=0;l<n_element_node;l++)
	  {

	    // Get the boundary equation (=node) number from the global one
	    unsigned source_node =
	      this->convert_global_to_boundary_equation_number(
							       elem_pt->node_pt(l)->eqn_number(8));
	    //??ds need to get real phi_2_micromag somehow... :(

	    //std::cout << "source node = " << source_node << std::endl;

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
int main(int argc, char* argv[])
{
  // Set number of elements in each direction
  unsigned n_x = 3;
  unsigned n_y = 3;

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
  doc_info.set_directory("results");

  // Run self tests on problem and boundary problem:
  problem.self_test();
  boundary_problem.self_test();

  // // Dump boundary mesh positions
  unsigned n_node = boundary_problem.mesh_pt()->nnode();
  std::cout << n_node << std::endl;
  // for(unsigned i_node=0; i_node < n_node; i_node++)
  //   {
  //     std::cout << boundary_problem.mesh_pt()->node_pt(i_node)->position(0)
  // 		<< ", " << boundary_problem.mesh_pt()->node_pt(i_node)->position(1)
  // 		<< std::endl;
  //   }

  // Setup the boundary equation numbering map
  boundary_problem.create_global_boundary_equation_number_map();

  DoubleVector dummy_vec;
  DenseDoubleMatrix jacobian;
  boundary_problem.get_boundary_matrix(jacobian);

  // dump for testing:
  std::ofstream matrix_file;
  char filename[100] = "fake_jacobian";
  matrix_file.open(filename);
  jacobian.output(matrix_file);
  matrix_file.close();

  return 0;
} //end of main

#endif
