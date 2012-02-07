#ifndef OOMPH_HYBRID_BOUNDARY_ELEMENT_DRIVER_H
#define OOMPH_HYBRID_BOUNDARY_ELEMENT_DRIVER_H

#include <map>
#include "generic.h"
#include "../micromagnetics_boundary_element.h"
#include "meshes/rectangular_quadmesh.h"

// My variable gauss quadrature header
#include "../variable_quadrature.h"

using namespace oomph;
using namespace MathematicalConstants;

namespace oomph
{
  //=====start_of_TwoDMicromagProblem===============================================
  /// A problem class to solve the LLG equation using a hybrid BEM/FEM.
  //========================================================================
  template<class BULK_ELEMENT, template<class> class FACE_ELEMENT>
  class TwoDMicromagProblem : public Problem
  {

  public:

    /// Constructor
    TwoDMicromagProblem(const unsigned& n_x, const unsigned& n_y);

    /// Destructor (empty -- all the cleanup is done in the base class)
    ~TwoDMicromagProblem(){};

    /// Doc the solution
    void doc_solution(DocInfo& doc_info);

    /// Build the mesh for this problem
    void build_face_mesh(Mesh* face_mesh_pt) const;

    /// Acess to the face mesh pointer
    Mesh* face_mesh_pt() {return Face_mesh_pt;}

    /// Get the boundary element matrix (similar to problem::get_jacobian)
    void get_boundary_matrix();

    /// \short Get the mapping between the global equation numbering and
    /// the boundary equation numbering.
    void create_global_boundary_equation_number_map();

    /// Get the boundary equation number from the global equation number
    unsigned convert_global_to_boundary_equation_number(const unsigned &global_num);

    /// Update the values of phi_2 on the boundary
    void update_boundary_phi_2();

    /// Return pointer to the boundary matrix
    DenseDoubleMatrix* boundary_matrix_pt()
    {return &Boundary_matrix;}

  private:

    /// The map between the global equation numbering and the boundary equation numbering.
    std::map<unsigned,unsigned> Global_boundary_equation_num_map;

    /// Pointer to control node at which the solution is documented
    Node* Control_node_pt;

    /// Mesh containing the face elements
    Mesh* Face_mesh_pt;

    /// Doc info object
    DocInfo Doc_info;

    /// Matrix to store the relationship between phi_1 and phi_2 on the boundary
    DenseDoubleMatrix Boundary_matrix;

    //  /// Trace file
    //  std::ofstream Trace_file;

    // Update the problem before Newton convergence check
    void actions_before_newton_convergence_check()
    {update_boundary_phi_2();}

    /// Update the problem specs before solve
    // Nothing to do here since no dirichlet boundaries.
    void actions_before_newton_solve(){};

    /// Update the problem specs after solve
    void actions_after_newton_solve(){};

    /// Update the problem specs after solve (empty)
    void actions_after_implicit_timestep(){};

    /// Update the problem specs before next timestep
    void actions_before_implicit_timestep(){};

    /// Set initial condition (incl previous timesteps) according to specified function.
    void set_initial_condition(){};

  }; // end of problem class

  //=====start_of_constructor===============================================
  ///
  //========================================================================
  template<class BULK_ELEMENT, template<class> class FACE_ELEMENT>
  TwoDMicromagProblem<BULK_ELEMENT,FACE_ELEMENT>::
  TwoDMicromagProblem(const unsigned& n_x, const unsigned& n_y)
  {
    // Set domain size
    double l_x = 2.0;
    double l_y = 2.0;

    // Allocate steady state timestepper
    add_time_stepper_pt(new Steady<2>);

    // Build mesh
    mesh_pt() = new RectangularQuadMesh<BULK_ELEMENT>
      (n_x,n_y,l_x,l_y,time_stepper_pt());

    // Create empty face mesh
    Face_mesh_pt = new Mesh;

    // Create face elements on all boundaries and add to face mesh.
    build_face_mesh(face_mesh_pt());

    // Loop over elements to set pointers
    // ??ds add this when we are actually solving things

    // Loop over elements to set mesh pointer
    unsigned n_element = face_mesh_pt()->nelement();
    for(unsigned i=0;i<n_element;i++)
      {
	// Upcast from GeneralisedElement to the face element
	FACE_ELEMENT<BULK_ELEMENT>* elem_pt =
	  dynamic_cast<FACE_ELEMENT<BULK_ELEMENT> *>(face_mesh_pt()->element_pt(i));

	// Set boundary mesh pointer in element
	elem_pt->set_mesh_pt(face_mesh_pt());
      }

    // Setup equation numbering scheme
    std::cout << "FEM number of equations: " << assign_eqn_numbers() << std::endl;

  } // end of constructor

  //=====================start_of_doc=======================================
  /// ??ds write this!
  //========================================================================
  template<class BULK_ELEMENT, template<class> class FACE_ELEMENT>
  void TwoDMicromagProblem<BULK_ELEMENT,FACE_ELEMENT>::
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

  //======start_of_convert_global_to_boundary_equation_number===============
  /// ??ds write this!
  //========================================================================
  template<class BULK_ELEMENT, template<class> class FACE_ELEMENT>
  unsigned TwoDMicromagProblem<BULK_ELEMENT,FACE_ELEMENT>::
  convert_global_to_boundary_equation_number(const unsigned &global_num)
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
			    "TwoDMicromagProblem::convert_global_to_boundary_equation_number",
			    OOMPH_EXCEPTION_LOCATION);
      }
#endif
    return ((*Global_boundary_equation_num_map.find(global_num)).second);
  }
  //======start_of_build_face_mesh==========================================
  /// \short Constuct a Mesh of FACE_ELEMENTs along the b-th boundary
  /// of the mesh (which contains elements of type BULK_ELEMENT)
  //========================================================================
  template<class BULK_ELEMENT, template<class> class FACE_ELEMENT>
  void TwoDMicromagProblem<BULK_ELEMENT,FACE_ELEMENT>::
  build_face_mesh(Mesh* face_mesh_pt) const
  {
    // Create a set to temporarily store the list of boundary nodes
    // (we do it via a set because sets automatically detect duplicates)
    std::set<Node*> node_set;
    std::set<Node*>::iterator it;

    // Loop over the boundaries
    unsigned n_boundary = mesh_pt()->nboundary();
    for(unsigned b=0; b<n_boundary; b++)
      {
	//Loop over the boundary nodes on boundary b
	unsigned n_bound_node = mesh_pt()->nboundary_node(b);
	for(unsigned n=0;n<n_bound_node;n++)
	  {
	    // Add the boundary node to a set
	    node_set.insert(mesh_pt()->boundary_node_pt(b,n));
	  }

	//Loop over the elements on boundary b
	unsigned n_bound_element = mesh_pt()->nboundary_element(b);
	for(unsigned e=0;e<n_bound_element;e++)
	  {
	    //Create the corresponding FaceElement
	    FACE_ELEMENT<BULK_ELEMENT>* face_element_pt = new FACE_ELEMENT<BULK_ELEMENT>
	      (mesh_pt()->boundary_element_pt(b,e),
	       mesh_pt()->face_index_at_boundary(b,e));

	    //Add the face element to the face mesh
	    face_mesh_pt->add_element_pt(face_element_pt);
	  }
      }

    // Iterate over all nodes in the set and add to the face mesh
    for(it=node_set.begin(); it!=node_set.end(); it++)
      face_mesh_pt->add_node_pt(*it);

    //??ds taken from mesh.h - no idea what this does but maybe useful later...
#ifdef OOMPH_HAS_MPI
    // If the bulk mesh has been distributed then the face mesh is too
    if (this->is_mesh_distributed())
      {
	face_mesh_pt->set_mesh_distributed();
      }
#endif
  }

  //======start_of_create_global_boundary_equation_number_map===============
  /// Create a map from the global equation numbers to the boundary node.
  /// Note since we use the global equation number as a key the map will
  /// automatically reject duplicate entries.
  //========================================================================
  template<class BULK_ELEMENT, template<class> class FACE_ELEMENT>
  void TwoDMicromagProblem<BULK_ELEMENT,FACE_ELEMENT>::
  create_global_boundary_equation_number_map()
  {
    // Initialise the map
    Global_boundary_equation_num_map.clear();

    // Loop over boundary nodes assigning a boundary equation number to each.
    unsigned n_boundary_node = this->face_mesh_pt()->nnode(), k=0;
    for(unsigned i_node=0; i_node<n_boundary_node; i_node++)
      {
	// Get global equation number for phi_2
	unsigned global_eqn_number = this->face_mesh_pt()->
	  node_pt(i_node)->eqn_number(0);

	// Set up the pair ready to input with key="global equation number" and
	// value ="boundary equation number"=k.
	std::pair<unsigned,unsigned> input_pair = std::make_pair(global_eqn_number,k);

	// Add entry to map and store whether this was a new addition
	bool new_addition = (Global_boundary_equation_num_map.insert(input_pair)
			     ).second;

	// Increment k if this was a new addition to the map
	if(new_addition) k++;
      }
  }

  //=============================================================================
  /// Get the fully assembled boundary matrix in dense storage.
  //=============================================================================
  template<class BULK_ELEMENT, template<class> class FACE_ELEMENT>
  void TwoDMicromagProblem<BULK_ELEMENT,FACE_ELEMENT>::
  get_boundary_matrix()
  {

    // get the number of nodes in the boundary problem
    unsigned long n_node = face_mesh_pt()->nnode();

    // Initialise and resize the boundary matrix
    Boundary_matrix.resize(n_node,n_node);
    Boundary_matrix.initialise(0.0);

    // Loop over all the elements
    unsigned long n_element = face_mesh_pt()->nelement();
    for(unsigned long e=0;e<n_element;e++)
      {
	// Get the pointer to the element (and cast to FiniteElement)
	//??ds might need only finite ele in the end? not sure how to get phi_2_index yet
	FiniteElement* elem_pt =
	  dynamic_cast < FiniteElement* > (face_mesh_pt()->element_pt(e));

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
	    unsigned l_number =
	      this->convert_global_to_boundary_equation_number
	      (elem_pt->node_pt(l)->eqn_number(0));

	    //std::cout << "target node = " << l << std::endl;

	    // Loop over all nodes in the mesh and add contributions from this element
	    for(unsigned long source_node=0; source_node<n_node; source_node++)
	      {
		unsigned source_number =
		  this->convert_global_to_boundary_equation_number
		  (face_mesh_pt()->node_pt(source_node)->eqn_number(0));

		Boundary_matrix(l_number,source_number)
		  += element_boundary_matrix(l,source_node);
	      }
	  }
      }
  }

  //=============================================================================
  ///
  // updating the boundary conditions on phi_2 from the values
  // of phi_1 on the boundary goes in here. This combined with
  // including it in the jacobian allows the solver to work as normal.
  //=============================================================================
  template<class BULK_ELEMENT, template<class> class FACE_ELEMENT>
  void TwoDMicromagProblem<BULK_ELEMENT,FACE_ELEMENT>::
  update_boundary_phi_2()
  {
    // Get the index of phi_2
    BULK_ELEMENT* elem_pt = (dynamic_cast<BULK_ELEMENT*>(mesh_pt()->element_pt(0)));
    unsigned phi_2_index = elem_pt->phi_2_index_micromag();

    // Loop over all (target) nodes on the boundary
    unsigned n_boundary_node = face_mesh_pt()->nnode();
    for(unsigned target_node=0; target_node<n_boundary_node; target_node++)
      {
	// Get a pointer to the target node
	Node* target_node_pt = face_mesh_pt()->node_pt(target_node);

    	// Get boundary equation number for this target node
    	unsigned target_number = convert_global_to_boundary_equation_number
    	  (target_node_pt->eqn_number(0));

	// Double to store the value of phi_2 during computation
	double target_phi_2_value = 0;

    	// Loop over all source nodes adding contribution from each
    	for(unsigned source_node=0; source_node<n_boundary_node; source_node++)
    	  {
    	    // Get a pointer to the source node
    	    Node* source_node_pt = face_mesh_pt()->node_pt(source_node);

    	    // Get boundary equation number for this source node
    	    unsigned source_number = convert_global_to_boundary_equation_number
    	      (source_node_pt->eqn_number(0));

    	    // Add the contribution to phi_2 at the target node due to
    	    // the source node (relationship is given by the boundary matrix).
    	    //??ds check consistency of boundary matrix numbering
	    target_phi_2_value += Boundary_matrix(target_number,source_number)
    	      * source_node_pt->value(0);
    	  }
	// Save the total into the target node
	target_node_pt->set_value(phi_2_index,target_phi_2_value);
      }
  }

} // End of oomph namespace



//==========start_of_main=================================================
/// ??ds
//========================================================================
int main(int argc, char* argv[])
{
  unsigned n_x, n_y;

  // Get inputs
  if(argc >= 2)
    {
      n_x = atoi(argv[1]);
      n_y = atoi(argv[2]);
    }
  else
  {
      // Set number of elements in each direction
      n_x = 20;
      n_y = 20;
  }

  // Set dimension of problem and number of nodes along each element edge
  const unsigned dim = 2;
  const unsigned nnode_1d = 2;

  // Set up the bulk problem
  TwoDMicromagProblem<QMicromagElement<dim,nnode_1d>, MicromagFaceElement>
    problem(n_x,n_y);

  // Setup doc info
  DocInfo doc_info;
  doc_info.set_directory("results");

  // Run self tests on problem and boundary problem:
  problem.self_test();

  // // // Dump boundary mesh positions
  // unsigned n_node = problem.face_mesh_pt()->nnode();
  // std::cout << n_node << std::endl;
  // // for(unsigned i_node=0; i_node < n_node; i_node++)
  // //   {
  // //     std::cout << boundary_problem.mesh_pt()->node_pt(i_node)->position(0)
  // // 		<< ", " << boundary_problem.mesh_pt()->node_pt(i_node)->position(1)
  // // 		<< std::endl;
  // //   }

  // Set up the boundary equation numbering system
  problem.create_global_boundary_equation_number_map();




  // // ??ds testing my variable gaussian scheme
  VariableClenshawCurtis variable_gauss;
  variable_gauss.set_dim(1);

  // Set all boundary elements to use the variable order gauss integration
  unsigned n_element = problem.face_mesh_pt()->nelement();
  for(unsigned i=0; i<n_element; i++)
    {
      FiniteElement* finite_element_pt =
  	dynamic_cast<FiniteElement*>(problem.face_mesh_pt()->element_pt(i));
      finite_element_pt->set_integration_scheme(&variable_gauss);
    }

  unsigned max_order = 30;

  for(unsigned order=2; order<max_order; order++)
    {
      // Set the integration scheme order
      variable_gauss.set_order(order);

      // Get the boundary matrix
      problem.get_boundary_matrix();

      // dump for testing
      std::ofstream matrix_file;
      matrix_file.setf(std::ios::fixed,std::ios::floatfield); // Set high precision output
      matrix_file.precision(16);
      char filename[100];
      sprintf(filename,"results/boundary_matrix_%u",order);
      matrix_file.open(filename);
      problem.boundary_matrix_pt()->output(matrix_file);
      matrix_file.close();
    }


  return 0;
} //end of main

#endif
