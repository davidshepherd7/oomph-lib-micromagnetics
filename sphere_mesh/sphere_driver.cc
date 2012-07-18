#ifndef OOMPH_COMBINED_HYBRID_DRIVER_H
#define OOMPH_COMBINED_HYBRID_DRIVER_H

/*
  description of file goes here
*/

// Include all my stuff
#include "../my_general_header.h"

// Mesh
#include "meshes/tetgen_mesh.h"

#include "meshes/rectangular_quadmesh.h"

using namespace oomph;
using namespace MathematicalConstants;


namespace Inputs
{
  // ??ds need to do this input stuff properly..

  double llg_damping = 1.5;
  double llg_precession = 1.0;
  double exchange_coeff = 0.5;
  double magnetostatic_coeff = 1;
  double crystal_anis_coeff = 0;
  double dt = 1e-5;
  double tmax = 50;

  bool adaptive_timestepping = 1;
  bool full_jacobian_fd = 0;
  bool GMRES = 1;
  bool midpointmethod = 0;
  const unsigned bdforder = 2;

  const unsigned dim = 2;
  // If changing dim you also probably need to swap element type
  // (QMicromagElement vs TMicromagElement) because the mesh changes type.



  void applied_field(const double& t, const Vector<double>& x, Vector<double>& h_app)
  {
    h_app.assign(3,0.0);
    h_app[0] = 0;
    h_app[1] = 0;
    h_app[2] = 0;
  }

  void cryst_anis_field(const double& t, const Vector<double>& x,
			const Vector<double>& m, Vector<double>& h_ca)
  {
    Vector<double> easy_axis(3,0.0); easy_axis[0] = 1.0;

    h_ca = easy_axis;
    double magnitude = dot(easy_axis,m);
    for(unsigned i=0; i<h_ca.size(); i++)
      h_ca[i] *= magnitude * crystal_anis_coeff;
  }

  void initial_m(const double& t, const Vector<double>& x,
		 Vector<double>& m)
  {
    m.assign(3,0.0);
    m[0] = -1;
    m[1] = -0.3;
    m[2] = -0.1;
    normalise(m);
  }
}

namespace oomph
{

  //======================================================================
  /// A problem class to test the combination of the hybrid BEM/FEM for
  /// magnetostatic fields and the LLG equations for micromagnetics.
  //======================================================================
  template<class BULK_ELEMENT,
	   template<class BULK_ELEMENT,unsigned DIM> class BEM_ELEMENT,
	   unsigned DIM>
  class ThreeDHybridProblem : public Problem
  {

  public:

    /// Constructor
    ThreeDHybridProblem(const std::string& node_file_name,
			const std::string& element_file_name,
			const std::string& face_file_name);

    /// Destructor (empty -- once the problem is done with the program is over)
    ~ThreeDHybridProblem(){};

    /// Doc the solution
    void doc_solution(DocInfo& doc_info);

    /// \short Create the face elements to apply flux boundary conditions to the
    /// potential on boundary b.
    void create_flux_elements(const unsigned& b, Mesh* const &bulk_mesh_pt,
			      Mesh* const &surface_mesh_pt);


    /// Build the meshes of bem elements
    void build_bem_mesh(Mesh* bem_mesh_pt) const;

    /// Get the boundary element matrix (similar to problem::get_jacobian)
    void build_boundary_matrix();

    /// Set initial condition (incl previous timesteps)
    void set_initial_condition();


    /// \short Get the mapping between the global equation numbering and
    /// the boundary equation numbering.
    void create_global_boundary_equation_number_maps();

    /// Get the boundary equation number from the global equation number
    unsigned convert_global_to_boundary_equation_number(const unsigned &global_num);

    /// Update the values of phi on the boundary
    void update_boundary_phi();


    /// Access to the pointer to the boundary element method mesh
    Mesh* bem_mesh_pt() const {return Bem_mesh_pt;}

    /// Access to the boundary matrix
    DenseDoubleMatrix* boundary_matrix_pt()
    {return &Boundary_matrix;}

    /// Dummy - always reduce timestep as much as possible if fails to converge.
    double global_temporal_error_norm()
    {
      return 1e-8;
    }

    // //?? check if node is corner, nasty!
    // bool corner(Vector<double> x)
    // {
    //   bool temp(1);
    //   for(unsigned i=0; i<DIM; i++)
    // 	{
    // 	  temp = temp && ( (x[i] == 0) || (x[i] == 1) );
    // 	}
    //   return temp;
    // }

    /// Overload get_jacobian to include the boundary matrix and re-arrange the
    /// sparse part of the Jacobian into block form (as well as getting the
    /// Jacobian as normal).  Note that this is ONLY called if using GMRES -
    /// otherwise the newton solver will just call the standard get_jacobian
    /// since it does not know about the SumOfMatrices format.
    void get_jacobian(DoubleVector& residuals, SumOfMatrices& jacobian)
    {
      std::cout << "Calling your get_jacobian function." << std::endl;

      // Get the sparse part of the Jacobian due to finite elements and set it
      // as the main (first) matrix of the sum.  ??ds this is a memory leak,
      // need a flag to say if sum can delete individual matrices?
      CRDoubleMatrix* sparse_jacobian_pt = new CRDoubleMatrix;
      Problem::get_jacobian(residuals,*sparse_jacobian_pt);
      jacobian.main_matrix_pt() = sparse_jacobian_pt;

      // Set the sparse part of the Jacobian to be deleted along with the sum of
      // the matrices - avoid a memory leak.
      jacobian.set_delete_main_matrix();

      // Add the boundary element matrix to the total Jacobian. It represents
      // the derivative of phi with respect to phi_1 so each entry goes in the
      // phi row and the phi_1 column of the respective element (this is done
      // via the two maps).
      jacobian.add_matrix(boundary_matrix_pt(),
			  &Global_phi_1_num_map,
			  &Global_boundary_equation_num_map,0);



      // //??temp dump jacobian for tests
      // std::ofstream matrix_file;
      // matrix_file.precision(16);
      // char filename[100];
      // sprintf(filename,"results/jacobian");
      // matrix_file.open(filename);
      // sparse_jacobian_pt->sparse_indexed_output(matrix_file);
      // matrix_file.close();

    }

    /// Get the index of phi for use in BEM mapping
    unsigned phi_index() {return Phi_index;}

    /// Get the index of phi_1 for use in BEM mapping
    unsigned phi_1_index() {return Phi_1_index;}


  private:

    /// The map between the global equation numbers for phi and the boundary
    /// equation/matrix numbering.
    std::map<long unsigned,long unsigned> Global_boundary_equation_num_map;

    /// The map between the global equation numbers for phi_1 and the boundary
    /// equation/matrix numbering.
    std::map<long unsigned,long unsigned> Global_phi_1_num_map;

    Mesh* Bulk_mesh_pt;

    Mesh* Flux_mesh_pt;

    /// The pointer to the boundary element method mesh
    Mesh* Bem_mesh_pt;

    /// Doc info object
    DocInfo Doc_info;

    /// Store the index of phi for use in BEM mapping
    unsigned Phi_index;

    /// Store the index of phi_1 for use in BEM mapping
    unsigned Phi_1_index;

    /// Matrix to store the relationship between phi_1 and phi on the boundary
    DenseDoubleMatrix Boundary_matrix;

    /// Update the problem before Newton convergence check (update boundary
    /// conditions on phi).
    void actions_before_newton_convergence_check()
    {
      if(Inputs::magnetostatic_coeff != 0.0) update_boundary_phi();
    }

    /// Update the problem specs before solve.
    void actions_before_newton_solve(){}

    /// Update the problem specs after solve
    void actions_after_newton_solve();

    /// Update the problem specs after solve (empty)
    void actions_after_implicit_timestep(){}
  }; // end of problem class

  //======================================================================
  /// Constructor
  //======================================================================
  template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
  ThreeDHybridProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
  ThreeDHybridProblem(const std::string& node_file_name,
		      const std::string& element_file_name,
		      const std::string& face_file_name)
  {

    if(Inputs::midpointmethod)
      {
	if(Inputs::adaptive_timestepping)
	  std::cerr << "No adaptive midpoint yet." << std::endl;

	add_time_stepper_pt(new BDF<1>);
      }
    else // use bdf
      {
	if(Inputs::adaptive_timestepping)
	  add_time_stepper_pt(new BDF<Inputs::bdforder>(true));
	else
	  add_time_stepper_pt(new BDF<Inputs::bdforder>);
      }

    if(Inputs::full_jacobian_fd)
      {
	linear_solver_pt() = new FD_LU;
      }
    else if (Inputs::GMRES)
      {
	linear_solver_pt() = new GMRES<SumOfMatrices>;
      }
    else
      {
	linear_solver_pt() = new SuperLUSolver;
      }


    if(Inputs::dim == 3)
      {
	// Build mesh from tetgen
	Bulk_mesh_pt = new TetgenMesh<BULK_ELEMENT>(node_file_name, element_file_name,
						    face_file_name, time_stepper_pt());
      }
    else if (Inputs::dim == 2)
      {
	//??temp try a 2d mesh
	unsigned nx = 5, ny = 5;
	Bulk_mesh_pt = new RectangularQuadMesh<BULK_ELEMENT>
	  (nx,ny,1.0,1.0,time_stepper_pt());
      }
    else
      {
	std::cerr << "dim must be 2 or 3" << std::endl;
	throw 123413;
      }

    // Bulk elements
    //------------------------------------------------------------

    // Loop over elements in bulk mesh to set function pointers
    for(unsigned i=0; i< Bulk_mesh_pt->nelement(); i++)
      {
	// Upcast from GeneralisedElement to the present element
	BULK_ELEMENT* elem_pt = dynamic_cast<BULK_ELEMENT*>(Bulk_mesh_pt->element_pt(i));

	// Set pointer to continuous time
	elem_pt->time_pt() = time_pt();

	// Set the function pointers for parameters
	//??ds fix this to use proper encapsulation asap
	elem_pt->applied_field_pt() = &Inputs::applied_field;
	elem_pt->cryst_anis_field_pt() = &Inputs::cryst_anis_field;
	// 	elem_pt->sat_mag_pt() = &Inputs::sat_mag;
	elem_pt->llg_damp_pt() = &Inputs::llg_damping;
	elem_pt->llg_precess_pt() = &Inputs::llg_precession;
	elem_pt->exchange_coeff_pt() = &Inputs::exchange_coeff;
	elem_pt->magnetostatic_coeff_pt() = &Inputs::magnetostatic_coeff;
      }

    // Get the indicies for phi and phi_1 via casting a pointer
    BULK_ELEMENT* elem_pt = dynamic_cast<BULK_ELEMENT*>(Bulk_mesh_pt->element_pt(0));
    Phi_index = elem_pt->phi_index_micromag();
    Phi_1_index = elem_pt->phi_1_index_micromag();

    // Pin the values of phi on the boundary nodes (since it is a Dirichlet //
    // boundary set from the value of phi_1 and the boundary matrix).
    // BULK_ELEMENT* some_el_pt = dynamic_cast< BULK_ELEMENT* >
    // (Bulk_mesh_pt->element_pt(0)); for(unsigned b=0;
    // b<Bulk_mesh_pt->nboundary(); b++) { for(unsigned nd=0; nd <
    // Bulk_mesh_pt->nboundary_node(b); nd++) {
    // Bulk_mesh_pt->boundary_node_pt(b,nd)->
    // pin(some_el_pt->phi_index_micromag()); } }

    // //??temp to help with testing phi_1 pin it's value to zero at r cos(azi) sin(polar) = 0,
    // // i.e when r =0.
    // bool found = false;
    // for(unsigned nd=0; nd< Bulk_mesh_pt->nnode(); nd++)
    //   {
    // 	Vector<double> nd_x(DIM,0.0);
    // 	for(unsigned j=0; j<DIM; j++)
    // 	  nd_x[j] = Bulk_mesh_pt->node_pt(nd)->x(j);
    // 	if ( small(nd_x[0]) && small(nd_x[1]) && small(1 - nd_x[2]) )
    // 	  {
    // 	    Bulk_mesh_pt->node_pt(nd)->pin(some_el_pt->phi_1_index_micromag());
    // 	    Bulk_mesh_pt->node_pt(nd)->set_value(some_el_pt->phi_1_index_micromag(),0.0);
    // 	    found = true;
    // 	  }
    //   }
    // if (!(found))
    //   throw OomphLibError("No node enar middle","",OOMPH_EXCEPTION_LOCATION);



    // Flux elements
    //------------------------------------------------------------

    // We want Neumann (flux) boundary condtions on phi_1 on all boundaries so
    // create the face elements needed.
    Flux_mesh_pt = new Mesh;
    for(unsigned b=0; b < Bulk_mesh_pt->nboundary(); b++)
      {
	create_flux_elements(b,Bulk_mesh_pt,Flux_mesh_pt);
      }

    // Build global finite element mesh (not including the bem mesh)
    add_sub_mesh(Bulk_mesh_pt);
    add_sub_mesh(Flux_mesh_pt);
    build_global_mesh();

    // Setup equation numbering scheme for all the finite elements
    std::cout << "FEM number of equations: " << assign_eqn_numbers() << std::endl;

    if(Inputs::magnetostatic_coeff != 0.0)
      {
	// BEM elements
	//------------------------------------------------------------

	// Create integration scheme in case it is needed.
	QVariableOrderGaussLegendre<DIM-1>* variable_scheme_pt
	  = new QVariableOrderGaussLegendre<DIM-1>;

	// Create BEM elements on all boundaries and add to BEM mesh
	Bem_mesh_pt = new Mesh;
	build_bem_mesh(bem_mesh_pt());

	// Set boundary mesh pointer in element (static memeber so only do once)
	BEM_ELEMENT<BULK_ELEMENT,DIM>::set_boundary_mesh_pt(bem_mesh_pt());

	// Set pointers in all elements
	for(unsigned i_ele = 0; i_ele < bem_mesh_pt()->nelement(); i_ele++)
	  {
	    // Cast pointer
	    BEM_ELEMENT<BULK_ELEMENT,DIM>* ele_pt
	      = dynamic_cast<BEM_ELEMENT<BULK_ELEMENT,DIM>* >
	      (bem_mesh_pt()->element_pt(i_ele));

	    // Set integration scheme
	    ele_pt->set_integration_scheme(variable_scheme_pt);
	  }

	// Make the boundary matrix (including setting up the numbering scheme)
	build_boundary_matrix();
      }

  } // end of constructor


  //======================================================================
  /// Create potential flux boundary condition elements on boundary b.
  //======================================================================
  template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
  void ThreeDHybridProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
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

	// Add the prescribed-flux element to the mesh
	surface_mesh_pt->add_element_pt(flux_element_pt);

      } // End of loop over bulk elements adjacent to boundary b
  }

  //======================================================================
  /// Output function
  //======================================================================
  template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
  void ThreeDHybridProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
  doc_solution(DocInfo& doc_info)
  {
    // Number of plot points
    unsigned npts=2;

    // File set up
    std::ofstream some_file;
    char filename[100];
    sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
	    doc_info.number());

    // Output
    some_file.open(filename);
    mesh_pt()->output(some_file,npts);
    some_file.close();
  } // end of doc


  //======================================================================
  /// Given the global equation number return the boundary equation number. Most
  /// of the function is an error check.
  //======================================================================
  template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
  unsigned ThreeDHybridProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
  convert_global_to_boundary_equation_number(const unsigned &global_num)
  {
#ifdef PARANOID
    // Get the location of the global_num key in an iterator
    std::map<long unsigned, long unsigned>::iterator it
      = Global_boundary_equation_num_map.find(global_num);

    // If the iterator is placed at the end the given global equation number is
    // not in the map, so return an error.
    if(it == Global_boundary_equation_num_map.end())
      {
	std::ostringstream error_stream;
	error_stream << "Global equation number " << global_num
		     << " is not in the global to boundary map.";
	throw OomphLibError(error_stream.str(),
			    "ThreeDHybridProblem::convert_global_to_boundary_equation_number",
			    OOMPH_EXCEPTION_LOCATION);
      }

    //??ds I should use something other than just "equation number 1" for my
    //boundary matrix numbering scheme. doesn't matter what so long as it's
    //consistent and unique to each node. Problems could occur if equation number 1
    //ever has pinned values...
    if(global_num < 0)
      throw OomphLibError("Pinned equation in equation no one, use a different eq num?",
			  "ThreeDHybridProblem::convert_global_to_boundary_equation_number",
			  OOMPH_EXCEPTION_LOCATION);
#endif
    return ((*Global_boundary_equation_num_map.find(global_num)).second);
  }


  //======================================================================
  /// Build the mesh of bem elements.
  //======================================================================
  template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
  void ThreeDHybridProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
  build_bem_mesh(Mesh* bem_mesh_pt) const
  {
    // Create a set to temporarily store the list of boundary nodes
    // (we do it via a set because sets automatically detect duplicates)
    std::set<Node*> node_set;
    std::set<Node*>::iterator it, c_it;

    // Loop over the boundaries
    unsigned n_boundary = Bulk_mesh_pt->nboundary();
    for(unsigned b=0; b<n_boundary; b++)
      {
	// Loop over the boundary nodes on boundary b making a set of nodes
	unsigned n_bound_node = Bulk_mesh_pt->nboundary_node(b);
	for(unsigned n=0;n<n_bound_node;n++)
	  node_set.insert(Bulk_mesh_pt->boundary_node_pt(b,n));

	// Loop over the elements on boundary b creating bem elements
	unsigned n_bound_element = Bulk_mesh_pt->nboundary_element(b);
	for(unsigned e=0;e<n_bound_element;e++)
	  {
	    // Create the corresponding BEM Element
	    BEM_ELEMENT<BULK_ELEMENT,DIM>* bem_element_pt = new BEM_ELEMENT<BULK_ELEMENT,DIM>
	      (Bulk_mesh_pt->boundary_element_pt(b,e),
	       Bulk_mesh_pt->face_index_at_boundary(b,e));

	    // Add the BEM element to the BEM mesh
	    bem_mesh_pt->add_element_pt(bem_element_pt);
	  }
      }

    // Iterate over all nodes in the set and add them to the BEM mesh
    for(it=node_set.begin(); it!=node_set.end(); it++)
      bem_mesh_pt->add_node_pt(*it);
  }


  //======================================================================
  /// Create the map between the global equation numbering system and the
  /// boundary equation numbering system.
  //======================================================================
  template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
  void ThreeDHybridProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
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

    // std::cout << Global_boundary_equation_num_map << std::endl;
  }

  //=============================================================================
  /// Get the fully assembled boundary matrix in dense storage.
  //=============================================================================
  template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
  void ThreeDHybridProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
  build_boundary_matrix()
  {
    // Create the mapping from global to boundary equations
    create_global_boundary_equation_number_maps();

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
	      this->convert_global_to_boundary_equation_number
	      (elem_pt->node_pt(l)->eqn_number(phi_index()));

	    // Loop over all nodes in the mesh and add contributions from this element
	    for(unsigned long source_node=0; source_node<n_node; source_node++)
	      {
		unsigned source_number =
		  this->convert_global_to_boundary_equation_number
		  (bem_mesh_pt()->node_pt(source_node)->eqn_number(phi_index()));

		Boundary_matrix(l_number,source_number)
		  -= element_boundary_matrix(l,source_node);
		//??Ds I think the sign here is negative because the lindholm
		//formula is for +ve but our final equation has negative
		//kernel...
	      }
	  }
      }

    // Lindholm formula does not contain the solid angle contribution: add in
    // here loop over the matrix diagonals adding the angle factor.
    for(unsigned long nd = 0; nd < bem_mesh_pt()->nnode(); nd++)
      {
	// For a sphere all points are smooth (at least in the limit of infinite
	// refinement) so the solid angle contribution is (2*pi)/(4*pi) = 0.5.
	Boundary_matrix(nd,nd) += 0.5; //??temp messing with this..
      }

  }

  //======================================================================
  /// get the new boundary condition on phi using the boundary element matrix
  /// and phi_1.
  //======================================================================
  template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
  void ThreeDHybridProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
  update_boundary_phi()
  {
    //??temp, maybe If magnetostatic coeff is zero then we are ignoring the
    // magnetostatic field so leave boundry pinned at zero?
    // Not a good way to ignore ms... should somehow remove eqns alltogether...
    if (Inputs::magnetostatic_coeff == 0)
      return;

    // In order to use DoubleVector (for matrix multiplications) we need to have
    // this thingy. In our serial case it just gives a number of rows.
    LinearAlgebraDistribution dummy(0,bem_mesh_pt()->nnode(),false);

    // Assemble a vector of phi_1 values on boundary nodes
    DoubleVector boundary_phi_1(dummy,0.0);
    DoubleVector exact_phi_1(dummy,0.0);
    for(unsigned i_nd=0; i_nd < bem_mesh_pt()->nnode(); i_nd++)
      {
	// Get boundary equation number
	unsigned nd_equation = convert_global_to_boundary_equation_number
	  ( bem_mesh_pt()->node_pt(i_nd)->eqn_number(phi_index()) );

	// Copy into appropriate vector entry
	boundary_phi_1[nd_equation] = bem_mesh_pt()->node_pt(i_nd)->value(phi_1_index());

	// //??temp
	// Vector<double> nd_x(DIM,0.0);
	// for(unsigned j=0; j<DIM; j++)
	//   nd_x[j] = bem_mesh_pt()->node_pt(i_nd)->x(j);

	// //??temp: compare with exact solution
	// // phi_1 = - r * cos(azimuthal) * sin (polar)
	// double r = mod(nd_x);
	// double azimuthal = atan2(nd_x[1],nd_x[0]);
	// double polar = acos(nd_x[2]/r);
	// exact_phi_1[nd_equation] = - r * cos(azimuthal) * sin(polar);


	// if (std::abs( bem_mesh_pt()->node_pt(i_nd)->value(phi_1_index())
	// 	      - exact_phi_1[nd_equation]) > 1e-5)
	//   {
	//     std::cout << nd_x << " " << bem_mesh_pt()->node_pt(i_nd)->value(phi_1_index())
	// 	      <<  " " << exact_phi_1[nd_equation] << std::endl;
	//   }
      }


    // std::cout << boundary_phi_1 << std::endl;

    // Dense matrix multiplication to calculate phi (result goes in boundary_phi)
    DoubleVector boundary_phi(dummy,0.0);
    Boundary_matrix.multiply(boundary_phi_1, boundary_phi);

    // Assign the boundary values of phi from the result of the multiplication
    for(unsigned i_nd=0; i_nd < bem_mesh_pt()->nnode(); i_nd++)
      {
	// Get boundary equation number
	unsigned nd_equation = convert_global_to_boundary_equation_number
	  ( bem_mesh_pt()->node_pt(i_nd)->eqn_number(phi_index()) );

	// Copy out appropriate vector entry
	bem_mesh_pt()->node_pt(i_nd)->set_value(phi_index(), boundary_phi[nd_equation]);

	// //??temp
	// Vector<double> nd_x(DIM,0.0);
	// for(unsigned j=0; j<DIM; j++)
	//   nd_x[j] = bem_mesh_pt()->node_pt(i_nd)->x(j);


	// if( ! small(nd_x[0] + boundary_phi_1[nd_equation]))
	//   std::cout << nd_x << " " << boundary_phi_1[nd_equation] << std::endl;

      }
  }


  //======================================================================
  /// Set up the initial conditions
  //======================================================================
  template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
  void ThreeDHybridProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
  set_initial_condition()
  {

    // Backup time in global Time object
    double backed_up_time=time_pt()->time();

    // Past history needs to be established for t=time0-deltat, ...
    // Then provide current values (at t=time0) which will also form
    // the initial guess for the first solve at t=time0+deltat

    // Get M indicies
    Vector<unsigned> m_index_micromag(3,0);
    BULK_ELEMENT* elem_pt = dynamic_cast< BULK_ELEMENT* >(Bulk_mesh_pt->element_pt(0));
    for(unsigned i=0; i<3; i++)
      m_index_micromag[i] = elem_pt->m_index_micromag(i);

    //Find number of nodes in mesh
    unsigned num_nod = mesh_pt()->nnode();

    // Set continuous times at previous timesteps:
    int nprev_steps=time_stepper_pt()->nprev_values();
    Vector<double> prev_time(nprev_steps+1);
    for (int t=nprev_steps;t>=0;t--)
      {
	prev_time[t]=time_pt()->time(t);
      }

    // Loop over current & previous timesteps
    for (int t=nprev_steps;t>=0;t--)
      {
	// Continuous time
	double time = prev_time[t];
	std::cout << "setting IC at time =" << time << std::endl;

	// Loop over the nodes to set initial values everywhere
	for (unsigned n=0;n<num_nod;n++)
	  {
	    // Get initial value of m from inputs
	    //??ds encapsulate properly
	    Vector<double> m(3,0.0), x(DIM,0.0);
	    mesh_pt()->node_pt(n)->position(t,x);
	    Inputs::initial_m(time,x,m);

	    // Set initial condition on m
	    for(unsigned i=0; i<3; i++)
	      mesh_pt()->node_pt(n)->set_value(t,m_index_micromag[i],m[i]);
	  }
      }

    // Reset backed up time for global timestepper
    time_pt()->time()=backed_up_time;
  }


  //======================================================================
  /// When using mid-point method we must update after each solveto go from
  /// [n+0.5] to [n+1].
  //======================================================================
  template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
  void ThreeDHybridProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
  actions_after_newton_solve()
  {
    // If we are using midpoint then apply the required update (Malidi2005)
    if(Inputs::midpointmethod)
      {
	// Get m indicies
	BULK_ELEMENT* bulk_elem_pt = dynamic_cast<BULK_ELEMENT*>(Bulk_mesh_pt->element_pt(0));
	Vector<unsigned> m_indices(3,0);
	for(unsigned j=0; j<3; j++)
	  m_indices[j] = bulk_elem_pt->m_index_micromag(j);

	for(unsigned i_nd=0; i_nd<mesh_pt()->nnode(); i_nd++)
	  {
	    Node* nd_pt = mesh_pt()->node_pt(i_nd);
	    for(unsigned j=0; j<3; j++)
	      {
		// m[n] --> 2*m[n] - m[n-1]
		// because we have finished this timestep and moved forward one
		double new_m = 2*(nd_pt->value(m_indices[j]))
		  - (nd_pt->value(1,m_indices[j]));
		nd_pt->set_value(m_indices[j],new_m);
	      }
	  }
      }
  }

} // End of oomph namespace

int main(int argc, char* argv[])
{

  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  // Store command line arguments
  CommandLineArgs::setup(argc,argv);

  if(argc < 4)
    {
      throw OomphLibError("Not enough args: needs mesh inputs",
			  "main()", OOMPH_EXCEPTION_LOCATION);
    }

  // Inputs
  const unsigned dim = Inputs::dim;
  const unsigned nnode_1d = 2;
  // const double dt = 1e-3;
  // const unsigned nstep = 10;
  double dt = Inputs::dt;
  const double tmax = Inputs::tmax;

  // Dummy error for timestepper - always be ok
  const double dummy_t_eps = 100;



  // Create the problem
  ThreeDHybridProblem< QMicromagElement <dim,nnode_1d>, MicromagFaceElement, dim >
    problem(argv[1],argv[2],argv[3]);

  // // dump mesh for testing
  // std::ofstream mesh_plot;
  // mesh_plot.open("./mesh_points");
  // for(unsigned nd=0; nd<problem.mesh_pt()->nnode(); nd++)
  //   {
  //     for(unsigned j=0; j<dim; j++)
  // 	mesh_plot << problem.mesh_pt()->node_pt(nd)->x(j) << " ";
  //     mesh_plot << std::endl;
  //   }
  // mesh_plot.close();


  // // dump boundary for testing
  // unsigned b = 0;
  // std::ofstream bound_plot;
  // bound_plot.open("./bound_points");
  // for(unsigned nd=0; nd<problem.Bulk_mesh_pt->nboundary_node(b); nd++)
  //   {
  //     for(unsigned j=0; j<dim; j++)
  // 	bound_plot << problem.Bulk_mesh_pt->boundary_node_pt(b,nd)->x(j) << " ";
  //     bound_plot << std::endl;
  //   }
  // bound_plot.close();


  // Initialise timestep, initial conditions
  problem.initialise_dt(dt);
  problem.set_initial_condition();

  // Set up outputs and output initial conditions
  DocInfo doc_info;
  doc_info.set_directory("results");
  doc_info.number()=0;
  problem.doc_solution(doc_info);
  doc_info.number()++;

  /// Check problem
  if(!(problem.self_test()==0))
    throw OomphLibError("Problem self_test failed","main",
  			OOMPH_EXCEPTION_LOCATION);

  std::cout << "constructor done, everything ready" << "\n" << std::endl;

  std::cout << std::endl;

  // Open a trace file
  std::ofstream trace_file;
  char trace_filename[100];
  sprintf(trace_filename,"%s/trace.dat",doc_info.directory().c_str());
  trace_file.open(trace_filename);


  // std::ofstream residual_file;
  // residual_file.precision(16);
  // char filename2[100];
  // sprintf(filename2,"results/residual");
  // residual_file.open(filename2);
  // residuals.output(residual_file);
  // residual_file.close();

  // std::ofstream bem_file;
  // bem_file.precision(16);
  // char bem_filename[100];
  // sprintf(bem_filename,"results/bem");
  // bem_file.open(bem_filename);
  // problem.boundary_matrix_pt()->output(bem_file);
  // bem_file.close();


  if(Inputs::adaptive_timestepping)
    {
      // Adaptive while loop

      while (problem.time_pt()->time()<tmax)
	{
	  std::cout << "Time is " << problem.time_pt()->time()<< std::endl
		    << "Current timestep is " << dt << std::endl << std::endl;


	  // Take an adaptive timestep -- the input dt is the suggested timestep.
	  // The adaptive timestepper will adjust dt until the required error
	  // tolerance is satisfied. The function returns a suggestion
	  // for the timestep that should be taken for the next step. This
	  // is either the actual timestep taken this time or a larger
	  // value if the solution was found to be "too accurate".
	  double dt_next=problem.adaptive_unsteady_newton_solve(dt,dummy_t_eps);

	  // Use dt_next as suggestion for the next timestep
	  dt=dt_next;

	  //Output solution
	  problem.doc_solution(doc_info);

	  trace_file << doc_info.number() << " " << problem.time_pt()->time()
		     << " " << dt_next << std::endl;

	  //Increment counter for solutions
	  doc_info.number()++;

	} // end of timestepping loop


    }
  else
    {
      unsigned nstep = int(tmax/dt);

      // Standard timestepping loop
      for(unsigned istep=0; istep<nstep; istep++)
	{
	  std::cout << "Timestep " << istep << std::endl;

	  // Take timestep
	  problem.unsteady_newton_solve(dt);

	  //Output solution
	  problem.doc_solution(doc_info);

	  // DenseDoubleMatrix jacobian;
	  // DoubleVector residuals;
	  // problem.get_jacobian(residuals,jacobian);

	  // std::ofstream residual_file;
	  // residual_file.precision(16);
	  // char filename2[100];
	  // sprintf(filename2,"results/residual");
	  // residual_file.open(filename2);
	  // residuals.output(residual_file);
	  // residual_file.close();

	  //Increment counter for solutions
	  doc_info.number()++;
	}

    }
}

#endif