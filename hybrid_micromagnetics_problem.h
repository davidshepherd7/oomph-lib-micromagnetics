#ifndef OOMPH_MICROMAGNETICS_PROBLEM_H
#define OOMPH_MICROMAGNETICS_PROBLEM_H


/*
  description of file goes here
*/

#include "generic.h"
#include "./magnetic_materials.h"
#include "./my_general_header.h"

using namespace oomph;

namespace oomph
{

  class MyDocInfo;

  template<class BULK_ELEMENT,
	   template<class BULK_ELEMENT,unsigned DIM> class BEM_ELEMENT,
	   unsigned DIM>
  class HybridMicromagneticsProblem : public Problem
  {

  public:

    HybridMicromagneticsProblem()
      : Flux_mesh_pt(0), Bem_mesh_pt(0),
	Bulk_mesh_pt(0), Bulk_mesh_magnetic_parameters_pt(0)
    {
      debug_doc().disable_doc();
      debug_doc().set_directory("debug");
    };

    /// Blank destructor -
    ~HybridMicromagneticsProblem() {}

    virtual void doc_solution(DocInfo& doc_info) = 0;

    virtual void set_initial_condition() = 0;

    void actions_after_newton_step()
    { debug_doc().next_newton_step(); }

    //??ds If we are using mid-point then update, if we are using bdf then
    // re-normalise
    void actions_after_newton_solve()
    {
      // renormalise_magnetisation();
    }

    void renormalise_magnetisation()
    {
      //??ds account for actual M_s?
      for(unsigned nd=0; nd<bulk_mesh_pt()->nnode(); nd++)
	{
	  Node* nd_pt = bulk_mesh_pt()->node_pt(nd);

	  Vector<double> m_values(3,0.0);
	  for(unsigned j=0; j<3; j++) m_values[j] = nd_pt->value(m_index(j));
	  VectorOps::normalise(m_values);
	  for(unsigned j=0; j<3; j++) nd_pt->set_value(m_index(j),m_values[j]);
	}
    }

    /// Finish off construction of problem (once mesh and magnetic properties
    /// have been set up).
    void finish_building_hybrid_problem();

    void build_bem_mesh(Mesh* bem_mesh_pt) const;

    void build_boundary_matrix();

    void create_flux_elements(const unsigned& b);

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

    std::map<Node*,double>* corners_map_pt() const {return Corners_map_pt;}
    std::map<Node*,double>*& corners_map_pt() {return Corners_map_pt;}

    MyDocInfo& debug_doc() const {return Debug_doc;}
    MyDocInfo& debug_doc() {return Debug_doc;}

    MagneticParameters* magnetic_parameters_pt() const
    {return Bulk_mesh_magnetic_parameters_pt;}

    MagneticParameters*& magnetic_parameters_pt()
    {return Bulk_mesh_magnetic_parameters_pt;}

    /// Get the index of phi for use in BEM mapping
    unsigned phi_index() const {return Phi_index;}

    /// Get the index of phi_1 for use in BEM mapping
    unsigned phi_1_index() const {return Phi_1_index;}

    unsigned m_index(const unsigned& i) const {return M_index[i];}

    /// Get the flux mesh pointer.
    Mesh* flux_mesh_pt() const {return Flux_mesh_pt;}

    /// Access to the boundary matrix
    DenseDoubleMatrix* boundary_matrix_pt()
    {return &Boundary_matrix;}

    void get_residuals(DoubleVector& residuals)
    {
      //std::cout << "Calling your get_residuals" << std::endl;
      Problem::get_residuals(residuals);

      //??dsparallel need to override the parallelism on residuals somehow...
      // alternatively we could properly parallelise dense matrix
      // multiplication.
      insert_bem_phi_residual_contribution(residuals);
    }

    void get_mean_bulk_values(Vector<double>& means)
    {
      means.clear();
      unsigned n_equation_numbers = 5;
      for(unsigned eq_num = 0; eq_num < n_equation_numbers; eq_num++)
	{
	  double sum =0;
	  unsigned n_node = bulk_mesh_pt()->nnode();
	  for(unsigned nd=0; nd< n_node; nd++)
	    {
	      sum += bulk_mesh_pt()->node_pt(nd)->value(eq_num);
	    }
	  means.push_back(sum/n_node);
	}
    }

    Node* get_non_boundary_node(Mesh* mesh_pt)
    {
      for(unsigned nd=0; nd< mesh_pt->nnode(); nd++)
	{
	  Node * node_pt = mesh_pt->node_pt(nd);
	  if(! (node_pt->is_on_boundary()))
	    return node_pt;
	}

      std::ostringstream error_msg;
      error_msg << "No non-boundary nodes in the mesh.";
      throw OomphLibError(error_msg.str(),
			  "HybridMicromagneticsProblem::get_non_boundary_node",
			  OOMPH_EXCEPTION_LOCATION);

      // Never get here!
      return 0;
    }


    // ============================================================
    /// Calculate boundary values for phi from BEM
    // ============================================================
    void insert_bem_phi_residual_contribution(DoubleVector& residuals) const
    {
      // Parallelisation: Do not distribute because dense matrix cannot be
      // distributed. ??dsparallel
      LinearAlgebraDistribution
	dist(communicator_pt(), bem_mesh_pt()->nnode(), false);

      DoubleVector bem_phi_values(dist);
      get_bem_phi_values(bem_phi_values);

      // for each node in bem mesh
      for(unsigned nd=0; nd< bem_mesh_pt()->nnode(); nd++)
      	{
      	  // Get a pointer
      	  Node* nd_pt = bem_mesh_pt()->node_pt(nd);

      	  // get the the bem equation number
      	  unsigned bem_eqn_num = get_boundary_equation_number(nd_pt);
      	  unsigned global_eqn_num = nd_pt->eqn_number(phi_index());

      	  // Insert appropriate value into residuals.
      	  double r = bem_phi_values[bem_eqn_num]
      	    - nd_pt->value(phi_index());

	  // Boundary values of phi are entirely determined by phi_1 so we
	  // use = rather than +=.
      	  residuals[global_eqn_num] = r;
      	}
    }

    /// Overload get_jacobian to include the boundary matrix in a sparse form.
    void get_jacobian(DoubleVector& residuals, CRDoubleMatrix& jacobian)
    {
      std::cout << "Calling your get_jacobian function with bem added directly." << std::endl;

      // Get the fem jacobian (in the same distribution pattern as original).
      CRDoubleMatrix sparse_jacobian(jacobian.distribution_pt());
      Problem::get_jacobian(residuals,sparse_jacobian);

      // Finish off the residual calculation
      insert_bem_phi_residual_contribution(residuals);

      // Finish off sparse jacobian
      overwrite_bem_sparse_block(&sparse_jacobian);

      // Storage for jacobian data
      Vector<int> row_starts, cols;
      Vector<double> vals;
      {
	// Storage for values of total Jacobian
	Vector<RowColVal> row_col_val;
	unsigned bem_size = Global_to_boundary_equation_map.size();
	unsigned total_nnz = sparse_jacobian.nnz() + bem_size
	  + bem_size*bem_size;
	row_col_val.reserve(total_nnz);

	// Get vector of row indicies for sparse J
	Vector<int> sparse_jac_rowstart(sparse_jacobian.nrow() + 1), sparse_jac_rows;
	std::copy(sparse_jacobian.row_start(),
		  sparse_jacobian.row_start() + sparse_jacobian.nrow() + 1,
		  sparse_jac_rowstart.begin());
	VectorOps::rowstart2rowindex(sparse_jac_rowstart,sparse_jac_rows);

	// Add sparse jacobian part rows/cols/values to our list
	for(unsigned i=0; i<sparse_jacobian.nnz(); i++)
	  {
	    RowColVal entry(sparse_jac_rows[i],
			    *(sparse_jacobian.column_index() + i),
			    *(sparse_jacobian.value() + i));
	    row_col_val.push_back(entry);
	  }

	// Add the boundary matrix entries to the list
	std::map<long unsigned, long unsigned>::
	  const_iterator row_it, col_it;
	for(row_it = Global_to_boundary_equation_map.begin();
	    row_it != Global_to_boundary_equation_map.end();
	    row_it++)
	  {
	    for(col_it = Global_phi_1_num_map.begin();
		col_it != Global_phi_1_num_map.end();
		col_it++)
	      {
		int this_row = row_it->first;
		int this_col = col_it->first;
		double this_value = boundary_matrix_pt()->
		  operator()(row_it->second,col_it->second);
		//double this_value =11234.0;
		RowColVal entry(this_row,this_col,this_value);
		row_col_val.push_back(entry);
	      }

	  }

	// Sort the list ready for making CRDoubleMatrix
	std::sort(row_col_val.begin(), row_col_val.end());

	// Convert back to vectors
	unsigned length = row_col_val.size();
	Vector<int> rows(length,0);
	cols.assign(length,0), vals.assign(length,0.0);
	for(unsigned i=0; i<length; i++)
	  {
	    rows[i] = row_col_val[i].row;
	    cols[i] = row_col_val[i].col;
	    vals[i] = row_col_val[i].val;
	  }

	// Convert rows back to row_start format
	VectorOps::rowindex2rowstart(rows, row_starts);


#ifdef PARANOID
	// Check for duplicates
	for(unsigned j=1; j<cols.size(); j++)
	  {
	    if((cols[j] == cols[j-1]) && ( rows[j] == rows[j-1]))
	      {
		std::cout << "dupe: [" << cols[j] << ", " << rows[j] << ", " << vals[j]
			  << "] and ["
			  << cols[j-1] << ", " << rows[j-1] << ", " << vals[j-1] << "]"
			  << std::endl;
	      }
	  }
#endif
      }

      // Rebuild jacobian with new values
      unsigned ncol = cols.size();
      jacobian.build(sparse_jacobian.distribution_pt(),ncol, vals, cols, row_starts);

      // This is probably very slow...
      if(debug_doc().is_doc_enabled())
    	{
	  // Parallelisation: distribute using the problems communcator,
	  // split rows evenly over all processors. ??dsparallel
	  LinearAlgebraDistribution
	    dist(communicator_pt(), residuals.nrow(), true);
    	  DoubleVector dofs(dist);
    	  get_dofs(dofs);

    	  // Get filenames
    	  char dof_filename[100], jac_filename[100], res_filename[100],
    	    spjac_filename[100];
    	  sprintf(dof_filename, "%s/dofs_%u_%u", debug_doc().directory().c_str(),
    		  debug_doc().timestep(), debug_doc().newton_step());
    	  sprintf(jac_filename, "%s/jac_%u_%u", debug_doc().directory().c_str(),
    		  debug_doc().timestep(), debug_doc().newton_step());
    	  sprintf(res_filename, "%s/res_%u_%u", debug_doc().directory().c_str(),
    		  debug_doc().timestep(), debug_doc().newton_step());
    	  sprintf(spjac_filename, "%s/spjac_%u_%u", debug_doc().directory().c_str(),
    		  debug_doc().timestep(), debug_doc().newton_step());

    	  // Output
    	  dofs.output(dof_filename);
    	  jacobian.Matrix::sparse_indexed_output(jac_filename);
    	  sparse_jacobian.Matrix::sparse_indexed_output(spjac_filename);
    	  residuals.output(res_filename);
    	}
    }


    /// Overload get_jacobian to include the boundary matrix in sumofmatrices
    /// format. Note that this will only work with iterative solvers since we
    /// can only multiply when using sumofmatrices.
    void get_jacobian(DoubleVector& residuals, SumOfMatrices& jacobian)
    {
      std::cout << "Calling your get_jacobian function using SumOfMatrices." << std::endl;

      // Create a matrix to store the sparse part of the Jacobian, put it into
      // the SumOfMatrices. Tell the SumOfMatrices to delete it when we are done
      // with this solve.
      CRDoubleMatrix* sparse_jacobian_pt = new CRDoubleMatrix;
      jacobian.main_matrix_pt() = sparse_jacobian_pt;
      jacobian.set_delete_main_matrix();

      // Fill in values from FEM part
      Problem::get_jacobian(residuals,*sparse_jacobian_pt);


      double extra_stuff_start_time = TimingHelpers::timer();

      // Overwrite the dphi_dphi boundary block with -1*I (sparse BEM part)
      overwrite_bem_sparse_block(sparse_jacobian_pt);

      // Add the boundary element matrix to the total Jacobian. It represents
      // the derivative of phi with respect to phi_1 so each entry goes in the
      // phi row and the phi_1 column of the respective element (this is done
      // via the two maps). 0 says don't delete when done.
      jacobian.add_matrix(boundary_matrix_pt(),
      			  &Global_to_boundary_equation_map,
      			  &Global_phi_1_num_map,
			  0);

      // Output for debugging purposes. This is probably very slow...
      if(debug_doc().is_doc_enabled())
      	{
	  // Parallelisation: distribute using the problems communcator, split
	  // rows evenly over all processors. ??dsparallel
	  LinearAlgebraDistribution
	    dist(communicator_pt(), residuals.nrow(), true);
    	  DoubleVector dofs(dist);
    	  get_dofs(dofs);

      	  // Get filenames
      	  char dof_filename[100], jac_filename[100], res_filename[100],
      	    spjac_filename[100];
      	  sprintf(dof_filename, "%s/dofs_%u_%u", debug_doc().directory().c_str(),
      		  debug_doc().timestep(), debug_doc().newton_step());
      	  sprintf(jac_filename, "%s/jac_%u_%u", debug_doc().directory().c_str(),
      		  debug_doc().timestep(), debug_doc().newton_step());
      	  sprintf(res_filename, "%s/res_%u_%u", debug_doc().directory().c_str(),
      		  debug_doc().timestep(), debug_doc().newton_step());
      	  sprintf(spjac_filename, "%s/spjac_%u_%u", debug_doc().directory().c_str(),
      		  debug_doc().timestep(), debug_doc().newton_step());

      	  // Output
      	  dofs.output(dof_filename);
      	  jacobian.Matrix::sparse_indexed_output(jac_filename);
      	  sparse_jacobian_pt->Matrix::sparse_indexed_output(spjac_filename);
      	  residuals.output(res_filename);
      	}

      double extra_stuff_stop_time = TimingHelpers::timer();
      std::cout << "Set up time increase due to BEM and/or debug output was " <<
	extra_stuff_stop_time - extra_stuff_start_time << std::endl;
    }

    void overwrite_bem_sparse_block(CRDoubleMatrix* const sparse_jacobian_pt) const;

  private:

    /// The map between the global equation numbers for phi and the boundary
    /// equation/matrix numbering.
    std::map<long unsigned,long unsigned> Global_to_boundary_equation_map;

    /// The map between the global equation numbers for phi_1 and the boundary
    /// equation/matrix numbering.
    std::map<long unsigned,long unsigned> Global_phi_1_num_map;

    Mesh* Flux_mesh_pt;

    /// The pointer to the boundary element method mesh
    Mesh* Bem_mesh_pt;

    ///
    Mesh* Bulk_mesh_pt;

    /// A map containing node pointers to nodes which are on sharp corners and
    /// the angle of their corner.
    std::map<Node*,double>* Corners_map_pt;

    /// Store the magnetic parameters of the bulk mesh region.
    //??ds eventually may need to be able to store multiple meshes with
    // different parameters
    MagneticParameters* Bulk_mesh_magnetic_parameters_pt;

    /// Store the index of phi for use in BEM mapping
    unsigned Phi_index;

    /// Store the index of phi_1 for use in BEM mapping
    unsigned Phi_1_index;

    /// Store the
    Vector<unsigned> M_index;

    /// Matrix to store the relationship between phi_1 and phi on the boundary
    DenseDoubleMatrix Boundary_matrix;

    MyDocInfo Debug_doc;

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

  M_index.assign(3,0);
  for(unsigned j=0; j<3; j++)
    M_index[j] = elem_pt->m_index_micromag(j);

  // Pin a value of phi_1 to zero somewhere in the bulk_mesh.
  // ============================================================

  //This is necessary to avoid having a free constant of integration (which
  // causes scaling problems). Just get the first boundary node ??ds not sure
  // this is ok...
  Node* pinned_phi_1_node_pt = get_non_boundary_node(bulk_mesh_pt());
  pinned_phi_1_node_pt->pin(Phi_1_index);
  pinned_phi_1_node_pt->set_value(Phi_1_index,0.0);

  // Create flux elements (in a seperate mesh so that block preconditioning
  // works).
  // ============================================================

  // We want Neumann (flux) boundary condtions on phi_1 on all boundaries, so
  // create the face elements needed.
  Flux_mesh_pt = new Mesh;
  for(unsigned b=0; b < bulk_mesh_pt()->nboundary(); b++)
    {
      create_flux_elements(b);
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
  add_sub_mesh(flux_mesh_pt());
  // add_sub_mesh(bem_mesh_pt());
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

  //std::cout << Global_to_boundary_equation_map << std::endl;

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
	      // I think the sign here is negative because the lindholm formula
	      // is for +ve but our final equation has negative kernel...
	    }
	}
    }

#ifdef PARANOID
  if(corners_map_pt() == 0)
    {
      std::ostringstream error_msg;
      error_msg << "Map listing sharp corners is not set up.";
      throw OomphLibError(error_msg.str(),
			  "HybridMicromagneticsProblem::build_boundary_matrix()",
			  OOMPH_EXCEPTION_LOCATION);
    }
#endif

  // Lindholm formula does not contain the solid angle contribution so add it
  // here: loop over the matrix diagonals adding the angle factor.
  for(unsigned long nd = 0; nd < bem_mesh_pt()->nnode(); nd++)
    {
      // Check if the node is in the map of corner nodes
      Node* node_pt = bem_mesh_pt()->node_pt(nd);
      std::map<Node*,double>::const_iterator it = corners_map_pt()->find(node_pt);
      if(it != corners_map_pt()->end())
	{
	  // Add the fractional angle for this node
	  Boundary_matrix(nd,nd) += it->second;
	}
      else
	{
	  // This accounts for points which are smooth (at least in the limit
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
  Global_to_boundary_equation_map.clear();
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

#ifdef PARANOID
      //Problems could occur if the index we are using ever has pinned values
      if((global_phi_number < 0) || (global_phi_1_number < 0))
	{
	  //std::cout << Global_to_boundary_equation_map << std::endl;
	  throw OomphLibError
	    ("Pinned equation found in one of the boundary phi values, this destroys the numbering system used to map between the boundary and finite element methods.",
	     "HybridMicromagneticsProblem::get_boundary_equation_number",
	     OOMPH_EXCEPTION_LOCATION);
	}
#endif

      // Set up the pair ready to input with key="global equation number" and
      // value ="boundary equation number"=k.
      std::pair<unsigned,unsigned> input_pair_phi
	= std::make_pair(global_phi_number,k_phi);
      std::pair<unsigned,unsigned> input_pair_phi_1
	= std::make_pair(global_phi_1_number,k_phi_1);

      // Add entry to map and store whether this was a new addition
      bool new_addition_phi = (Global_to_boundary_equation_map.insert(input_pair_phi)
			       ).second;
      bool new_addition_phi_1 = (Global_phi_1_num_map.insert(input_pair_phi_1)
				 ).second;

      // Increment k if this was a new addition to the map
      if(new_addition_phi) k_phi++;
      if(new_addition_phi_1) k_phi_1++;
    }

  //std::cout << Global_to_boundary_equation_map << std::endl;
  // std::cout << Global_phi_1_num_map << std::endl;
}

//======================================================================
/// Given a pointer to a boundary node find the corresponding row/column in the
/// BEM matrix.
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
  if(Global_to_boundary_equation_map.find(global_equation_number)
     == Global_to_boundary_equation_map.end())
    {
      std::ostringstream error_stream;
      error_stream << "Global equation number " << global_equation_number
		   << " is not in the global to boundary map.";
      throw OomphLibError(error_stream.str(),
			  "HybridMicromagneticsProblem::get_boundary_equation_number",
			  OOMPH_EXCEPTION_LOCATION);
    }
#endif

  return Global_to_boundary_equation_map.find(global_equation_number)->second;
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
create_flux_elements(const unsigned& b)
{
  // How many bulk elements are adjacent to boundary b?
  unsigned n_element = bulk_mesh_pt()->nboundary_element(b);

  // Loop over the bulk elements adjacent to boundary b
  for(unsigned e=0;e<n_element;e++)
    {
      // Get pointer to the bulk element that is adjacent to boundary b
      BULK_ELEMENT* bulk_elem_pt = dynamic_cast<BULK_ELEMENT*>
	(bulk_mesh_pt()->boundary_element_pt(b,e));

      // What is the index of the face of the bulk element at the boundary
      int face_index = bulk_mesh_pt()->face_index_at_boundary(b,e);

      // Build the corresponding prescribed-flux element
      MicromagFluxElement<BULK_ELEMENT>* flux_element_pt =
	new MicromagFluxElement<BULK_ELEMENT>(bulk_elem_pt,face_index);

      // Pass a pointer to the flux element to the bulk element
      bulk_elem_pt->add_face_element_pt(flux_element_pt);

      // Add the prescribed-flux element to the mesh
      flux_mesh_pt()->add_element_pt(flux_element_pt);

    } // End of loop over bulk elements adjacent to boundary b
}

template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
void HybridMicromagneticsProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
get_boundary_phi_1(DoubleVector& boundary_phi_1) const
{
#ifdef PARANOID
  if(!(boundary_phi_1.built()))
    {
      std::ostringstream error_msg;
      error_msg << "Distribution should be set up for boundary_phi_1 before"
		<< " passing it into this function.";
      throw OomphLibError(error_msg.str(),
			  "HybridMicromagneticsProblem::get_boundary_phi_1",
			  OOMPH_EXCEPTION_LOCATION);
    }
#endif

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
  // Parallelisation: Do not distribute because dense matrix cannot be
  // distributed. ??dsparallel
  LinearAlgebraDistribution
    dist(communicator_pt(), bem_mesh_pt()->nnode(), false);

  // Assemble a vector of phi_1 values on boundary nodes
  DoubleVector boundary_phi_1(dist);
  get_boundary_phi_1(boundary_phi_1);

  // Dense matrix multiplication to calculate phi (result goes in Bem_phi_values)
  Boundary_matrix.multiply(boundary_phi_1, bem_phi_values);
}

// ============================================================
///
// ============================================================
template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
void HybridMicromagneticsProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
overwrite_bem_sparse_block(CRDoubleMatrix* const sj_pt) const
{
  std::map<long unsigned, long unsigned>::const_iterator it;
  for(it = Global_to_boundary_equation_map.begin();
      it != Global_to_boundary_equation_map.end();
      it++)
    {
      // Find where this row is in the CRDoubleMatrix
      long int jacobian_row = it->first;
      long unsigned row_start = sj_pt->row_start()[jacobian_row];
      long unsigned row_end = sj_pt->row_start()[jacobian_row+1];

      for(unsigned i=row_start; i<row_end; i++)
	{
	  // Sort of check if anything other than the dummy value has been added
	  // to the entry by checking if there is any remainder. Not perfect
	  // because it is possible (but unlikely) that we could have no
	  // remainder by fluke.
	  if(!(std::fmod(sj_pt->value()[i],
			 MicromagEquations<DIM>::DummyBEMControlledEntry)
	       < 1e-10))
	    {
	      std::ostringstream error_msg;
	      error_msg << "Trying to overwrite a value to which something else has (probably) been added.";
	      throw OomphLibError(error_msg.str(),
				  "HybridMicromagneticsProblem::overwrite_bem_sparse_block",
				  OOMPH_EXCEPTION_LOCATION);
	    }


	  if (sj_pt->column_index()[i] == jacobian_row)
	    sj_pt->value()[i] = -1;
	  else
	    {
	      sj_pt->value()[i] = 0;
	      std::cerr << "Overwrote non-diagonal values, this is ok if you are FD-ing the jacobian but otherwise might be a problem." << std::endl;
	    }
	}
    }
}


} // End of oomph namespace


#endif
