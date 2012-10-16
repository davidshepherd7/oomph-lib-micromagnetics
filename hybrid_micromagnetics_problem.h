#ifndef OOMPH_MICROMAGNETICS_PROBLEM_H
#define OOMPH_MICROMAGNETICS_PROBLEM_H


/*
  description of file goes here
*/

#include "generic.h"
#include "./magnetic_materials.h"
#include "./my_general_header.h"

#include "./boundary_element_handler.h"

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
   : Flux_mesh_pt(0),
     Bulk_mesh_pt(0), Bulk_mesh_magnetic_parameters_pt(0)
  {
   debug_doc().disable_doc();
   debug_doc().set_directory("debug");
  }

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

  /// ??ds
  void create_flux_elements(const unsigned& b);

  /// ??ds
  void get_residuals(DoubleVector& residuals)
  {
   //std::cout << "Calling your get_residuals" << std::endl;
   Problem::get_residuals(residuals);

   //??dsparallel need to override the parallelism on residuals somehow...
   // alternatively we could properly parallelise dense matrix
   // multiplication.
   insert_bem_phi_residual_contribution(residuals);
  }

  /// ??ds
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

  /// ??ds
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

  /// Fill in residuals that the hybrid method part is responsible for.
  void insert_bem_phi_residual_contribution(DoubleVector& residuals) const
  {
   // Get the bem output
   DoubleVector bem_values;
   bem_handler_pt()->get_bem_values(bem_values);

   // Insert into residuals
   for(unsigned nd=0; nd< bem_mesh_pt()->nnode(); nd++)
    {
     Node* nd_pt = bem_mesh_pt()->node_pt(nd);

     // Calculate residual: difference between bem value and actual value.
     double r = bem_values[nd] - nd_pt->value(phi_index());

     // Put into residual in position. Note the residual is entirely
     // determined by a single node (rather than integrated over element)
     // so we can use =.
     residuals[nd_pt->eqn_number(phi_index())] = r;
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
    unsigned bem_size = bem_mesh_pt()->nnode();
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

    // // Set up iterators and maps
    // std::map<long unsigned, long unsigned>::const_iterator row_it, col_it;
    // const std::map<long unsigned, long unsigned> *phi_map_pt, *phi1_map_pt;
    // phi_map_pt = bem_handler_pt()->global_to_boundary_equation_map_pt();
    // phi1_map_pt = bem_handler_pt()->global_phi_1_num_map_pt();

    // // Add the boundary matrix entries to the list
    // for(row_it = phi_map_pt->begin(); row_it != phi_map_pt->end(); row_it++)
    //  {
    //   for(col_it = phi1_map_pt->begin(); col_it != phi1_map_pt->end();
    //       col_it++)
    //    {

    // Double loop over all bem nodes - once for phi once for phi1.
    for(unsigned pnd=0; pnd< bem_mesh_pt()->nnode(); pnd++)
     {
      for(unsigned p1nd=0; p1nd< bem_mesh_pt()->nnode(); p1nd++)
       {
	int global_row = bem_mesh_pt()->node_pt(pnd)->eqn_number(phi_index());
	int global_col = bem_mesh_pt()->node_pt(p1nd)->eqn_number(phi_1_index());
	double value = bem_handler_pt()-> boundary_matrix_pt()->operator()(pnd,p1nd);

	RowColVal entry(global_row,global_col,value);
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
     jacobian.sparse_indexed_output(jac_filename);
     sparse_jacobian.sparse_indexed_output(spjac_filename);
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

   sparse_jacobian_pt->sparse_indexed_output("pre-change");

   // Overwrite the dphi_dphi boundary block with -1*I (sparse BEM part)
   overwrite_bem_sparse_block(sparse_jacobian_pt);

   sparse_jacobian_pt->sparse_indexed_output("post-change");

   // Add the boundary element matrix to the total Jacobian. It represents
   // the derivative of phi with respect to phi_1 so each entry goes in the
   // phi row and the phi_1 column of the respective element (this is done
   // via the two lookup schemes). 0 says don't delete when done.
   jacobian.add_matrix(boundary_matrix_pt(),
		       bem_handler_pt()->row_lookup_pt(),
		       bem_handler_pt()->col_lookup_pt(),
		       0);

   //  std::cout << *(bem_handler_pt()->col_lookup_pt()->node_to_global_mapping_pt()) << std::endl;
   // std::cout << *(bem_handler_pt()->col_lookup_pt()->global_to_node_mapping_pt()) << std::endl;

   // jacobian.sparse_indexed_output("test");

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
     jacobian.sparse_indexed_output(jac_filename);
     sparse_jacobian_pt->sparse_indexed_output(spjac_filename);
     residuals.output(res_filename);
    }

  }

  /// ??ds
  void overwrite_bem_sparse_block(CRDoubleMatrix* const sparse_jacobian_pt) const;

  // Access functions:
  // ============================================================

  /// ??ds
  Mesh* bulk_mesh_pt() const {return Bulk_mesh_pt;}

  /// ??ds
  Mesh*& bulk_mesh_pt() {return Bulk_mesh_pt;}

  /// ??ds
  MyDocInfo& debug_doc() const {return Debug_doc;}

  /// ??ds
  MyDocInfo& debug_doc() {return Debug_doc;}

  /// ??ds
  MagneticParameters* magnetic_parameters_pt() const
  {return Bulk_mesh_magnetic_parameters_pt;}

  /// ??ds
  MagneticParameters*& magnetic_parameters_pt()
  {return Bulk_mesh_magnetic_parameters_pt;}

  /// Get the index of phi for use in BEM mapping
  unsigned phi_index() const {return Phi_index;}

  /// Get the index of phi_1 for use in BEM mapping
  unsigned phi_1_index() const {return Phi_1_index;}

  /// ??ds
  unsigned m_index(const unsigned& i) const {return M_index[i];}

  /// Get the flux mesh pointer.
  Mesh* flux_mesh_pt() const {return Flux_mesh_pt;}

  /// Get a pointer to the bem handler
  const BoundaryElementHandler<BEM_ELEMENT<BULK_ELEMENT, DIM> >*
  bem_handler_pt() const {return &Bem_handler;}

  /// Get a non-const pointer to the bem handler
  BoundaryElementHandler<BEM_ELEMENT<BULK_ELEMENT, DIM> >*
  bem_handler_pt() {return &Bem_handler;}

  /// Get mesh from bem handler
  const Mesh* bem_mesh_pt() const {return bem_handler_pt()->bem_mesh_pt();}

  /// Non-const access to bem matrix (const can be got via bem_handler if needed).
  DoubleMatrixBase* boundary_matrix_pt() {return Bem_handler.boundary_matrix_pt();}

 private:

  /// ??ds
  BoundaryElementHandler<BEM_ELEMENT<BULK_ELEMENT, DIM> >
  Bem_handler;

  /// ??ds
  Mesh* Flux_mesh_pt;

  /// ??ds
  Mesh* Bulk_mesh_pt;

  /// Store the magnetic parameters of the bulk mesh region.
  //??ds eventually may need to be able to store multiple meshes with
  // different parameters
  MagneticParameters* Bulk_mesh_magnetic_parameters_pt;

  /// Store the index of phi for use in BEM mapping
  unsigned Phi_index;

  /// Store the index of phi_1 for use in BEM mapping
  unsigned Phi_1_index;

  /// ??ds
  Vector<unsigned> M_index;

  /// ??ds
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

 // Get the indicies for phi, phi_1 and m via casting a pointer
 BULK_ELEMENT* elem_pt = dynamic_cast< BULK_ELEMENT* >(bulk_mesh_pt()->element_pt(0));
 Phi_index = elem_pt->phi_index_micromag();
 Phi_1_index = elem_pt->phi_1_index_micromag();
 M_index.assign(3,0);
 for(unsigned j=0; j<3; j++) M_index[j] = elem_pt->m_index_micromag(j);


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

 // Build global finite element mesh
 add_sub_mesh(bulk_mesh_pt());
 add_sub_mesh(flux_mesh_pt());
 // add_sub_mesh(bem_mesh_pt());
 build_global_mesh();

 // Setup equation numbering scheme for all the finite elements
 std::cout << "FEM number of equations: " << assign_eqn_numbers() << std::endl;


 // Set up the boundary element method
 // ============================================================

 // Set all boundaries of bulk mesh to use bem
 Bem_handler.set_bem_all_boundaries(bulk_mesh_pt());

 // Tell it the dof numbering
 Bem_handler.input_index() = phi_1_index();
 Bem_handler.output_index() = phi_index();

 // Build everything needed by BEM
 Bem_handler.build();
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

// ============================================================
/// Write a -I block into the sparse Jacobian in the dphi/dphi_1 block. The
/// way this is done is a little hacky because we are working with
/// CRDoubleMatrices...
// ============================================================
template<class BULK_ELEMENT, template<class,unsigned> class BEM_ELEMENT, unsigned DIM>
void HybridMicromagneticsProblem<BULK_ELEMENT,BEM_ELEMENT,DIM>::
overwrite_bem_sparse_block(CRDoubleMatrix* const sparse_jac_pt) const
{
 for(unsigned nd=0; nd< bem_mesh_pt()->nnode(); nd++)
  {
   // Find where this row is in the CRDoubleMatrix
   int jacobian_row = bem_handler_pt()->row_lookup_pt()->node_to_global(nd);
   long unsigned row_start = sparse_jac_pt->row_start()[jacobian_row];
   long unsigned row_end = sparse_jac_pt->row_start()[jacobian_row+1];

   for(unsigned i=row_start; i<row_end; i++)
    {
#ifdef PARANOID
     // Sort of check if anything other than the dummy value has been added
     // to the entry by checking if there is any remainder. Not perfect
     // because it is possible (but unlikely) that we could have no
     // remainder by fluke.
     if(!(std::fmod(sparse_jac_pt->value()[i],
		    MicromagEquations<DIM>::DummyBEMControlledEntry)
	  < 1e-10))
      {
       std::ostringstream error_msg;
       error_msg << "Trying to overwrite a value to which something"
		 << " else has (probably) been added.";
       throw OomphLibError(error_msg.str(),
			   "HybridMicromagneticsProblem::overwrite_bem_sparse_block",
			   OOMPH_EXCEPTION_LOCATION);
      }
#endif

     // Overwrite the diagonal with -1. The diagonal should already exist
     // because of some nasty hacks in the elemental get_jacobian
     // functions.
     if(sparse_jac_pt->column_index()[i] == jacobian_row)
      {
       sparse_jac_pt->value()[i] = -1;
      }
     else
      {
       // There should only be diagonal entries in these rows, but if there
       // is anything else (maybe we messed with the jacobian generating
       // functions?) then overwrite it with a zero and give a warning.
       sparse_jac_pt->value()[i] = 0;
       std::cerr << "Overwrote non-diagonal values, this is ok if you are"
		 << " FD-ing the jacobian but otherwise might be a problem."
		 << std::endl;
      }
    }
  }
}


} // End of oomph namespace


#endif
