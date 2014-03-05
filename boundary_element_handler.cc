

#include "boundary_element_handler.h"
#include "hmatrix.h"

namespace oomph
{

//==========================================================================
/// Get the fully assembled boundary matrix in dense storage.
//==========================================================================
void BoundaryElementHandler::build_bem_matrix()
{

#ifdef PARANOID
  // Check the corner list has been set up
  if((corner_list_pt() == 0) || !(corner_list_pt()->is_set_up()))
    {
      std::ostringstream error_msg;
      error_msg << "List of the sharp corners of the mesh has not been set up.";
      throw OomphLibError(error_msg.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

  // Check that mesh pointer in elements are equal to mesh pointer here
  for(unsigned e=0, ne=bem_mesh_pt()->nelement(); e < ne; e++)
    {
      MicromagBEMElementEquations* ele_pt =
        checked_dynamic_cast<MicromagBEMElementEquations*>(bem_mesh_pt()->element_pt(e));
      if (ele_pt->boundary_mesh_pt() != bem_mesh_pt())
        {
          std::ostringstream error_msg;
          error_msg << "Wrong mesh pointers in elements.";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }

  if(Bem_matrix_pt != 0)
    {
      std::string err = "Already have a bem matrix, delete before rebuilding.";
      throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                          OOMPH_CURRENT_FUNCTION);
    }
#endif

  // Create the matrix
  DenseDoubleMatrix* dense_matrix_pt = new DenseDoubleMatrix;
  Bem_matrix_pt = dense_matrix_pt;


  // Get the number of nodes in the boundary problem
  unsigned long n_node = bem_mesh_pt()->nnode();

  // Initialise and resize the boundary matrix ??ds this wastes some time
  // copying over old values. Write a new resize function?
  dense_matrix_pt->resize(n_node, n_node);
  dense_matrix_pt->initialise(0.0);

  // Loop over all elements in the BEM mesh
  unsigned long n_bem_element = bem_mesh_pt()->nelement();
  for(unsigned long e=0;e<n_bem_element;e++)
    {
      // Get the pointer to the element and cast it.
      MicromagBEMElementEquations* elem_pt =
        checked_dynamic_cast<MicromagBEMElementEquations*>
        (bem_mesh_pt()->element_pt(e));

      // Find number of nodes in the element
      unsigned long n_element_node = elem_pt->nnode();

      // Set up and initialise matrix
      DenseMatrix<double> element_boundary_matrix(n_element_node,n_node,0.0);

      // Fill the matrix
      elem_pt->fill_in_contribution_to_boundary_matrix(element_boundary_matrix,
                                                       Numerical_int_bem);

      //??ds some of this stuff should be done using iterators to avoid
      // lookups in a map.

      // Loop over the nodes in this element (to copy results into final matrix)
      for(unsigned l=0;l<n_element_node;l++)
        {
          // Get the node number (in the bem mesh) from the global equation number.
          unsigned l_number = input_equation_number(elem_pt->node_pt(l));

          // Loop over all nodes in the mesh and add contributions from this element
          for(unsigned long s_nd=0; s_nd<n_node; s_nd++)
            {
              unsigned s_number
                = output_equation_number(bem_mesh_pt()->node_pt(s_nd));

              // Rows are indexed by output (source node) number, columns
              // are indexed by input (l) number.
              dense_matrix_pt->operator()(s_number, l_number)
                += element_boundary_matrix(l,s_nd);

              //??ds I think elemental bem matrices are currently actually
              //the transpose, so we needed to switch them back here. Fix this?
            }
        }
    }

  if(!Debug_disable_corner_contributions)
    {
      // Lindholm formula/adaptive integral does not contain the solid angle
      // contribution so add it now.
      corner_list_pt()->add_corner_contributions(*dense_matrix_pt);
    }

  // dense_matrix_pt->output("dense_oomph_mat");

}

/// Build the bem_matrix as a sum of the hierarchical matrix and the corner
/// contributions. If hlib is not installed and linked then throw an error
/// instead.
void BoundaryElementHandler::build_hierarchical_bem_matrix()
{
#ifdef OOMPH_HAS_HLIB
  // Note: we can't include the corner contributions directly into the
  // H-matrix because it would no longer be purely a discretised double
  // layer potential, and hence some of the hierarchical properties
  // may/would be lost.

  // Generate the hierarchical part of the bem matrix.
  HMatrix* h_matrix_pt = new HMatrix;
  h_matrix_pt->build(*Bem_mesh_pt, this);

  // Put h matrix in total matrix, set it to be deleted when the total
  // matrix is.
  SumOfMatrices* total_bem_matrix_pt = new SumOfMatrices;
  total_bem_matrix_pt->main_matrix_pt() = h_matrix_pt;
  total_bem_matrix_pt->set_delete_main_matrix();

  // Create the lookup scheme between hlib indicies and our numbering
  // scheme.
  Hmatrix_dof2idx_mapping_pt =
    new AddedMainNumberingLookup(h_matrix_pt->cluster_tree_pt()->dof2idx,
                                 h_matrix_pt->nrow());

  // Put matrix containing the bem sharp corner contributions into the
  // total matrix, also set it to be deleted when the main matrix is. The
  // mapping between the two is the H matrix's dof2idx mapping.
  CRDoubleMatrix* corners_matrix_pt = new CRDoubleMatrix;
  corner_list_pt()->make_diagonal_corner_matrix(*corners_matrix_pt);
  total_bem_matrix_pt->add_matrix(corners_matrix_pt,
                                  Hmatrix_dof2idx_mapping_pt,
                                  Hmatrix_dof2idx_mapping_pt,
                                  true);

  // Set the total as the bem matrix
  Bem_matrix_pt = total_bem_matrix_pt;
#else
  throw OomphLibError("Cannot build H matrix without Hlib library installed",
                      OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
#endif
}

//======================================================================
/// Build the mesh of bem elements.
//======================================================================
void BoundaryElementHandler::build_bem_mesh()
{
  build_bem_mesh_helper(Bem_boundaries, Bem_element_factory_fpt,
                        integration_scheme_pt(), *Bem_mesh_pt);
}

// =================================================================
/// If the boundary matrix is distributed then get its
/// distribution. Otherwise return a "non-distributed distribution" with
/// the correct number of rows.
// =================================================================
void BoundaryElementHandler::
get_bm_distribution(LinearAlgebraDistribution& dist) const
{
  // Try to cast to a distributed object.
  const DistributableLinearAlgebraObject* dist_bm_pt =
    dynamic_cast<const DistributableLinearAlgebraObject* >
    (bem_matrix_pt());

  // If it's not distributable (i.e. if the cast failed) then make a dummy
  // one, otherwise copy the distribution.
  if(dist_bm_pt == 0)
    {
      dist.build(MPI_Helpers::communicator_pt(),
                 bem_matrix_pt()->nrow(), false);
    }
  else
    {
      dist.build(dist_bm_pt->distribution_pt());
    }
}


// =================================================================
/// Put the output values from the boundary element method into a
/// DoubleVector.
// =================================================================
void BoundaryElementHandler::
get_bem_values(DoubleVector &bem_output_values) const
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
      unsigned in_eqn = input_equation_number(bem_mesh_pt()->node_pt(nd));
      input_values[in_eqn] = bem_mesh_pt()->node_pt(nd)->value(input_index());
    }

  // Matrix multiply to get output values
  bem_matrix_pt()->multiply(input_values, bem_output_values);

}

// =================================================================
/// Put the output values from the boundary element method into vectors
/// (one per boundary). Not exceptionally efficient....
// =================================================================
void BoundaryElementHandler::
get_bem_values(const Vector<DoubleVector*> &bem_output_values) const
{

  // Get as one big vector
  DoubleVector full_vector;
  get_bem_values(full_vector);

#ifdef PARANOID
  if(bem_output_values.size() != Bem_boundaries.size())
    {
      std::string error_msg = "Wrong number of doublevectors in output vector.";
      throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif


  // Now split it up into vectors on each boundary
  for(unsigned i=0, ni=Bem_boundaries.size(); i < ni; i++)
    {
      // Get info on this boundary
      unsigned b = Bem_boundaries[i].first;
      const Mesh* m_pt = Bem_boundaries[i].second;
      unsigned nnode = m_pt->nboundary_node(b);

      // Check output DoubleVector is the correct size
#ifdef PARANOID
      if(bem_output_values[i]->nrow() != nnode)
        {
          std::string error_msg = "Output DoubleVector " + to_string(i)
            + " is the wrong size.";
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      // Fill it in
      for(unsigned nd=0; nd<nnode; nd++)
        {
          unsigned out_eqn = output_equation_number(m_pt->boundary_node_pt(b,nd));
          (*bem_output_values[i])[nd] = full_vector[out_eqn];
        }
    }

}

void BoundaryElementHandler::get_bem_values_and_copy_into_values() const
  {
    // Get as one big vector
    DoubleVector full_vector;
    get_bem_values(full_vector);


    for(unsigned i=0, ni=Bem_boundaries.size(); i < ni; i++)
      {
        // Get info on this boundary
        unsigned b = Bem_boundaries[i].first;
        const Mesh* m_pt = Bem_boundaries[i].second;
        unsigned nnode = m_pt->nboundary_node(b);

        // Set the entry corresponding to output_value_index in this node
        // to the corresponding value from the doublevector.
        for(unsigned nd=0; nd<nnode; nd++)
          {
            unsigned out_eqn = output_equation_number(m_pt->boundary_node_pt(b,nd));
            m_pt->boundary_node_pt(b, nd)->set_value(output_index(),
                                                     full_vector[out_eqn]);
          }
      }
  }

/// If we are using H-matrix for bem then write out some data on it.
void BoundaryElementHandler::maybe_write_h_matrix_data(const std::string& outdir) const
{
#ifdef OOMPH_HAS_HLIB
  if(Hierarchical_bem)
    {
      // Fun with types...
      SumOfMatrices* sum_pt = checked_dynamic_cast<SumOfMatrices*>(Bem_matrix_pt);
      HMatrix* hm_pt = checked_dynamic_cast<HMatrix*>(sum_pt->main_matrix_pt());

      // Dump the data
      std::string rank_filename = outdir + "/h_matrix_rank.ps";
      outputrank_supermatrix(hm_pt->supermatrix_pt(), rank_filename.c_str());
    }
  // Else do nothing
#endif
  // Do nothing if we don't have hlib
}

} // End of oomph namespace
