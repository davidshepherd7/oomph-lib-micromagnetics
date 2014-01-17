

#include "boundary_element_handler.h"

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

#endif

  // Get the number of nodes in the boundary problem
  unsigned long n_node = bem_mesh_pt()->nnode();

  // Initialise and resize the boundary matrix ??ds this wastes some time
  // copying over old values. Write a new resize function?
  Bem_matrix.resize(n_node, n_node);
  Bem_matrix.initialise(0.0);

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
                                                       Use_numerical_integration);

      //??ds some of this stuff should be done using iterators to avoid
      // lookups in a map.

      // Loop over the nodes in this element (to copy results into final matrix)
      for(unsigned l=0;l<n_element_node;l++)
        {
          // Get the node number (in the bem mesh) from the global equation number.
          int global_l_number = elem_pt->node_pt(l)->eqn_number(input_index());
          unsigned l_number = input_lookup_pt()->global_to_node(global_l_number);

          // Loop over all nodes in the mesh and add contributions from this element
          for(unsigned long s_nd=0; s_nd<n_node; s_nd++)
            {
              int global_s_number
                = bem_mesh_pt()->node_pt(s_nd)->eqn_number(output_index());
              unsigned s_number = output_lookup_pt()->
                global_to_node(global_s_number);

              // Rows are indexed by output (source node) number, columns
              // are indexed by input (l) number.
              Bem_matrix(s_number, l_number) += element_boundary_matrix(l,s_nd);

              //??ds I think elemental bem matrices are currently actually
              //the transpose, so we needed to switch them back here. Fix this?
            }
        }
    }

  if(!Debug_disable_corner_contributions)
    {
      // Lindholm formula/adaptive integral does not contain the solid angle
      // contribution so add it now.
      corner_list_pt()->add_corner_contributions(Bem_matrix);
    }

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
      unsigned geqn = bem_mesh_pt()->node_pt(nd)->eqn_number(input_index());
      unsigned in_eqn = input_lookup_pt()->global_to_node(geqn);

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
          unsigned g_eqn = m_pt->boundary_node_pt(b,nd)->eqn_number(output_index());
          unsigned out_eqn = output_lookup_pt()->global_to_node(g_eqn);

          (*bem_output_values[i])[nd] = full_vector[out_eqn];
        }
    }

}

void BoundaryElementHandler::
get_bem_values_and_copy_into_values() const
  {
    // Get as one big vector
    DoubleVector full_vector;
    get_bem_values(full_vector);


    // Now split it up into vectors on each boundary
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
            unsigned g_eqn = m_pt->boundary_node_pt(b,nd)->eqn_number(output_index());
            unsigned out_eqn = output_lookup_pt()->global_to_node(g_eqn);

            m_pt->boundary_node_pt(b, nd)->set_value(output_index(),
                                                     full_vector[out_eqn]);
          }
      }
  }

} // End of oomph namespace
