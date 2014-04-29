

#include "pinned_boundary_element_handler.h"
#include "hmatrix.h"

namespace oomph
{

//==========================================================================
/// Get the fully assembled boundary matrix in dense storage.
//==========================================================================
void PinnedBoundaryElementHandler::build_bem_matrix()
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


  if(Pinned_bem_matrix_pt != 0)
    {
      std::string err = "Already have a pinned bem matrix, delete before rebuilding.";
      throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                          OOMPH_CURRENT_FUNCTION);
    }
#endif

  // Create the matrices, store pointers to dense matrices (class variables
  // are DoubleMatrixBase*).
  DenseDoubleMatrix* d_bem_matrix_pt = new DenseDoubleMatrix;
  DenseDoubleMatrix* d_pinned_bem_matrix_pt = new DenseDoubleMatrix;
  Bem_matrix_pt = d_bem_matrix_pt;
  Pinned_bem_matrix_pt = d_pinned_bem_matrix_pt;



  // Get the number of nodes in the boundary problem
  unsigned long nunpinned = Lookup_unpinned_input.size();
  unsigned long npinned = Lookup_pinned_input.size();
  unsigned long nunpinned_output = Lookup_unpinned_output.size();

  // Initialise and resize the boundary matrices
  d_bem_matrix_pt->resize(nunpinned_output, nunpinned);
  d_bem_matrix_pt->initialise(0.0);
  d_pinned_bem_matrix_pt->resize(nunpinned_output, npinned);
  d_pinned_bem_matrix_pt->initialise(0.0);

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

      // Set up and initialise the element bem matrix
      DenseMatrix<double> element_boundary_matrix
        (n_element_node, bem_mesh_pt()->nnode(), 0.0);

      // Fill the matrix
      elem_pt->fill_in_contribution_to_boundary_matrix(element_boundary_matrix,
                                                       Numerical_int_bem);

      //??ds optimise using iterators instead of lookups?

      // Loop over the nodes in this element with unpinned input dof (to copy results into final matrix)
      for(unsigned l=0;l<n_element_node;l++)
        {
          // Get the bem equation number from the node pt.
          unsigned l_bemeq = 0;
          if(elem_pt->node_pt(l)->is_pinned(input_index()))
            {
              l_bemeq = Lookup_pinned_input.node_to_bemeq(elem_pt->node_pt(l));
            }
          else
            {
              l_bemeq = Lookup_unpinned_input.node_to_bemeq(elem_pt->node_pt(l));
            }

          // Loop over all nodes in the mesh with unpinned output dof and
          // add contributions from this element.
          for(unsigned long s_nd=0; s_nd<nunpinned_output; s_nd++)
            {
              // Skip if pinned
              if(elem_pt->node_pt(l)->is_pinned(output_index())) continue;

              unsigned s_bemeq
                = Lookup_unpinned_output.node_to_bemeq(bem_mesh_pt()->node_pt(s_nd));

              // Rows are indexed by output (source node) bemeq, columns
              // are indexed by input (l) bemeq.

              //??ds I think elemental bem matrices are currently actually
              //the transpose, so we needed to switch them back here. Fix this?

              // If input is pinned then it goes in the pinned matrix
              if(elem_pt->node_pt(l)->is_pinned(input_index()))
                {
                  d_pinned_bem_matrix_pt->operator()(s_bemeq, l_bemeq)
                    += element_boundary_matrix(l,s_nd);
                }
              else
                {
                  d_bem_matrix_pt->operator()(s_bemeq, l_bemeq)
                    += element_boundary_matrix(l,s_nd);
                }

            }
        }
    }

//??ds now create pinned matrix: loop over pinned inputs then unpinned outputs

  if(!Debug_disable_corner_contributions)
    {
      // Lindholm formula/adaptive integral does not contain the solid angle
      // contribution so add it now.
      corner_list_pt()->add_corner_contributions(*d_bem_matrix_pt,
                                                 *d_pinned_bem_matrix_pt,
                                                 Lookup_unpinned_input,
                                                 Lookup_pinned_input,
                                                 Lookup_unpinned_output,
                                                 input_index());
    }

  // d_bem_matrix_pt->output("dense_oomph_mat");

}

/// Build the bem_matrix as a sum of the hierarchical matrix and the corner
/// contributions. If hlib is not installed and linked then throw an error
/// instead.
void PinnedBoundaryElementHandler::build_hierarchical_bem_matrix()
{
  throw OomphLibError("Not implemented (yet?).", OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
}


// =================================================================
/// Put the output values from the boundary element method into a
/// DoubleVector. Output is numbered by bem lookup unpinned outputs.
// =================================================================
void PinnedBoundaryElementHandler::
get_bem_values(DoubleVector &bem_output_values) const
{
  // Get the boundary matrix linear algebra distribution (if there is one).
  LinearAlgebraDistribution dist;
  get_bm_distribution(dist);

  // Build output vector
  bem_output_values.build(dist);


  // Extract unpinned phis (inputs)
  unsigned nunpinned = Lookup_unpinned_input.size();
  DoubleVector input_values(LinearAlgebraDistribution(dist.communicator_pt(),
                                                      nunpinned, false));
  for(unsigned nd=0; nd<nunpinned; nd++)
    {
      input_values[nd] = Lookup_unpinned_input.bemeq_to_node(nd)
        ->value(input_index());
    }

  // Matrix multiply to get output values
  bem_matrix_pt()->multiply(input_values, bem_output_values);


  // Extract pinned phis
  const unsigned npinned = Lookup_pinned_input.size();
  DoubleVector pinned_input(LinearAlgebraDistribution(dist.communicator_pt(),
                                                      npinned, false));
  for(unsigned i=0; i<npinned; i++)
    {
      pinned_input[i] = Lookup_pinned_input.bemeq_to_node(i)
        ->value(input_index());
    }

  //multiply to get contribution from pinned input values
  DoubleVector bem_output_values_from_pinned_input(dist);
  Pinned_bem_matrix_pt->multiply(pinned_input,
                                 bem_output_values_from_pinned_input);

  // add the two to get the final result
  bem_output_values += bem_output_values_from_pinned_input;
}

// =================================================================
/// Put the output values from the boundary element method into vectors
/// (one per boundary). Not exceptionally efficient....
// =================================================================
void PinnedBoundaryElementHandler::
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

      //??ds use iterator?

      // Fill it in
      for(unsigned nd=0; nd<nnode; nd++)
        {
          unsigned out_eqn = Lookup_unpinned_output
            .node_to_bemeq(m_pt->boundary_node_pt(b,nd));
          (*bem_output_values[i])[nd] = full_vector[out_eqn];
        }
    }

}

void PinnedBoundaryElementHandler::get_bem_values_and_copy_into_values() const
  {
    // Get as one big vector
    DoubleVector full_vector;
    get_bem_values(full_vector);

    //??ds use new lookups


    for(unsigned i=0, ni=Bem_boundaries.size(); i < ni; i++)
      {
        // Get info on this boundary
        unsigned b = Bem_boundaries[i].first;
        const Mesh* m_pt = Bem_boundaries[i].second;
        unsigned nnode = m_pt->nboundary_node(b);

        // Set the entry corresponding to output_index in this node
        // to the corresponding value from the doublevector.
        for(unsigned nd=0; nd<nnode; nd++)
          {
            unsigned out_eqn = Lookup_unpinned_output
              .node_to_bemeq(m_pt->boundary_node_pt(b,nd));
            m_pt->boundary_node_pt(b, nd)->set_value(output_index(),
                                                     full_vector[out_eqn]);
          }
      }
  }

} // End of oomph namespace
