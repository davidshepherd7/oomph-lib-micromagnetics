/*
  description of file goes here
*/


using namespace oomph;
using namespace MathematicalConstants;


namespace oomph
{

  /// Assign the static boundary mesh pointer for the bem elements. We have to do
  /// this seperately to the class because it is a static variable.
  template<class ELEMENT,unsigned DIM>
  Mesh* MicromagFaceElement<ELEMENT,DIM>::Boundary_mesh_pt=0;


  //======================================================================
  /// Compute element residual Vector and/or element Jacobian matrix
  ///
  /// flag=1: compute both
  /// flag=0: compute only residual Vector
  //======================================================================
  template<class ELEMENT, unsigned DIM>
  void MicromagFaceElement<ELEMENT,DIM>::fill_in_generic_residual_contribution_micromag_boundary
  (Vector<double> &residuals, DenseMatrix<double> &jacobian, const unsigned& flag) const
  {
    // The Jacobians and residuals added by this element are somewhat different
    // to the typical oomph-lib case:

    // There is no integration to be done here, the residual is simply this
    // difference. This is because the boundary element method part is done
    // using the co-location method rather than a Galerkin method.

    // So the residual at this node is the difference between the value phi
    // should have according to the BEM and the actual value.

    // We use the opposite sign convention here to elsewhere (here: predicted
    // value minus actual value) because it simplifies the addition of the
    // boundary matrix to the Jacobian. Should be ok.


    // ??ds if something does go wrong these calculation are something to double
    // check...

    for(unsigned l=0; l<nnode(); l++)
      {
    	int phi_eqn = nodal_local_eqn(l,phi_index_micromag());
	if(phi_eqn >= 0)
	  {
	    residuals[phi_eqn] = local_bem_phi_value(l)
	      - nodal_value(l,phi_index_micromag());
	  }
      }

    if(flag)
      {
	for(unsigned  l=0; l<nnode(); l++)
	  {
	    int phi_eqn = nodal_local_eqn(l,phi_index_micromag());
	    if(phi_eqn >= 0)
	      {
		jacobian(phi_eqn,phi_eqn) = -1;
	      }
	  }
      }

  }


} // End of oomph namespace
