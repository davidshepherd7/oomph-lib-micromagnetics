#ifndef OOMPH_MICROMAGNETICS_FLUX_ELEMENT_H
#define OOMPH_MICROMAGNETICS_FLUX_ELEMENT_H

/*
  A face element to impose the flux boundary condition on the potential for
  micromagnetics. It is not possible to use the already existing Poisson flux
  elements becuase they assume that:

  1) The "parent" element will be called PoissonElement.

  2) The Jacobian contribution due to the flux elements is always zero. This is
  not true for our case since the flux is fixed as m.n (where m is the
  magnetisation). So differentiation wrt m gives a non-zero result.

*/

#include "generic.h"

using namespace oomph;

namespace oomph
{

  //======================================================================
  ///
  //======================================================================
  template <class ELEMENT>
  class MicromagFluxElement : public virtual FaceGeometry<ELEMENT>,
			      public virtual FaceElement
  {

  public:

    /// \short Constructor, takes the pointer to the "bulk" element and the
    /// index of the face to which the element is attached.
    MicromagFluxElement(ELEMENT* const &bulk_el_pt,
			const int& face_index);

    ///\short  Broken empty constructor
    MicromagFluxElement()
    {
      throw OomphLibError
	("Don't call empty constructor for MicromagFluxElement",
	 "MicromagFluxElement::MicromagFluxElement()",
	 OOMPH_EXCEPTION_LOCATION);
    }

    /// Broken copy constructor
    MicromagFluxElement(const MicromagFluxElement& dummy)
    {
      BrokenCopy::broken_copy("MicromagFluxElement");
    }

    /// Broken assignment operator
    void operator=(const MicromagFluxElement&)
    {
      BrokenCopy::broken_assign("MicromagFluxElement");
    }

    // /// \short Specify the value of nodal zeta from the face geometry
    // /// The "global" intrinsic coordinate of the element when
    // /// viewed as part of a geometric object should be given by
    // /// the FaceElement representation, by default (needed to break
    // /// indeterminacy if bulk element is SolidElement)
    // //??ds don't think I need this....
    // double zeta_nodal(const unsigned &n, const unsigned &k,
    // 		      const unsigned &i) const
    // {return FaceElement::zeta_nodal(n,k,i);}

    /// Add the element's contribution to its residual vector
    void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
      // Call the generic residuals function with flag set to 0 using a dummy
      // matrix argument
      fill_in_generic_residual_contribution_fluxes
    	(residuals,GeneralisedElement::Dummy_matrix,0);
    }

    /// \short Add the element's contribution to its residual vector and its
    /// Jacobian matrix
    void fill_in_contribution_to_jacobian(Vector<double> &residuals,
    					  DenseMatrix<double> &jacobian)
    {
      //Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution_fluxes(residuals,jacobian,1);
    }

    void fill_in_bulk_contribution_to_face_jacobian
    (DenseMatrix<double>& jacobian) const;

    /// \short Output function -- forward to broken version in FiniteElement
    /// until somebody decides what exactly they want to plot here...
    void output(std::ostream &outfile, const unsigned &n_plot=5)
    {FiniteElement::output(outfile,n_plot);}

    /// \short C-style output function -- forward to broken version in
    /// FiniteElement until somebody decides what exactly they want to plot
    /// here...
    void output(FILE* file_pt, const unsigned &n_plot=5)
    {FiniteElement::output(file_pt,n_plot);}

    /// Access function for the pointer the to bulk element
    ELEMENT* bulk_element_pt() const {return Bulk_element_pt;}

  protected:

    /// \short Function to compute the shape and test functions and to return
    /// the Jacobian of mapping between local and global (Eulerian)
    /// coordinates
    double shape_test(const Vector<double> &s, Shape &psi, Shape &test)
      const
    {
      //Get the shape functions
      shape(s,psi);

      //Set the test functions to be the same as the shape functions
      test = psi;

      //Return the value of the jacobian
      return J_eulerian(s);
    }

    double dshape_dtest(const Vector<double>& s, Shape& psi, DShape& dpsidx,
			Shape& test, DShape& dtestdx)
      const
    {
      // Get s in bulk
      Vector<double> s_bulk(3,0.0);
      get_local_coordinate_in_bulk(s,s_bulk);

      // get dshape dtest from bulk
      unsigned bulk_n_node = bulk_element_pt()->nnode();
      Shape dummy1(bulk_n_node), dummy2(bulk_n_node);
      DShape dpsidx_bulk(bulk_n_node,Dim), dtestdx_bulk(bulk_n_node,Dim);
      bulk_element_pt()->dshape_dtest(s_bulk,dummy1,dpsidx_bulk,dummy2,dtestdx_bulk);

      // Find shape derivatives which are actually in the face element
      for(unsigned l=0; l<this->nnode(); l++)
	{
	  unsigned l_bulk = bulk_node_number(l);
	  for(unsigned j=0; j<Dim; j++)
	    {
	      dpsidx(l,j) = dpsidx_bulk(l_bulk,j);
	      dtestdx(l,j) = dtestdx_bulk(l_bulk,j);
	    } //??ds should check this...
	}

      // Get shape, test and J
      return shape_test(s,psi,test);
    }

    /// \short Add the element's contribution to its residual vector. When
    /// flag=1 (or 0): do (or don't) compute the contribution to the Jacobian as
    /// well.
    void fill_in_generic_residual_contribution_fluxes
    (Vector<double> &residuals, DenseMatrix<double> &jacobian,
     const unsigned& flag);

  private:

    ///The spatial dimension of the problem
    unsigned Dim;

    /// Pointer to the attached bulk element
    ELEMENT* Bulk_element_pt;
  };

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////



  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, the
  /// index of the fixed local coordinate and its value represented
  /// by an integer (+/- 1), indicating that the face is located
  /// at the max. or min. value of the "fixed" local coordinate
  /// in the bulk element.
  //===========================================================================
  template<class ELEMENT>
  MicromagFluxElement<ELEMENT>::
  MicromagFluxElement(ELEMENT* const &bulk_el_pt,
		      const int &face_index) :
    FaceGeometry<ELEMENT>(),
    FaceElement(),
    Bulk_element_pt(bulk_el_pt)
  {
    //??ds left in from PoissonFluxElement
// #ifdef PARANOID
//     {
//       //Check that the element is not a refineable 3d element
//       ELEMENT* elem_pt = new ELEMENT;
//       //If it's three-d
//       if(elem_pt->dim()==3)
// 	{
// 	  //Is it refineable
// 	  if(dynamic_cast<RefineableElement*>(elem_pt))
// 	    {
// 	      //Issue a warning
// 	      OomphLibWarning
// 		("This flux element will not work correctly if nodes are hanging\n",
// 		 "MicromagFluxElement::Constructor",
// 		 OOMPH_EXCEPTION_LOCATION);
// 	    }
// 	}
//     }
// #endif

    // Let the bulk element build the FaceElement, i.e. setup the pointers
    // to its nodes (by referring to the appropriate nodes in the bulk
    // element), etc.
    bulk_element_pt()->build_face_element(face_index,this);

    // Extract the dimension of the problem from the dimension of
    // the first node
    Dim = this->node_pt(0)->ndim();

    // // Fill in the index values
    // Phi_1_index_micromag = bulk_element_pt()->phi_1_index_micromag();
    // M_index_micromag.push_back(bulk_element_pt()->m_index_micromag(0));
    // M_index_micromag.push_back(bulk_element_pt()->m_index_micromag(1));
    // M_index_micromag.push_back(bulk_element_pt()->m_index_micromag(2));

    // Exchange_index_micromag.push_back(bulk_element_pt()->exchange_index_micromag(0));
    // Exchange_index_micromag.push_back(bulk_element_pt()->exchange_index_micromag(1));
    // Exchange_index_micromag.push_back(bulk_element_pt()->exchange_index_micromag(2));

  }


  //===========================================================================
  /// Compute the element's residual vector and the Jacobian matrix.
  //===========================================================================
  template<class ELEMENT>
  void MicromagFluxElement<ELEMENT>::
  fill_in_generic_residual_contribution_fluxes
  (Vector<double> &residuals,
   DenseMatrix<double> &jacobian,
   const unsigned& flag)
  {
    const unsigned n_node = nnode();
    const unsigned n_intpt = integral_pt()->nweight();

    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node,Dim), dtestdx(n_node,Dim);

    const unsigned phi_1_index = bulk_element_pt()->phi_1_index_micromag();
    Vector<unsigned> m_index(3);
    for(unsigned i=0; i<3; i++)
      m_index[i] = bulk_element_pt()->m_index_micromag(i);

    // Get coefficients
    const double exch_c = bulk_element_pt()->get_exchange_coeff();
    const double llg_precess_c = bulk_element_pt()->get_llg_precession_coeff();

    //Loop over the integration points
    //--------------------------------
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
      {
	// Local coords
	Vector<double> s(Dim-1);
	for(unsigned i=0;i<(Dim-1);i++)
	  s[i] = integral_pt()->knot(ipt,i);

	const double J = shape_test(s,psi,test);
	const double W = integral_pt()->weight(ipt) * J;

	// Get the shape/test functions and derrivatives
	dshape_dtest(s,psi,dpsidx,test,dtestdx);

	Vector<double> itp_x(Dim,0.0), itp_m(3,0.0);
	for(unsigned l=0; l<n_node; l++)
	  {
	    for(unsigned j=0; j<Dim; j++)
	      itp_x[j] += nodal_position(l,j)*psi(l);
	    for(unsigned j=0; j<3; j++)
	      itp_m[j] += nodal_value(l,m_index[j])*psi(l);
	  }

	// We have to ask the bulk element to calculate dmdx for us because the
	// derivatives of shape functions for nodes in the bulk but not in the
	// face could be non-zero on the boundary.
	Vector<double> s_bulk(3,0.0);
	get_local_coordinate_in_bulk(s,s_bulk);
	DenseDoubleMatrix itp_dmdx(3,3,0.0);
	bulk_element_pt()->interpolated_dmdx_micromag(s_bulk,itp_dmdx);

	Vector<double> low_dim_normal(Dim,0.0), normal(3,0.0);
	outer_unit_normal(s,low_dim_normal);
	for(unsigned j=0; j<Dim; j++) normal[j] = low_dim_normal[j];

	// Some pre-calculations:
	Vector<double> ndotgradmi(3,0.0);
	for(unsigned i=0; i<3; i++)
	  for(unsigned j=0; j<Dim; j++)
	    ndotgradmi[i] += normal[j] * itp_dmdx(i,j);
	Vector<double> mxndotgradmi(3,0.0);
	VectorOps::cross(itp_m,ndotgradmi,mxndotgradmi);
	double mdotn = VectorOps::dot(normal,itp_m);


	// Loop over the test functions doing residual and Jacobian
	// contributions.
	for(unsigned l=0;l<n_node;l++)
	  {
	    // Get indicies for phi_1 and m equations
	    int phi_1_eqn = nodal_local_eqn (l,phi_1_index);
	    Vector<int> m_eqn(3,0);

	    // Phi contribution (if not a dirichlet b.c.)
	    if(phi_1_eqn >= 0)
	      residuals[phi_1_eqn] += mdotn*test(l)*W;

	    // Exchange contribution: (m x (dm/dn)) * test
	    for(unsigned i=0; i<3; i++)
	      {
	    	m_eqn[i] = nodal_local_eqn(l,m_index[i]);
	    	if(m_eqn[i] >=0)
	    	  residuals[m_eqn[i]] += llg_precess_c * exch_c
	    	    * mxndotgradmi[i] * test(l) * W;
	      }

	    //??ds still need to deal with surface anisotropy...


	    // Jacobian (phi_1 w.r.t m_i)
	    if(!flag) continue;
	    for(unsigned l2=0; l2<n_node; l2++)
	      {
		// Loop over which m we are differentiating w.r.t
		Vector<unsigned> m_unknown(3,0.0);
		for(unsigned j=0; j<3; j++)
		  {
		    m_unknown[j] = nodal_local_eqn(l2,m_index[j]);
		    if(!(m_unknown[j] >= 0)) continue;

		    // phi_1 residual w.r.t m_j
		    if(phi_1_eqn >= 0)
		      {
			jacobian(phi_1_eqn,m_unknown[j])
			  += psi(l2) * test(l) * normal[j] * W;
		      }

		    // m residuals w.r.t m_j due to exchange boundary effects
		    // ============================================================

		    // Jacobian contribution is (using product rule):
		    // (ihat x dm/dn) * psi * ...
		    // + (m x ihat) * dpsidn * ...
		    Vector<double> jhat(3,0.0); jhat[j] = 1.0;
		    Vector<double> jxdmdn(3,0.0);
		    VectorOps::cross(jhat,ndotgradmi,jxdmdn);

		    for(unsigned i=0; i<3; i++)
		      {
		    	if(!(m_eqn[i] >= 0)) continue;

		    	jacobian(m_eqn[i],m_unknown[j]) +=
		    	  exch_c * llg_precess_c * W * test(l) *
		    	  psi(l2) * jxdmdn[i];

		    	// Unfortuantely the second term of this equation is
		    	// non-zero for m_unknown in bulk nodes in the same bulk
		    	// element. This means we must also include a Jacobian
		    	// contribution due to the surface integral in the bulk
		    	// element. This is implemented in the function ??ds
		    	// called by the bulk element.
		      }
		  }
	      }
	  }
      }
  }

  template<class ELEMENT>
  void MicromagFluxElement<ELEMENT>::
  fill_in_bulk_contribution_to_face_jacobian(DenseMatrix<double>& jacobian) const
  {
    // Get indices
    Vector<unsigned> m_index(3);
    for(unsigned i=0; i<3; i++)
      m_index[i] = bulk_element_pt()->m_index_micromag(i);

    const unsigned face_n_node = nnode();
    const unsigned n_node = bulk_element_pt()->nnode();
    const unsigned n_intpt = integral_pt()->nweight();

    // Get coefficients
    double exch_c = bulk_element_pt()->get_exchange_coeff();
    double llg_precess_c = bulk_element_pt()->get_llg_precession_coeff();

    for(unsigned ipt=0;ipt<n_intpt;ipt++)
      {
	Vector<double> s(Dim-1);
	for(unsigned i=0;i<(Dim-1);i++)
	  s[i] = integral_pt()->knot(ipt,i);

	// Get s in bulk
	Vector<double> s_bulk(3,0.0);
	get_local_coordinate_in_bulk(s,s_bulk);

	// Get shape, test and derivatives from bulk
	Shape psi_bulk(n_node), test_bulk(n_node);
	DShape dpsidx_bulk(n_node,Dim), dtestdx_bulk(n_node,Dim);
	const double J = bulk_element_pt()->
	  dshape_dtest(s_bulk, psi_bulk, dpsidx_bulk,
		       test_bulk, dtestdx_bulk);
	const double W = integral_pt()->weight(ipt) * J;

	// Get interpolated values
	Vector<double> itp_x(Dim,0.0), itp_m(3,0.0);
	for(unsigned l=0; l<face_n_node; l++)
	  {
	    for(unsigned j=0; j<Dim; j++)
	      itp_x[j] += nodal_position(l,j)*psi_bulk(l);
	    for(unsigned j=0; j<3; j++)
	      itp_m[j] += nodal_value(l,m_index[j])*psi_bulk(l);
	  }
	DenseDoubleMatrix itp_dmdx(3,3,0.0);
	bulk_element_pt()->interpolated_dmdx_micromag(s_bulk,itp_dmdx);

	// Get normal
	Vector<double> low_dim_normal(Dim,0.0), normal(3,0.0);
	outer_unit_normal(s,low_dim_normal);
	for(unsigned j=0; j<Dim; j++) normal[j] = low_dim_normal[j];

	// Now calculate Jacobian
	for(unsigned l=0;l<face_n_node;l++)
	  {
	    for(unsigned l2=0;l2<n_node;l2++)
	      {

		// Storage for equation numbers, note that we need to
		// differentiate by ALL the bulk nodes.
		Vector<int> m_eqn(3);
		Vector<int> m_unknown(3);

		// Loop over m value that we are differentiating by
		for(unsigned j=0; j<3; j++)
		  {
		    m_unknown[j] = bulk_element_pt()->nodal_local_eqn(l2,m_index[j]);
		    if(!(m_unknown[j] >= 0)) continue;

		    Vector<double> jhat(3,0.0), mxj(3,0.0); jhat[j] = 1.0;
		    VectorOps::cross(itp_m,jhat,mxj);

		    double ndotgradpsil2 = 0;
		    for(unsigned i=0; i<3; i++)
		      ndotgradpsil2 += normal[i] * dpsidx_bulk(l2,i);

		    // Loop over residual that we are differentiating
		    for(unsigned i=0; i<3; i++)
		      {
			m_eqn[i] = nodal_local_eqn(l,m_index[i]);
			if(!(m_eqn[i] >= 0)) continue;

			jacobian(m_eqn[i],m_unknown[j]) +=
			  W * test_bulk(l) * exch_c * llg_precess_c
			  * ndotgradpsil2 * mxj[i];
		      }
		  }
	      }
	  }
      }
  }

} // End of oomph namespace

#endif
