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

    /// \short Specify the value of nodal zeta from the face geometry
    /// The "global" intrinsic coordinate of the element when
    /// viewed as part of a geometric object should be given by
    /// the FaceElement representation, by default (needed to break
    /// indeterminacy if bulk element is SolidElement)
    //??ds don't think I need this....
    double zeta_nodal(const unsigned &n, const unsigned &k,
		      const unsigned &i) const
    {return FaceElement::zeta_nodal(n,k,i);}

    /// Add the element's contribution to its residual vector
    inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
      // Call the generic residuals function with flag set to 0 using a dummy
      // matrix argument
      fill_in_generic_residual_contribution_poisson_flux
	(residuals,GeneralisedElement::Dummy_matrix,0);
    }


    /// \short Add the element's contribution to its residual vector and its
    /// Jacobian matrix
    inline void fill_in_contribution_to_jacobian(Vector<double> &residuals,
						 DenseMatrix<double> &jacobian)
    {
      //Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution_poisson_flux(residuals,jacobian,1);
    }

    /// Output function -- forward to broken version in FiniteElement
    /// until somebody decides what exactly they want to plot here...
    void output(std::ostream &outfile) {FiniteElement::output(outfile);}

    /// \short Output function -- forward to broken version in FiniteElement
    /// until somebody decides what exactly they want to plot here...
    void output(std::ostream &outfile, const unsigned &n_plot)
    {FiniteElement::output(outfile,n_plot);}


    /// C-style output function -- forward to broken version in FiniteElement
    /// until somebody decides what exactly they want to plot here...
    void output(FILE* file_pt) {FiniteElement::output(file_pt);}

    /// \short C-style output function -- forward to broken version in
    /// FiniteElement until somebody decides what exactly they want to plot
    /// here...
    void output(FILE* file_pt, const unsigned &n_plot)
    {FiniteElement::output(file_pt,n_plot);}

    /// Access function for the pointer the to bulk element
    ELEMENT* bulk_element_pt() const {return Bulk_element_pt;}

  protected:

    /// \short Function to compute the shape and test functions and to return
    /// the Jacobian of mapping between local and global (Eulerian)
    /// coordinates
    inline double shape_and_test(const Vector<double> &s, Shape &psi, Shape &test)
      const
    {
      //Find number of nodes
      unsigned n_node = nnode();

      //Get the shape functions
      shape(s,psi);

      //Set the test functions to be the same as the shape functions
      for(unsigned i=0;i<n_node;i++) {test[i] = psi[i];}

      //Return the value of the jacobian
      return J_eulerian(s);
    }


    /// \short Function to compute the shape and test functions and to return
    /// the Jacobian of mapping between local and global (Eulerian)
    /// coordinates
    inline double shape_and_test_at_knot(const unsigned &ipt,
					 Shape &psi, Shape &test)
      const
    {
      //Find number of nodes
      unsigned n_node = nnode();

      //Get the shape functions
      shape_at_knot(ipt,psi);

      //Set the test functions to be the same as the shape functions
      for(unsigned i=0;i<n_node;i++) {test[i] = psi[i];}

      //Return the value of the jacobian
      return J_eulerian_at_knot(ipt);
    }


  private:


    /// \short Add the element's contribution to its residual vector.
    /// flag=1(or 0): do (or don't) compute the contribution to the
    /// Jacobian as well.
    void fill_in_generic_residual_contribution_poisson_flux
    (Vector<double> &residuals,
     DenseMatrix<double> &jacobian,
     const unsigned& flag);

    ///The spatial dimension of the problem
    unsigned Dim;

    ///The index at which phi_1 is stored
    unsigned Phi_1_index_micromag;

    /// The indices at which m is stored
    Vector<unsigned> M_index_micromag;

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
    M_index_micromag(3,0),
    Bulk_element_pt(bulk_el_pt)
  {
    //??ds left in from PoissonFluxElement
#ifdef PARANOID
    {
      //Check that the element is not a refineable 3d element
      ELEMENT* elem_pt = new ELEMENT;
      //If it's three-d
      if(elem_pt->dim()==3)
	{
	  //Is it refineable
	  if(dynamic_cast<RefineableElement*>(elem_pt))
	    {
	      //Issue a warning
	      OomphLibWarning
		("This flux element will not work correctly if nodes are hanging\n",
			      "MicromagFluxElement::Constructor",
			      OOMPH_EXCEPTION_LOCATION);
	    }
	}
    }
#endif

    // Let the bulk element build the FaceElement, i.e. setup the pointers
    // to its nodes (by referring to the appropriate nodes in the bulk
    // element), etc.
    bulk_element_pt()->build_face_element(face_index,this);

    // Extract the dimension of the problem from the dimension of
    // the first node
    Dim = this->node_pt(0)->ndim();

    // Fill in the index values
    Phi_1_index_micromag = bulk_element_pt()->phi_1_index_micromag();
    M_index_micromag[0] = bulk_element_pt()->m_index_micromag(0);
    M_index_micromag[1] = bulk_element_pt()->m_index_micromag(1);
    M_index_micromag[2] = bulk_element_pt()->m_index_micromag(2);
  }


  //===========================================================================
  /// Compute the element's residual vector and the (zero) Jacobian matrix.
  //===========================================================================
  template<class ELEMENT>
  void MicromagFluxElement<ELEMENT>::
  fill_in_generic_residual_contribution_poisson_flux
  (Vector<double> &residuals,
   DenseMatrix<double> &jacobian,
   const unsigned& flag)
  {
    const unsigned n_node = nnode();
    Shape psif(n_node), testf(n_node);
    const unsigned n_intpt = integral_pt()->nweight();

    //Integer to hold the local equation number
    int local_eqn=0;

    //Loop over the integration points
    //--------------------------------
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
      {
	// Get value of s
	Vector<double> s(Dim-1);
	for(unsigned i=0;i<(Dim-1);i++)
	  {s[i] = integral_pt()->knot(ipt,i);}

	// Get value of s in bulk element
	Vector<double> s_bulk(Dim,0.0);
	get_local_coordinate_in_bulk(s,s_bulk);

	//Get the integral weight
	double w = integral_pt()->weight(ipt);

	//Find the shape and test functions and return the Jacobian
	//of the mapping
	double J = shape_and_test(s,psif,testf);

	//Premultiply the weights and the Jacobian
	double W = w*J;

	// Get outward unit normal at this integration point
	Vector<double> normal(Dim,0.0);
	outer_unit_normal(s,normal);

	// Get m at this integration point from the bulk element
	Vector<double> interpolated_m(3,0.0);
	bulk_element_pt()->interpolated_m_micromag(s_bulk,interpolated_m);

	// Calculate prescribed flux =  m.n
	double flux = 0;
	for(unsigned i=0; i<Dim; i++)
	  {
	    flux += normal[i]*interpolated_m[i];
	  }

	//Loop over the test functions
	for(unsigned l=0;l<n_node;l++)
	  {
	    local_eqn = nodal_local_eqn(l,Phi_1_index_micromag);

	    if(local_eqn >= 0) 	//If it's not a Dirichlet boundary condition
	      {
		//Add the prescribed flux terms
		residuals[local_eqn] -= flux*testf[l]*W;

		//??ds fill in Jacobian when possible
	      }
	  }
      }
  }





} // End of oomph namespace

#endif
