#ifndef OOMPH_MICROMAGNETICS_BOUNDARY_ELEMENT_H
#define OOMPH_MICROMAGNETICS_BOUNDARY_ELEMENT_H


#include "generic.h"
#include "../micromagnetics_element.h"
#include "../micromagnetics_element.cc"
#include "./variable_quadrature.h"
#include <functional>

using namespace oomph;
using namespace MathematicalConstants;

namespace oomph
{

  //===========================================================================
  ///
  //===========================================================================
  template<class ELEMENT>
  class MicromagFaceElement : public virtual FaceGeometry<ELEMENT>,
			      public virtual FaceElement
  {

  public:

    /// \short Constructor, takes the pointer to the bulk element and the
    /// index of the face to which the element is attached.
    MicromagFaceElement(FiniteElement* const &bulk_el_pt,
			const int& face_index);

    ///\short  Broken empty constructor
    MicromagFaceElement()
    {
      throw OomphLibError("Don't call empty constructor for MicromagFaceElement",
			  "MicromagFaceElement::MicromagFaceElement()",
			  OOMPH_EXCEPTION_LOCATION);
    }

    /// Broken copy constructor
    MicromagFaceElement(const MicromagFaceElement& dummy)
    {
      BrokenCopy::broken_copy("MicromagFaceElement");
    }

    /// Broken assignment operator
    void operator=(const MicromagFaceElement&)
    {
      BrokenCopy::broken_assign("MicromagFaceElement");
    }

    /// \short Specify the value of nodal zeta from the face geometry
    /// The "global" intrinsic coordinate of the element when
    /// viewed as part of a geometric object should be given by
    /// the FaceElement representation, by default (needed to break
    /// indeterminacy if bulk element is SolidElement)
    double zeta_nodal(const unsigned &n, const unsigned &k,
		      const unsigned &i) const
    {return FaceElement::zeta_nodal(n,k,i);}

    /// Add the element's contribution to its residual vector - do
    /// nothing (no residuals to add).
    inline void fill_in_contribution_to_residuals(Vector<double> &dummy) {}

    /// \short Add the element's contribution to its residual vector and its
    /// Jacobian matrix
    inline void fill_in_contribution_to_jacobian(Vector<double> &dummy,
						 DenseMatrix<double> &boundary_matrix)
    {
      fill_in_be_contribution_adaptive(boundary_matrix);
    }

    /// Output function -- forward to broken version in FiniteElement
    /// until somebody decides what exactly they want to plot here...
    void output(std::ostream &outfile) const
    {//FiniteElement::output(outfile);}
    }

    /// \short Output function -- forward to broken version in FiniteElement
    /// until somebody decides what exactly they want to plot here...
    void output(std::ostream &outfile, const unsigned &n_plot) const
    {//FiniteElement::output(outfile,n_plot);}
    }


    // /// C-style output function -- forward to broken version in FiniteElement
    // /// until somebody decides what exactly they want to plot here...
    // void output(FILE* file_pt) const
    // {FiniteElement::output(file_pt);}

    // /// \short C-style output function -- forward to broken version in
    // /// FiniteElement until somebody decides what exactly they want to plot
    // /// here...
    // void output(FILE* file_pt, const unsigned &n_plot) const
    // {FiniteElement::output(file_pt,n_plot);}

    /// \short Function giving the normal derivative of the Green's function.
    /// The parameters x and y are the positions of source and point of interest
    /// (interchangable), n is the normal unit vector out of the surface.
    double green_normal_derivative(const Vector<double>& x,
				   const Vector<double>& y,
				   const Vector<double>& n) const;

    /// Const access function for mesh pointer
    inline Mesh* boundary_mesh_pt() const {return Boundary_mesh_pt;}

    /// Set function for mesh pointer
    inline void set_boundary_mesh_pt(Mesh* boundary_mesh_pointer)
    {Boundary_mesh_pt = boundary_mesh_pointer;}

  protected:

    /// \short Function to compute the shape and test functions and to return
    /// the Jacobian of mapping between local and global (Eulerian)
    /// coordinates
    inline double shape_and_test(const Vector<double> &s,
				 Shape &psi, Shape &test) const
    {
      // Get the shape function and set test = shape
      shape(s,psi);
      for(unsigned i=0;i<nnode();i++) {test[i] = psi[i];}

      //Return the value of the jacobian
      return J_eulerian(s);
    }


    /// \short Function to compute the shape and test functions and to return
    /// the Jacobian of mapping between local and global (Eulerian)
    /// coordinates
    inline double shape_and_test_at_knot(const unsigned &ipt,
					 Shape &psi, Shape &test) const
    {
      // Get the shape function and set test = shape
      shape_at_knot(ipt,psi);
      for(unsigned i=0;i<nnode();i++) {test[i] = psi[i];}

      //Return the value of the jacobian
      return J_eulerian_at_knot(ipt);
    }

  private:

    /// Add the element's contribution to the boundary element matrix using adaptive quadrature.
    void fill_in_be_contribution_adaptive(DenseMatrix<double> &boundary_matrix) const;

    /// \short Pointer to the boundary mesh (needed to access nodes
    /// outside of this element for calculation of boundary matrix).
    Mesh* Boundary_mesh_pt;

    /// The index at which phi_1 is stored
    unsigned Phi_1_index_micromag;

    /// The index at which phi_2 is stored
    unsigned Phi_2_index_micromag;

    /// The number of values to be stored at each boundary element node
    inline unsigned required_nvalue(const unsigned &n) const
    {return 0;}

  };

  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, the
  /// index of the fixed local coordinate and its value represented
  /// by an integer (+/- 1), indicating that the face is located
  /// at the max. or min. value of the "fixed" local coordinate
  /// in the bulk element.
  //===========================================================================
  template<class ELEMENT>
  MicromagFaceElement<ELEMENT>::
  MicromagFaceElement(FiniteElement* const &bulk_el_pt, const int &face_index)
    : FaceGeometry<ELEMENT>(), FaceElement()
  {
    // Let the bulk element build the FaceElement, i.e. setup the pointers
    // to its nodes (by referring to the appropriate nodes in the bulk
    // element), etc.
    bulk_el_pt->build_face_element(face_index,this);
  }

  /// Get boundary element matrix contributions for this element using
  /// an adaptive scheme.
  template<class ELEMENT>
  void MicromagFaceElement<ELEMENT>::
  fill_in_be_contribution_adaptive(DenseMatrix<double> &boundary_matrix)
    const
  {
     // Find out dimension of element
    const unsigned el_dim = dim(), node_dim = nodal_dimension();

    //Find out how many nodes there are
    const unsigned n_element_node = nnode();

 // Cast pointer to the final integration scheme
#ifdef PARANOID
    // Dynamic casts are slow but type checked
    //??ds put try/catch in here and call oomphlib error if fails?
    VariableGaussLegendre* variable_int_pt =
      dynamic_cast<VariableGaussLegendre*>(integral_pt());
#else
    VariableGaussLegendre* variable_int_pt =
      static_cast<VariableGaussLegendre*>(integral_pt());
#endif

    // Start of adaptive integration scheme
    //====================================

    // Set parameters for adaptive integration
    unsigned max_order = variable_int_pt->max_order();
    unsigned highest_order_used = 0;
    double reltol = 1e-4, reldiff;

    // Set up storage for data at each knot for each order
    Vector< Vector< Vector<double> > > x_kn, normal;
    Vector< Vector<double> > jw;
    Vector< Vector<Shape> > psi, test;

    // Loop over all source nodes on the boundary
    unsigned n_boundary_node = boundary_mesh_pt()->nnode();
    for(unsigned i_sn=0; i_sn<n_boundary_node; i_sn++)
      {
	// Get coordinates of source node
	Vector<double> source_node_x(node_dim,0.0);
	boundary_mesh_pt()->node_pt(i_sn)->position(source_node_x);

	// Reset the quadrature scheme order (magic number: 0 = automatically
	// select min_order when adaptive_scheme_next_order is used)
	unsigned adapt_order = 0;

	Vector<double> temp_bm_prev(n_element_node,0.0), temp_bm(n_element_node,5000);
	do
	  {
	    // Copy the previous result into another vector (for comparison later)
	    // and initialise current result vector.
	    //??ds do the initilisation better?
	    temp_bm_prev = temp_bm;
	    for(unsigned l=0; l<n_element_node; l++)
	      {
		temp_bm[l] = 0.0;
	      }

	    // Move to next adaptive scheme order
	    adapt_order = variable_int_pt->adaptive_scheme_next_order(adapt_order);

	    ////////////////////////////////////////
	    // Precalculations (if needed)

	    // Check if this order of integration has been used yet, if not calculate things
	    if( !(adapt_order <= highest_order_used))
	      {
		// Get the number of knots at this order
		unsigned n_knot = variable_int_pt->nweight(adapt_order);

		// Resize the vectors to store the data for the new order
		x_kn.resize(adapt_order+1);
		normal.resize(adapt_order+1);
		psi.resize(adapt_order+1);
		test.resize(adapt_order+1);
		jw.resize(adapt_order+1);

		for(unsigned kn=0; kn<n_knot; kn++)
		  {
		    // Get the local coordinate
		    Vector<double> s(el_dim,0.0);
		    for(unsigned j=0; j<el_dim; j++)
		      s[j] = variable_int_pt->knot(kn,j,adapt_order);

		    // Get the unit normal
		    Vector<double> local_normal(node_dim,0.0);
		    outer_unit_normal(s,local_normal);
		    normal[adapt_order].push_back(local_normal);

		    // Get the Jacobian*weight
		    jw[adapt_order].push_back(J_eulerian(s)
					      * variable_int_pt->weight(kn,adapt_order));

		    // Get shape and test(=shape) functions
		    Shape psi_local(n_element_node);
		    shape(s,psi_local);
		    psi[adapt_order].push_back(psi_local);
		    test[adapt_order].push_back(psi_local);

		    // Get the global coordinate
		    Vector<double> x_local(node_dim,0.0);
		    for(unsigned l=0; l<n_element_node; l++)
		      {
			for(unsigned j=0; j<node_dim; j++)
			  x_local[j] += nodal_position(l,j)*psi[adapt_order][kn][l];
		      }
		    x_kn[adapt_order].push_back(x_local);

		  }

		// Set the new highest order
		highest_order_used = adapt_order;
	      }

	    ////////////////////////////////////
	    // The calculation itself

		double n_knot = variable_int_pt->nweight(adapt_order);
		for(unsigned kn=0;kn<n_knot;kn++)
		  {
		    double dgreendn =
		      green_normal_derivative(x_kn[adapt_order][kn],source_node_x,
					      normal[adapt_order][kn]);

		    // Loop over test functions, i.e. local/target nodes, adding contributions
		    for(unsigned l=0; l<n_element_node; l++)
		      temp_bm[l] -= dgreendn * test[adapt_order][kn][l]*jw[adapt_order][kn];
		  }

	    ////////////////////////////////////////////
	    // Evaluate the results

	    // Get the differences between the the results from the last two
	    // quadrature schemes tried.
	    Vector<double> diff_bm(n_element_node);
	    for(unsigned l=0; l<n_element_node; l++)
	      diff_bm[l] = fabs(temp_bm_prev[l] - temp_bm[l]);

	    // Get the worst case error (the maximum) and calculate the relative
	    // error (divide error by the appropriate integral value).
	    Vector<double>::iterator worst_error_it
	      = max_element(diff_bm.begin(),diff_bm.end());
	    reldiff = *worst_error_it
	      / temp_bm[worst_error_it - diff_bm.begin()];
	    //??ds need to put in something special to handle zero integrals

	  }
	while(((reldiff>reltol) && (adapt_order<max_order))
	      || (adapt_order==variable_int_pt->adaptive_scheme_next_order(0)));
	// Repeat unless the difference is small (i.e. quadrature has
	// converged), terminate if max_order has been reached, continue anyway
	// if we are still on the first order.

	// If we hit the order limit without being accurate enough give an error
	if (adapt_order >= max_order)
	  {
	    throw OomphLibError("Quadrature order not high enough.",
				"MicromagFaceElement::fill_in_be_contribution_adaptive",
				OOMPH_EXCEPTION_LOCATION);
	  }

	// When done add the values in the temp vector to the real boundary matrix
	for(unsigned l=0; l<n_element_node; l++)
	  boundary_matrix(l,i_sn) += temp_bm[l];

      } // End of loop over source nodes

  } // End of function

  //======================================================================
  ///
  //======================================================================
  template<class ELEMENT>
  double MicromagFaceElement<ELEMENT>::
  green_normal_derivative(const Vector<double>& x,
			  const Vector<double>& y,
			  const Vector<double>& n) const
  {
    // Get dimensions
    const unsigned node_dim = nodal_dimension();

#ifdef PARANOID
    //??ds check that n is a unit vector
#endif

    // Calculate the distance r from x to y.
    double r(0.0), r_sq(0.0);
    for(unsigned i=0; i<node_dim; i++)
      r_sq += pow((x[i] - y[i]),2);

    r = sqrt(r_sq);

    //??ds if x is a potentially singular point CHEAT
    //??ds could go horribly wrong, probably fine so long as not actually singular
    //??ds only works in 2d anyway
    if(r < 1e-14)
      {
	// This is approximately true because the r and n unit vectors
	// become perpendicular as x approaches y. Hence n.r = 0 and
	// the whole function is zero.
	return 0;
      }

    // Calculate the unit vector in direction r = x - y
    Vector<double> r_unit(node_dim,0.0);
    for(unsigned i=0; i<node_dim; i++)
      r_unit[i] = (x[i] - y[i])/r;


#ifdef PARANOID
    //??ds check that r is not ~machine error (division by zero)
#endif

    // Calculate dot product of n with the unit vector in direction r.
    double ndotr(0.0);
    for(unsigned i=0; i<node_dim; i++)
      ndotr += r_unit[i]*n[i];

    // std::cout << "ndotr = " << ndotr << std::endl;
    // std::cout << "r = " << r << std::endl;
    // std::cout << "greens/ndotr = " << -1/Pi * pow((2*r),-2) * 2 << std::endl;
    // std::cout << "greens = " <<-1/Pi * ndotr * pow((2*r),-2) * 2 << std::endl;

    // dgreendn = -n dot r * 1/pi * (1/2)^(node_dim-1) * (1/r)^node_dim
    // See write up for details of calculation.
    return -1/Pi * ndotr * pow((2*r),-2) * 2;
  }

}

#endif
