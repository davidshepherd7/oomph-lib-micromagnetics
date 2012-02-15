#ifndef OOMPH_MICROMAGNETICS_BOUNDARY_ELEMENT_H
#define OOMPH_MICROMAGNETICS_BOUNDARY_ELEMENT_H


#include "generic.h"
#include "../micromagnetics_element.h"
#include "../micromagnetics_element.cc"
#include "./variable_quadrature.h"

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
      // fill_in_be_contribution(boundary_matrix);
      // fill_in_be_contribution_adaptive(boundary_matrix);
      fill_in_be_contribution_quadpack(boundary_matrix);
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
    inline Mesh* mesh_pt() const {return Mesh_pt;}

    /// Set function for mesh pointer
    inline void set_mesh_pt(Mesh* mesh_pointer) {Mesh_pt = mesh_pointer;}

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

    /// Add the element's contribution to the boundary element matrix.
    void fill_in_be_contribution(DenseMatrix<double> &boundary_matrix) const;

    /// Add the element's contribution to the boundary element matrix using adaptive quadrature.
    void fill_in_be_contribution_adaptive(DenseMatrix<double> &boundary_matrix) const;

    /// Add the element's contribution to the boundary element matrix using QUADPACK algorithms.
    void fill_in_be_contribution_quadpack(DenseMatrix<double> &boundary_matrix) const;

    /// \short Pointer to the boundary mesh (needed to access nodes
    /// outside of this element for calculation of boundary matrix).
    Mesh* Mesh_pt;

    /// The dimension of the element surface/volume (i.e. one less
    /// than the dimension of the nodes.
    unsigned Node_dim;

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

    // Extract the nodal dimension of the problem from the dimension of
    // the first (face) node.
    Node_dim = this->node_pt(0)->ndim();
  }

  /// Get boundary element matrix contributions for this element using
  /// QUADPACK routines (probably very inefficient).
  template<class ELEMENT>
  void MicromagFaceElement<ELEMENT>::
  fill_in_be_contribution_quadpack(DenseMatrix<double> &boundary_matrix)
    const
  {

  }

  /// Get boundary element matrix contributions for this element using
  /// an adaptive scheme.
  template<class ELEMENT>
  void MicromagFaceElement<ELEMENT>::
  fill_in_be_contribution_adaptive(DenseMatrix<double> &boundary_matrix)
    const
  {
    // Find out dimension of element
    const unsigned el_dim = Node_dim - 1;

    //Find out how many nodes there are
    const unsigned n_element_node = nnode();

    // Cast pointer to the Clenshaw-Curtis integration scheme
#ifdef PARANOID
    // Dynamic casts are slow but type checked
    //??ds put try/catch in here and call oomphlib error if fails?
    VariableFejerSecond* variable_int_pt =
      dynamic_cast<VariableFejerSecond*>(integral_pt());
#else
    VariableFejerSecond* variable_int_pt =
      static_cast<VariableFejerSecond*>(integral_pt());
#endif

    // Set parameters for adaptive integration
    unsigned max_order = 49;
    unsigned min_order = 2;
    double abstol = 1e-3;

    variable_int_pt->set_order(max_order);

    // Then pre-calculate everything at these points:

    // Find out how many knots there are
    const unsigned n_knot_max = variable_int_pt->nweight();

    // Set up vectors to store data at all knots
    Vector<double> J(n_knot_max,0.0);
    Vector<Vector<double> > s(n_knot_max,Vector<double>(el_dim,0.0)),
      interpolated_x(n_knot_max,Vector<double>(Node_dim,0.0)),
      normal(n_knot_max,Vector<double>(Node_dim,0.0));

    //??ds how can I make vectors of shape functions work?
    // Vector<Shape> psiv(n_knot_max, Shape(n_element_node)),
    //   testv(n_knot_max, Shape(n_element_node));

    // Loop over the knots pre-calculating values
    for(unsigned kn=0;kn<n_knot_max;kn++)
      {
	// Get the local coordinate at this knot
	for(unsigned j=0; j<el_dim; j++)
	  s[kn][j] = variable_int_pt->knot(kn,j);

	// Get the derivatives of the shape and test functions
	Shape psi(n_element_node), test(n_element_node);
	J[kn] = shape_and_test(s[kn],psi,test);

  	// Compute the normal vector
  	// ??ds not sure how this works - might be far too slow
  	outer_unit_normal(s[kn],normal[kn]);

  	// Get values of x (global coordinate)
  	for(unsigned j=0; j<el_dim; j++)
	  {
	    for(unsigned l=0; l<n_element_node; l++)
	      interpolated_x[kn][j] += nodal_position(l,j)*psi[l];
	  }
      }

    // Loop over all source nodes on the boundary
    unsigned n_boundary_node = mesh_pt()->nnode();
    for(unsigned i_sn=0; i_sn<n_boundary_node; i_sn++)
      {
  	// Get coordinates of source node
  	Vector<double> source_node_x(Node_dim,0.0);
  	mesh_pt()->node_pt(i_sn)->position(source_node_x);

	// Use an adaptive scheme: calculate at two lowest orders allowed
	// and compare. If they are close accept otherwise calculate the next
	// and compare... etc.
	unsigned order = min_order;
	double diff = 1e6;
	Vector<double> temp_bm_prev(n_element_node,0.0),
	  temp_bm(n_element_node,0.0);

	do
	  {
	    // Copy the previous result into another vector (for comparison later)
	    for(unsigned l=0; l<n_element_node; l++)
	      temp_bm_prev[l] = temp_bm[l];

	    // Get and set the next order to use
	    order = variable_int_pt->adaptive_scheme_next_order();
	    variable_int_pt->set_order(order);

	    // Get the shape and test functions at this knot
	    //??ds improve this by pre-calculating them
	    Shape psi(n_element_node), test(n_element_node);

	    // Loop over the knots used in this quadrature order
	    for(unsigned kn=0; kn<variable_int_pt->nweight(); kn++)
	      {
		// Get the knot corresponding to this one in high order scheme
		unsigned i =
		  variable_int_pt->find_corresponding_knot(kn,max_order);

		// Calculate dGreendn between source node and integration point
		double dgreendn = green_normal_derivative
		  (interpolated_x[i],source_node_x,normal[i]);

		// get the shape and test functions
		shape_and_test(s[i],psi,test);

		// Loop over test functions (i.e. local nodes)
		for(unsigned l=0; l<n_element_node; l++)
		  {
		    // Add contribution to integral (note dGreendn is
		    // negative in our definition). See write up for
		    // details.
		    temp_bm[l] = -dgreendn * test(l)
		      * J[i] * variable_int_pt->weight(kn);
		  }
	      }

	    // Get the worst case difference between the the results
	    // from the last two quadrature schemes tried.
	    Vector<double> diff_bm(n_element_node);
	    for(unsigned l=0; l<n_element_node; l++)
	      diff_bm[l] = fabs(temp_bm_prev[l] - temp_bm[l]);
	    diff = *max_element(diff_bm.begin(),diff_bm.end());
	    std::cout << diff << std::endl;
	  }
	while(((diff>abstol) && (order<max_order))
	      || (order==min_order));
	// Repeat unless the difference is small (i.e. quadrature has converged)
	// terminate if max_order has been reached
	// continue anyway if we are still on the first order.

	std::cout << std::endl;

	// If we hit the order limit without being accurate enough give an error
	if (variable_int_pt->order() >= max_order)
	  {
	    throw OomphLibError("Quadrature order not high enough.",
				"....",
				OOMPH_EXCEPTION_LOCATION);
	  }

	// When done add the values in the temp vector to the real boundary matrix
	for(unsigned l=0; l<n_element_node; l++)
	  boundary_matrix(l,i_sn) += temp_bm[l];
      }

  } // End of function

  //=======================================================================
  /// Not actually getting residuals - getting boundary element matrix but using
  /// machinery used to assemble jacobian normally.
  /// Compute the effect of phi_1 in this element on ALL nodes on the boundary.
  /// Note that this depends on all other elements but not on the values of phi_1.
  //??ds may have confused test and "psi" functions but it doesn't matter for now
  // since they are the same (Galerkin method)
  //=======================================================================
  template<class ELEMENT>
  void MicromagFaceElement<ELEMENT>::
  fill_in_be_contribution(DenseMatrix<double> &boundary_matrix)
    const
  {
    // Find out dimension of element
    const unsigned el_dim = Node_dim - 1;

    //Find out how many nodes there are
    const unsigned n_element_node = nnode();

    //Set up memory for the shape and test functions
    Shape psi(n_element_node), test(n_element_node);

    // // Get current time
    // double time = time_pt()->time();

    //??ds adaptive quadrature: given acceptable error choose integration
    // order/method and return integral_pt

    //Set the value of n_knot
    const unsigned n_knot = integral_pt()->nweight();

    //Loop over the integration points
    for(unsigned kn=0;kn<n_knot;kn++)
      {
    	// Get values of s (local coordinate)
    	Vector<double> s(el_dim,0.0);
    	for(unsigned j=0; j<el_dim; j++) {s[j] = integral_pt()->knot(kn,j);}

    	//Get the integral weight
    	double w = integral_pt()->weight(kn);

    	//Call the derivatives of the shape and test functions
    	double J = shape_and_test(s,psi,test);

    	//Premultiply the weights and the Jacobian
    	double W = w*J;

    	// Get values of x (global coordinate)
    	Vector<double> interpolated_x(Node_dim,0.0);
    	for(unsigned l=0; l<n_element_node; l++)
    	  {
    	    for(unsigned i=0; i<Node_dim; i++)
    	      interpolated_x[i] += nodal_position(l,i)*psi[l];
    	  }

    	// std::cout << "Gauss point at: (" << interpolated_x[0]
    	// 	  << ", " << interpolated_x[1] << ")" << std::endl;

    	// Compute the normal vector
    	//??ds not sure how this works - might not work for curved boundaries?
    	Vector<double> normal(Node_dim,0.0);
    	outer_unit_normal(s,normal);

    	// Loop over ALL nodes on in boundary mesh (except the current ones?) (source_node)
    	unsigned n_boundary_node = mesh_pt()->nnode();
    	for(unsigned i_sn=0; i_sn<n_boundary_node; i_sn++)
    	  {
    	    // Get coordinates of source node
    	    Vector<double> source_node_x(Node_dim,0.0);
    	    mesh_pt()->node_pt(i_sn)->position(source_node_x);

    	    // // Debugguging output:
    	    // std::cout << std::endl;
    	    // std::cout << ". Source node at: (" << source_node_x[0]
    	    // 	      << ", " << source_node_x[1] << ")"
    	    // 	      << std::endl;

    	    // Calculate dGreendn between source node and integration point
    	    double dgreendn = green_normal_derivative(interpolated_x,source_node_x,normal);

    	    // Loop over test functions, i.e. local/target nodes, adding contributions
    	    for(unsigned l=0; l<n_element_node; l++)
    	      {
    		// Add contribution to integral (note dGreendn is negative in our definition)
      		// See write up for details.
    		boundary_matrix(l,i_sn) -= dgreendn * test(l) * W;
    		// boundary_matrix(l,i_sn) += 1;
    	      }

    	  }
      }

    //??ds need to seperately add the contribution at each node from angles
  }

  //======================================================================
  ///
  //======================================================================
  template<class ELEMENT>
  double MicromagFaceElement<ELEMENT>::
  green_normal_derivative(const Vector<double>& x,
			  const Vector<double>& y,
			  const Vector<double>& n) const
  {

#ifdef PARANOID
    //??ds check that n is a unit vector
#endif

    // Calculate the distance r from x to y.
    double r(0.0), r_sq(0.0);
    for(unsigned i=0; i<Node_dim; i++)
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
    Vector<double> r_unit(Node_dim,0.0);
    for(unsigned i=0; i<Node_dim; i++)
      r_unit[i] = (x[i] - y[i])/r;


#ifdef PARANOID
    //??ds check that r is not ~machine error (division by zero)
#endif

    // Calculate dot product of n with the unit vector in direction r.
    double ndotr(0.0);
    for(unsigned i=0; i<Node_dim; i++)
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
