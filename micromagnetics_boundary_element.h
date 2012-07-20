#ifndef OOMPH_MICROMAGNETICS_BOUNDARY_ELEMENT_H
#define OOMPH_MICROMAGNETICS_BOUNDARY_ELEMENT_H


#include "generic.h"
#include "./micromagnetics_element.h"
#include "./micromagnetics_element.cc"
#include "./variable_order_quadrature.h"
#include <functional>

// magpar definitions
#include "./magpar_requirements.h"

// easier to use vector functions for now (probably not optimal...)
#include "./vector_helpers.h"

using namespace oomph;
using namespace MathematicalConstants;
using namespace VectorOps;

namespace oomph
{

  //===========================================================================
  ///
  //===========================================================================
  template<class ELEMENT,unsigned DIM>
  class MicromagFaceElement : public virtual FaceGeometry<ELEMENT>,
			      public virtual FaceElement
  {

  private:

    /// \short Pointer to the boundary mesh (needed to access nodes
    /// outside of this element for calculation of boundary matrix).
    static Mesh* Boundary_mesh_pt;

    ELEMENT* Bulk_element_pt;

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

    unsigned self_test()
    {
      // Just pass it (we can't let this element call the normal self-tests
      // because it messes with the integration scheme).
      return 0;
    }

    unsigned phi_index_micromag() const
    {return bulk_element_pt()->phi_index_micromag();}

    unsigned phi_1_index_micromag() const
    {return bulk_element_pt()->phi_1_index_micromag();}

    unsigned m_index_micromag(const unsigned& i) const
    {return bulk_element_pt()->m_index_micromag(i);}

    /// Pointer to higher-dimensional "bulk" element
    ELEMENT*& bulk_element_pt() {return Bulk_element_pt;}

    /// Pointer to higher-dimensional "bulk" element (const version)
    ELEMENT* bulk_element_pt() const {return Bulk_element_pt;}

    // double& local_bem_phi_value(const unsigned& l)
    // {return Local_bem_phi_value[l];}

    // double local_bem_phi_value(const unsigned& l) const
    // {return Local_bem_phi_value[l];}

    /// \short Specify the value of nodal zeta from the face geometry
    /// The "global" intrinsic coordinate of the element when
    /// viewed as part of a geometric object should be given by
    /// the FaceElement representation, by default (needed to break
    /// indeterminacy if bulk element is SolidElement)
    double zeta_nodal(const unsigned &n, const unsigned &k,
		      const unsigned &i) const
    {return FaceElement::zeta_nodal(n,k,i);}

    // /// Add the element's contribution to its residual vector and Jacobian - all
    // /// contributions due to BEM are done here.
    // void fill_in_generic_residual_contribution_micromag_boundary
    // (Vector<double> &residuals, DenseMatrix<double> &jacobian,
    //  const unsigned& flag) const;

    // /// Add the element's contribution to its residual vector (wrapper)
    // void fill_in_contribution_to_residuals(Vector<double> &residuals)
    // {
    //   //Call the generic residuals function with flag set to 0 using a dummy matrix argument
    //   fill_in_generic_residual_contribution_micromag_boundary
    // 	(residuals,GeneralisedElement::Dummy_matrix, 0);
    // }

    // /// \short Add the element's contribution to its residual vector and element
    // /// Jacobian matrix (wrapper)
    // void fill_in_contribution_to_jacobian(Vector<double> &residuals,
    // 					  DenseMatrix<double> &jacobian)
    // {
    //   // Call the generic routine with the flag set to 1
    //   fill_in_generic_residual_contribution_micromag_boundary
    // 	(residuals,jacobian,1);
    // }

    /// \short Add the element's contribution to its residual vector and its
    /// Jacobian matrix
    void fill_in_contribution_to_boundary_matrix(DenseMatrix<double> &boundary_matrix) const
    {
      fill_in_be_contribution_analytic(boundary_matrix);
      //fill_in_be_contribution_adaptive(boundary_matrix);
    }

    /// Function determining how to block the Jacobian.
    void get_dof_numbers_for_unknowns(std::list<std::pair<unsigned long,unsigned> >&
				      block_lookup_list)
    {
      unsigned phi_1_index = phi_1_index_micromag();
      unsigned phi_index = phi_index_micromag();
      unsigned bulk_dof_types = bulk_element_pt()->ndof_types();

      // We just need to move the boundary values of phi into their own blocks
      for(unsigned nd=0; nd<nnode(); nd++)
	{
	  int phi_1_local = bulk_element_pt()->nodal_local_eqn(nd,phi_1_index);
	  if(phi_1_local >=0)
	    {
	      std::pair<unsigned,unsigned> block_lookup;
	      block_lookup.first = bulk_element_pt()->eqn_number(phi_1_local);
	      block_lookup.second = bulk_dof_types + 1;
	      block_lookup_list.push_front(block_lookup);
	    }

	  int phi_local = bulk_element_pt()->nodal_local_eqn(nd,phi_index);
	  if(phi_local >=0)
	    {
	      std::pair<unsigned,unsigned> block_lookup;
	      block_lookup.first = bulk_element_pt()->eqn_number(phi_local);
	      block_lookup.second = bulk_dof_types + 1;
	      block_lookup_list.push_front(block_lookup);
	    }
	}
    }

    unsigned ndof_types()
    {
      return bulk_element_pt()->ndof_types() + additional_ndof_types();
    }

    unsigned additional_ndof_types()
    { return 2; }

    /// Output function - do nothing
    void output(std::ostream &outfile, const unsigned &n_plot=5) const {}

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
    static Mesh* boundary_mesh_pt() {return Boundary_mesh_pt;}

    /// Set function for mesh pointer
    static void set_boundary_mesh_pt(Mesh* boundary_mesh_pointer)
    {Boundary_mesh_pt = boundary_mesh_pointer;}

    /// Get the max difference between two vectors relative to that element of vector1
    double max_rel_error(const Vector<double> &vector1,
			 const Vector<double> &vector2) const;

    /// dump out important values for testing
    void dump_values(const Vector< Vector<double> > &x_kn,
		     const Vector<double> &source_node_x,
		     const Vector<Vector<double> > &normal) const;

  protected:

    /// \short Function to compute the shape and test functions and to return
    /// the Jacobian of mapping between local and global (Eulerian)
    /// coordinates
    inline double shape_and_test(const Vector<double> &s,
				 ShapeWithDeepCopy &psi, ShapeWithDeepCopy &test) const
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
					 ShapeWithDeepCopy &psi, ShapeWithDeepCopy &test) const
    {
      // Get the shape function and set test = shape
      shape_at_knot(ipt,psi);

      for(unsigned i=0;i<nnode();i++) {test[i] = psi[i];}

      //Return the value of the jacobian
      return J_eulerian_at_knot(ipt);
    }

    /// The number of values to be stored at each boundary element node
    inline unsigned required_nvalue(const unsigned &n) const
    {return 0;}

  private:

    /// Add the element's contribution to the boundary element matrix using adaptive quadrature.
    void fill_in_be_contribution_adaptive(DenseMatrix<double> &boundary_matrix) const;

    /// Add the element's contribution to the boundary element matrix using
    /// analytic calculations by calling analytic_integral_dgreendn with
    /// appropriate input for the element.
    void fill_in_be_contribution_analytic(DenseMatrix<double> &boundary_matrix) const;

    /// Calculate the contribution of a triangular region to the boundary
    /// element matrix using analytic calculations from Lindholm1984.
    void analytic_integral_dgreendn_triangle(const Vector<Vector<double> >& x_nds,
					     const Vector<unsigned>& l,
					     DenseMatrix<double>& boundary_matrix) const;

  };

  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, the
  /// index of the fixed local coordinate and its value represented
  /// by an integer (+/- 1), indicating that the face is located
  /// at the max. or min. value of the "fixed" local coordinate
  /// in the bulk element.
  //===========================================================================
  template<class ELEMENT,unsigned DIM>
  MicromagFaceElement<ELEMENT,DIM>::
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
  template<class ELEMENT,unsigned DIM>
  void MicromagFaceElement<ELEMENT,DIM>::
  fill_in_be_contribution_adaptive(DenseMatrix<double> &boundary_matrix)
    const
  {
    // Find out dimension of element
    const unsigned el_dim = DIM-1, node_dim = nodal_dimension();

    //Find out how many nodes there are
    const unsigned n_element_node = nnode();

    // Cast pointer to the final integration scheme
#ifdef PARANOID
    // Dynamic casts are slow but type checked
    //??ds put try/catch in here and call oomphlib error if fails?
    // don't think this error check is working...
    // need to generalise to pick the right scheme (T or Q)
    QVariableOrderGaussLegendre<el_dim>* v_int_pt =
      dynamic_cast<QVariableOrderGaussLegendre<el_dim>* >(integral_pt());
#else
    QVariableOrderGaussLegendre<el_dim>* v_int_pt =
      static_cast<QVariableOrderGaussLegendre<el_dim>*>(integral_pt());
#endif

    // Start of adaptive integration scheme
    //====================================

    // Setup error parameters for adaptive integration
    double reltol = 1e-10, reldiff;

    // Set up storage to keep track of orders used:
    //std::vector<unsigned> order_list = {2,4,8,16,32,64,128,256,};
    unsigned array_list[] = {2,4,8,16,32,64,128,256,};
    std::vector<unsigned> order_list(array_list,array_list+7);
    unsigned n_order = order_list.size();
    Vector<unsigned> order_is_calculated(n_order,0);

    // Set up storage for data at each knot for each order
    Vector< Vector< Vector<double> > > x_kn(n_order), normal(n_order);
    Vector< Vector<double> > jw(n_order);
    Vector< Vector<ShapeWithDeepCopy> > test(n_order), psi(n_order);

    // Loop over all source nodes on the boundary
    unsigned n_boundary_node = boundary_mesh_pt()->nnode();
    for(unsigned i_sn=0; i_sn<n_boundary_node; i_sn++)
      {
	// Get coordinates of source node
	Vector<double> source_node_x(node_dim,0.0);
	boundary_mesh_pt()->node_pt(i_sn)->position(source_node_x);

	// Reset the quadrature order
	unsigned i_order = 0;

	// Temporary storage for the values computed
	Vector<double> temp_bm_prev(n_element_node,0.0), temp_bm(n_element_node,5000);
	do
	  {
	    // Copy the previous result into another vector (for comparison later)
	    // and initialise current result vector.
	    temp_bm_prev = temp_bm;
	    temp_bm.assign(n_element_node,0.0);

	    ////////////////////////////////////////
	    // Precalculations (if needed)

	    // Check if this order of integration has been used yet, if not calculate things
	    if( order_is_calculated[i_order] != 1)
	      {
		// Get new order
		unsigned new_order = order_list[i_order];

		// Get the number of knots at this order
		unsigned n_knot = v_int_pt->nweight(new_order);

		// Calculate and push back values at each knot
		for(unsigned kn=0; kn<n_knot; kn++)
		  {
		    // Get the local coordinate
		    Vector<double> s(el_dim,0.0);
		    for(unsigned j=0; j<el_dim; j++)
		      s[j] = v_int_pt->knot(kn,j,new_order);

		    // Get the unit normal
		    Vector<double> local_normal(node_dim,0.0);
		    outer_unit_normal(s,local_normal);
		    normal[i_order].push_back(local_normal);

		    // Get the Jacobian*weight
		    jw[i_order].push_back(J_eulerian(s)
					  * v_int_pt->weight(kn,new_order));

		    // Get shape and test(=shape) functions
		    ShapeWithDeepCopy psi_local(n_element_node);
		    shape(s,psi_local);
		    psi[i_order].push_back(psi_local);
		    test[i_order].push_back(psi_local);

		    // Get the global coordinate
		    Vector<double> x_local(node_dim,0.0);
		    for(unsigned l=0; l<n_element_node; l++)
		      {
			for(unsigned j=0; j<node_dim; j++)
			  x_local[j] += nodal_position(l,j)*psi[i_order][kn][l];
		      }
		    x_kn[i_order].push_back(x_local);

		  }

		// Record that the new order has been calculated
		order_is_calculated[i_order] = 1;
	      }

	    ////////////////////////////////////////////////////////////
	    // The calculation itself

	    unsigned order = order_list[i_order];

	    double n_knot = v_int_pt->nweight(order);
	    for(unsigned kn=0;kn<n_knot;kn++)
	      {
		double dgreendn =
		  green_normal_derivative(x_kn[i_order][kn],source_node_x,
					  normal[i_order][kn]);

		// Loop over test functions, i.e. local/target nodes, adding contributions
		for(unsigned l=0; l<n_element_node; l++)
		  temp_bm[l] -= dgreendn *test[i_order][kn][l] *jw[i_order][kn];
	      }
	    ////////////////////////////////////////////////////////////

	    // Output values at knots if we are (probably) about to fail
	    if (i_order > 4)
	      {
		std::cout << i_order << std::endl;
		dump_values(x_kn[i_order],source_node_x,normal[i_order]);
	      }


	    // Compare with previous result
	    reldiff = max_rel_error(temp_bm,temp_bm_prev);

	    // if(i_order > 4)
	    //   {
	    // 	std::cout << "i_order = " << i_order << std::endl;
	    // 	std::cout << "i_sn = " << i_sn << std::endl;

	    // 	std::cout << "source node location = ";
	    // 	for(unsigned i=0; i<source_node_x.size(); i++)
	    // 	  std::cout << source_node_x[i] << ",";
	    // 	std::cout << std::endl;

	    // 	std::cout << "rough element location = ";
	    // 	for(unsigned i=0; i<x_kn[0][0].size(); i++)
	    // 	  std::cout << x_kn[0][0][i] << ",";
	    // 	std::cout << std::endl;

	    // 	std::cout << "previous estimates: ";
	    // 	for(unsigned i=0; i<temp_bm_prev.size(); i++)
	    // 	  std::cout << temp_bm_prev[i] << ",";
	    // 	std::cout << std::endl;

	    // 	std::cout << "current estimates : ";
	    // 	for(unsigned i=0; i<temp_bm.size(); i++)
	    // 	  std::cout << temp_bm[i] << ",";
	    // 	std::cout << std::endl;

	    // 	std::cout << "reldiff = " << reldiff << std::endl  << std::endl << std::endl;
	    //   }

	    // Go to the next order
	    i_order++;
	  }
	// Repeat while the difference is too large and the order is small enough
	while(( (reldiff>reltol) && (i_order < n_order) )
	      // Always repeat at least once to get a real reldiff value.
	      || (i_order == 0));

	// If we hit the order limit without being accurate enough give an error
	if (i_order >= n_order)
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

  // Given all the values for a certain order, dump them out
  template<class ELEMENT,unsigned DIM>
  void MicromagFaceElement<ELEMENT,DIM>::
  dump_values(const Vector< Vector<double> > &x_kn,
	      const Vector<double> &source_node_x,
	      const Vector<Vector<double> > &normal) const
  {
    // Get the values
    Vector<double> gnd(x_kn.size(),0.0);
    Vector< Vector<double> > s_kn(x_kn.size());
    for(unsigned kn=0; kn<x_kn.size(); kn++)
      {
	gnd[kn] = green_normal_derivative(x_kn[kn],source_node_x,normal[kn]);
      }

    // Output
    std::cout << "For source node at " << source_node_x << std::endl;

    for(unsigned kn=0; kn<x_kn.size(); kn++)
      std::cout << "s = " << s_kn[kn] << ", green normal deriv = " << gnd[kn] << std::endl;

  }


  //============================================================
  /// Get the maximum elementwise difference between two vectors relative to the
  /// appropriate element of v1.
  //============================================================
  template<class ELEMENT,unsigned DIM>
  double MicromagFaceElement<ELEMENT,DIM>::
  max_rel_error(const Vector<double> &v1, const Vector<double> &v2) const
  {
#ifdef PARANOID
    if(v1.size() != v2.size())
      {
	std::cerr << "Vectors are of different sizes." << std::endl;
      }
#endif

    // Get the element-wise difference between the vectors
    Vector<double> diff(v1.size());
    for(unsigned l=0; l<v1.size(); l++)
      diff[l] = std::abs(v1[l] - v2[l]);

    // Find the location of the maximum difference
    Vector<double>::iterator worst_error_it
      = std::max_element(diff.begin(),diff.end());

    // If the worst absolute error is close to zero return zero. This is sort of
    // a hack - if the actual values are tiny then this is wrong, but if the
    // absolute values are tiny things will go wrong anyway and/or the values
    // are unimportant.
    if( v1[worst_error_it - diff.begin()] < 1e-9)
      return 0.0;

    // Divide the maximum difference by the corresponding element of v1.
    return (*worst_error_it) / v1[worst_error_it - diff.begin()];
  }

  //======================================================================
  ///
  //======================================================================
  template<class ELEMENT,unsigned DIM>
  double MicromagFaceElement<ELEMENT,DIM>::
  green_normal_derivative(const Vector<double>& x,
			  const Vector<double>& y,
			  const Vector<double>& n) const
  {
    // Get dimensions
    const unsigned node_dim = nodal_dimension();

#ifdef PARANOID
    // Check that n is a unit vector
    double n_sq = 0.0;
    for(unsigned i=0; i<node_dim; i++)
      n_sq += n[i]*n[i];
    if(std::abs(n_sq - 1.0) > 1e-10)
      {
	throw OomphLibError("n is not a unit vector",
			    "MicromagFaceElement::green_normal_derivative",
			    OOMPH_EXCEPTION_LOCATION);
      }
#endif

    // Calculate the distance r from x to y.
    double r(0.0), r_sq(0.0);
    for(unsigned i=0; i<node_dim; i++)
      r_sq += pow((x[i] - y[i]),2);

    r = sqrt(r_sq);

    //??ds if x is a potentially singular point CHEAT
    //??ds could go horribly wrong, probably fine so long as not actually singular
    //??ds only works in 2d anyway
    if(r < 1e-6)
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

    // Calculate dot product of n with the unit vector in direction r.
    double ndotr(0.0);
    for(unsigned i=0; i<node_dim; i++)
      ndotr += r_unit[i]*n[i];

    if(std::abs(ndotr) < 1e-8)
      return 0.0;

    // std::cout << "ndotr = " << ndotr << std::endl;
    // std::cout << "r = " << r << std::endl;
    // std::cout << "greens/ndotr = " << -1/Pi * std::pow((2*r),-1) * 2 << std::endl;
    // std::cout << "greens = " <<  -1/Pi * ndotr * std::pow((2*r),-1*(dim()-1)) * 2 << std::endl;

    // dgreendn = -n dot r * 1/pi * (1/2)^(node_dim-1) * (1/r)^node_dim
    // See write up for details of calculation.
    double exponent = - (double(node_dim) - 1);
    return (1.0/Pi) * ndotr * std::pow((2.0*r),exponent); //??ds had *2 here for some reason...
  }

  /*
    hybrid FEM/BEM method from:
    T. R. Koehler, D. R. Fredkin, "Finite Element Methods for
    Micromagnetics", IEEE Trans. Magn. 28 (1992) 1239-1244

    integrations for matrix elements calculated according to:
    D. A. Lindholm, "Three-Dimensional Magnetostatic Fields from
    Point-Matched Integral Equations with Linearly Varying Scalar
    Sources", IEEE Trans. Magn. MAG-20 (1984) 2025-2032.
  */
  int Bele(const std::vector<double>& bvert, const std::vector<double>& facv1,
	   const std::vector<double>& facv2, const std::vector<double>& facv3,
	   std::vector<double>& matele)
  {
    using namespace magpar;

    std::vector<double> rr(ND,0.0), zeta(ND,0.0);
    std::vector<double> rho1(ND,0.0),rho2(ND,0.0),rho3(ND,0.0);
    std::vector<double> eta1(ND,0.0),eta2(ND,0.0),eta3(ND,0.0);
    std::vector<double> xi1(ND,0.0),xi2(ND,0.0),xi3(ND,0.0);
    std::vector<double> gamma1(ND,0.0),gamma2(ND,0.0),gamma3(ND,0.0);
    std::vector<double> p(ND,0.0);
    double a(0.0), omega(0.0), eta1l,eta2l,eta3l, s1,s2,s3, rho1l,rho2l,rho3l, zetal;

    matele[0]=matele[1]=matele[2]=0.0;

    /* get coordinates of face's vertices */
    my_dcopy(ND,facv1,1,rho1,1);
    my_dcopy(ND,facv2,1,rho2,1);
    my_dcopy(ND,facv3,1,rho3,1);

    /* calculate edge vectors and store them in xi_j */
    my_dcopy(ND,rho2,1,xi1,1);
    my_daxpy(ND,-1.0,rho1,1,xi1,1);
    my_dcopy(ND,rho3,1,xi2,1);
    my_daxpy(ND,-1.0,rho2,1,xi2,1);
    my_dcopy(ND,rho1,1,xi3,1);
    my_daxpy(ND,-1.0,rho3,1,xi3,1);

    /* calculate zeta direction */
    douter(ND,xi1,xi2,zeta);

    /* calculate area of the triangle */
    zetal=my_dnrm2(ND,zeta,1);
    a=0.5*zetal;

    /* renorm zeta */
    my_dscal(ND,1.0/zetal,zeta,1);

    /* calculate s_j and normalize xi_j */
    s1=my_dnrm2(ND,xi1,1);
    my_dscal(ND,1.0/s1,xi1,1);
    s2=my_dnrm2(ND,xi2,1);
    my_dscal(ND,1.0/s2,xi2,1);
    s3=my_dnrm2(ND,xi3,1);
    my_dscal(ND,1.0/s3,xi3,1);

    douter(ND,zeta,xi1,eta1);
    douter(ND,zeta,xi2,eta2);
    douter(ND,zeta,xi3,eta3);

    gamma1[0]=gamma3[1]=my_ddot(ND,xi2,1,xi1,1);
    gamma1[1]=my_ddot(ND,xi2,1,xi2,1);
    gamma1[2]=gamma2[1]=my_ddot(ND,xi2,1,xi3,1);

    gamma2[0]=gamma3[2]=my_ddot(ND,xi3,1,xi1,1);
    gamma2[2]=my_ddot(ND,xi3,1,xi3,1);

    gamma3[0]=my_ddot(ND,xi1,1,xi1,1);

    /* get R=rr, the location of the source vertex (and the singularity) */
    rr=bvert;

    // If very close to the current surface element (triangle) then result is zero, done.
    double d = PointFromPlane(rr,rho1,rho2,rho3);
    if (fabs(d)<D_EPS) return(0);

    /* calculate rho_j */
    my_daxpy(ND,-1.0,rr,1,rho1,1);
    my_daxpy(ND,-1.0,rr,1,rho2,1);
    my_daxpy(ND,-1.0,rr,1,rho3,1);

    /* zetal gives ("normal") distance of R from the plane of the triangle */
    zetal=my_ddot(ND,zeta,1,rho1,1);

    /* skip the rest if zetal==0 (R in plane of the triangle)
       -> omega==0 and zetal==0 -> matrix entry=0
    */
    if (std::abs(zetal)<=D_EPS) {
      return(0);
    }

    rho1l=my_dnrm2(ND,rho1,1);
    rho2l=my_dnrm2(ND,rho2,1);
    rho3l=my_dnrm2(ND,rho3,1);

    double t_nom,t_denom;
    t_nom=
      rho1l*rho2l*rho3l+
      rho1l*my_ddot(ND,rho2,1,rho3,1)+
      rho2l*my_ddot(ND,rho3,1,rho1,1)+
      rho3l*my_ddot(ND,rho1,1,rho2,1);
    t_denom=
      std::sqrt(2.0*
		(rho2l*rho3l+my_ddot(ND,rho2,1,rho3,1))*
		(rho3l*rho1l+my_ddot(ND,rho3,1,rho1,1))*
		(rho1l*rho2l+my_ddot(ND,rho1,1,rho2,1))
		);

    /* catch special cases where the argument of acos
       is close to -1.0 or 1.0 or even outside this interval

       use 0.0 instead of D_EPS?
       fixes problems with demag field calculation
       suggested by Hiroki Kobayashi, Fujitsu
    */
    if (t_nom/t_denom<-1.0) {
      omega=(zetal >= 0.0 ? 1.0 : -1.0)*2.0*PI;
    }
    /* this case should not occur, but does - e.g. examples1/headfield */
    else if (t_nom/t_denom>1.0) {
      return(0);
    }
    else {
      omega=(zetal >= 0.0 ? 1.0 : -1.0)*2.0*std::acos(t_nom/t_denom);
    }

    eta1l=my_ddot(ND,eta1,1,rho1,1);
    eta2l=my_ddot(ND,eta2,1,rho2,1);
    eta3l=my_ddot(ND,eta3,1,rho3,1);

    p[0]=log((rho1l+rho2l+s1)/(rho1l+rho2l-s1));
    p[1]=log((rho2l+rho3l+s2)/(rho2l+rho3l-s2));
    p[2]=log((rho3l+rho1l+s3)/(rho3l+rho1l-s3));

    matele[0]=(eta2l*omega-zetal*my_ddot(ND,gamma1,1,p,1))*s2/(8.0*PI*a);
    matele[1]=(eta3l*omega-zetal*my_ddot(ND,gamma2,1,p,1))*s3/(8.0*PI*a);
    matele[2]=(eta1l*omega-zetal*my_ddot(ND,gamma3,1,p,1))*s1/(8.0*PI*a);

    return(0);
  }

  //======================================================================
  ///
  //======================================================================
  template<class ELEMENT,unsigned DIM>
  void MicromagFaceElement<ELEMENT,DIM>::
  fill_in_be_contribution_analytic(DenseMatrix<double> &boundary_matrix) const
  {
#ifdef PARANOID
    // Check 3d
    if(nodal_dimension() !=3)
      throw OomphLibError("Analytic calculation of boundary matrix only works in 3D.",
			  "MicromagFaceElement<ELEMENT,DIM>::fill_in_be_contribution_analytic",
			  OOMPH_EXCEPTION_LOCATION);

    if(nnode_1d() != 2)
      throw OomphLibError("Analytic calculation of boundary matrix only works for linear (i.e. flat) elements.",
			  "MicromagFaceElement<ELEMENT,DIM>::fill_in_be_contribution_analytic",
			  OOMPH_EXCEPTION_LOCATION);

    // ??ds check no hanging nodes - not sure what to do there yet
#endif

    // List of the node numbers for the three nodes that we are taking to be the
    // three corners of the triangle.
    Vector<unsigned> l(3,0);
    l[0] = 0; l[1] = 1; l[2] = 2; // Start with the first three nodes

    // Get global nodal positions for first three nodes
    Vector<Vector<double> > x_nds(3);
    x_nds[0].assign(3,0.0); node_pt(0)->position(x_nds[0]);
    x_nds[1].assign(3,0.0); node_pt(1)->position(x_nds[1]);
    x_nds[2].assign(3,0.0); node_pt(2)->position(x_nds[2]);

    //    std::cout << x_nds[0] << " " << x_nds[1] << " " << x_nds[2] << " " << std::endl;

    // Convert element into triangle sub-elements (do nothing if elements are
    // already triangular).
    if(this->nvertex_node() == 3)
      {
	// Evaluate the integral
	analytic_integral_dgreendn_triangle(x_nds, l, boundary_matrix);

	// 	for(unsigned i=0; i<3; i++)
	// 	  {
	// 	    for(unsigned j=0; j<boundary_mesh_pt()->nnode(); j++)
	// 	      std::cout << boundary_matrix(i,j) << " ";
	// 	    std::cout << std::endl;
	// 	  }
	// std::cout << std::endl;
	// std::cout << std::endl;
	// std::cout << std::endl;

      }
    else if(this->nvertex_node() == 4)
      {
	// Evaluate for the triangle formed by the first three nodes
	analytic_integral_dgreendn_triangle(x_nds, l, boundary_matrix);

	// Change a node to the opposite vertex to form a new triangle and
	// repeat
	node_pt(3)->position(x_nds[1]);
	// std::cout << x_nds[0] << " " << x_nds[1] << " " << x_nds[2] << " " << std::endl;
	l[1] = 3; //??ds I think this might go wrong...
	analytic_integral_dgreendn_triangle(x_nds, l, boundary_matrix);
      }
    else
      {
	throw OomphLibError("Unhandled number of vertex nodes in boundary element. Could be the wrong dimension, hanging nodes or higher order elements. Hanging nodes could be implemented but higher order elements or lower dimension cannot (afaik).",
			    "MicromagFaceElement<ELEMENT,DIM>::fill_in_be_contribution_analytic",
			    OOMPH_EXCEPTION_LOCATION);
      }
  }

  //======================================================================
  ///
  //??ds could pass in unit normal since same for all triangle sub-elements
  //======================================================================
  template<class ELEMENT,unsigned DIM>
  void MicromagFaceElement<ELEMENT,DIM>::
  analytic_integral_dgreendn_triangle(const Vector<Vector<double> >& x_tn,
				      const Vector<unsigned>& l,
				      DenseMatrix<double>& boundary_matrix) const
  {
    // Only works in 3D and for triangles (3 nodes)
    const unsigned node_dim = 3;
    const unsigned n_node = 3;

    /* First some pre-calculations to get everything ready. */

    // Calculate the length of the triangle sides and the unit vectors along the
    // triangle sides.
    Vector<double> side_length(n_node,0.0);
    Vector<Vector<double> > side_direction(n_node);
    for(unsigned i=0; i < n_node; i++)
      {
	// Get the next node around the triangle (i.e. 0 -> 1 -> 2 -> 0 etc.).
	unsigned next_node = (i+1)%n_node;

	// Get the vector along this side (xi in Lindholm1984)
	vector_diff(x_tn[next_node],x_tn[i],side_direction[i]);

	// Get the length of this side (s in Lindholm1984)
	side_length[i] = mod(side_direction[i]);
      }

    // Calculate the non-unit normal to triangle. Assuming flat element => same
    // everywhere so we can just take the cross product of any two vectors in
    // the plane. Use the first two edge vectors.
    Vector<double> unit_normal(node_dim,0.0);
    cross(side_direction[0],side_direction[1],unit_normal);

    // Calculate area of triangle using the cross product of the (unnormalised)
    // side_direction vectors, already calculated since it is the normal.
    double area = mod(unit_normal)/2.0;

    // Normalise the unit normal and side direction vectors now that we have the
    // area.
    normalise(unit_normal);
    for(unsigned i=0; i<n_node; i++) normalise(side_direction[i]);

    // Calculate gamma
    Vector< Vector<double> > gamma(n_node);
    for(unsigned i=0; i<3; i++)
      {
	gamma[i].assign(n_node,0.0); // Initialise gamma[i]
	unsigned next_node = (i+1)%n_node; // Get next triangle vertex
	for(unsigned j=0; j<n_node; j++)
	  gamma[i][j] = dot(side_direction[next_node],side_direction[j]);
      }

    /* Now evaluate the integral for every boundary node and add to
       boundary matrix. */
    for(unsigned i_sn=0; i_sn < boundary_mesh_pt()->nnode(); i_sn++)
      {
	// Get position of this node
	Vector<double> x_sn(node_dim,0.0);
	boundary_mesh_pt()->node_pt(i_sn)->position(x_sn);

	// Calculate rho (vector from each node in the triangle to source node)
	// and it's length.
	Vector<Vector<double> > rho(n_node); Vector<double> rhol(n_node);
	for(unsigned i_nd=0; i_nd<n_node; i_nd++)
	  {
	    vector_diff(x_tn[i_nd],x_sn,rho[i_nd]);
	    rhol[i_nd] = mod(rho[i_nd]);
	  }

	// Calculate zeta: the distance between the element and the source node
	// in the normal direction to the element.
	double zeta = dot(unit_normal,rho[0]);

	// If source node is in the plane of the element (i.e. zeta ~ 0) then
	// n.r = 0, nothing to calculate or add so we can move on to the next
	// source node.
	double tol = 1e-8;
	if( (std::abs(zeta) < tol) )
	  continue;

	// Calculate "P" (see paper Lindholm1984) for each node in the triangle
	Vector<double> P(n_node,0.0);
	for(unsigned i=0; i<n_node; i++)
	  {
	    unsigned next_node = (i+1)%n_node;
	    P[i] = std::log( (rhol[i] + rhol[next_node] + side_length[i])
			     /(rhol[i] + rhol[next_node] - side_length[i]) );
	  }


	// Calculate the solid angle (see Lindholm 1984)
	// I think it's impossible for the denominator to be zero since that
	// would imply either a zero rho vector or two anti-parallel rho vectors
	// which only happens for a source within the triangle. Both of these
	// cases would give zeta < 1e-8.
	double ratio = (rhol[0]*rhol[1]*rhol[2]
			+ rhol[0] * dot(rho[1],rho[2])
			+ rhol[1] * dot(rho[0],rho[2])
			+ rhol[2] * dot(rho[1],rho[0]))
	  /
	  std::sqrt(2
		    *(rhol[1]*rhol[2] + dot(rho[1],rho[2]) )
		    *(rhol[1]*rhol[0] + dot(rho[1],rho[0]) )
		    *(rhol[0]*rhol[2] + dot(rho[0],rho[2]) )
		    );

	// Round-off errors can cause ratio to be out of range for inverse cos
	// so we need to check it.
	if(ratio > 1.0) ratio = 1.0;
	else if(ratio < -1.0) ratio = -1.0;

	int sign = (zeta > 0.0 ? +1.0 : -1.0);
	double omega = sign * 2 * std::acos(ratio);

	// Calculate eta: the unit vector normal to the side. Use it to
	// calculate etal: the distance in the diretion of eta to the source
	// node.
	Vector<Vector<double> > eta(3); Vector<double> etal(3,0.0);
	for(unsigned i=0; i<n_node; i++)
	  {
	    eta[i].assign(3,0.0);
	    cross(unit_normal,side_direction[i],eta[i]);
	    etal[i] = dot(eta[i],rho[i]);
	  }

	// //??temp calculate with magpar code for testing
	// Vector<double> magpar_contribution(3,0.0);
	// Bele(x_sn, x_tn[0], x_tn[1], x_tn[2], magpar_contribution);
	// Vector<double> oomph_contribution(3,0.0);

	/* Now put it all together and add the contribution to the boundary
	   element matrix */
	for(unsigned i_tn=0; i_tn < n_node; i_tn++)
	  {
	    unsigned next_node = (i_tn+1)%n_node;

	    // Add contribution to the appropriate value in the boundary matrix
	    boundary_matrix(l[i_tn],i_sn) += (side_length[next_node]/(8*Pi*area))
	      *((etal[next_node] * omega) - (zeta * dot(gamma[i_tn],P)));

	    // Note that Lindholm does include the factor of 1/(4*pi) in his
	    // operator. Also note that he does NOT include the solid angle
	    // factor on surfaces, the solid angle factor above is unrelated.
	  }

      }
  }


  //  //======================================================================
  //   /// Point element to sit at sharp corners and add the angle/solid angle of the
  //   /// corner to the boundary element matrix. DIM is the dimension of the entire
  //   /// problem not the dimension of the point element (which is always 0).
  //   //======================================================================
  //   template<class ELEMENT, unsigned DIM>
  //   class MicromagCornerAngleElement :
  //     public virtual FaceGeometry<FaceGeometry<ELEMENT> >,
  //     public virtual PointElement
  //   {};


  //   //======================================================================
  //   /// In 1D there are no sharp corners so no 1D version is needed.
  //   //======================================================================

  //   //======================================================================
  //   /// 2D specialisation: calculate angles.
  //   //======================================================================
  //   template< class ELEMENT>
  //   class MicromagCornerAngleElement<ELEMENT,2> :
  //     public virtual FaceGeometry<FaceGeometry<ELEMENT> >,
  //     public virtual PointElement
  //   {
  //   private:

  //     /// Pointer to the face elements to which this point element is
  //     /// attached.
  //     Vector<FaceGeometry<ELEMENT>* > Face_element_pt;

  //     /// The index of the node in face elements to which this point element is
  //     /// attached.
  //     Vector<unsigned> Node_number;

  //   public:

  //     MicromagCornerAngleElement() : Face_element_pt(2)
  //     {
  //       // Face_element_pt[0] = e1;
  //       // Face_element_pt[1] = e2;

  //       // Node_number[0] = n1;
  //       // Node_number[1] = n2;
  //     }

  //     /// Calculate the angle between the two attached face elements
  //     double calculate_corner_fractional_angle() const
  //     {
  //       Vector<Vector<double> > t;

  //       // For each attached face element (two of them) get the tangent vector
  //       for(unsigned fe=0; fe<2; fe++)
  // 	{
  // 	  // Find out the number of nodes in the face element
  // 	  unsigned n_node_face = Face_element_pt[fe]->nnode();

  // 	  // Get the value of the shape function derivatives at the node
  // 	  Shape psi(n_node_face); // We have to calculate shape functions as well...
  // 	  DShape dpsids(n_node_face,1);
  // 	  Face_element_pt[fe]->dshape_local(s_face,psi,dpsids);

  // 	  // Calculate all derivatives of the spatial coordinates wrt local
  // 	  // coordinates
  // 	  Vector<double> interpolated_dxds(2,0.0);
  // 	  // for(unsigned j=0;j<2;j++)
  // 	  //   {
  // 	  //     interpolated_dxds[j] +=
  // 	  // 	Face_element_pt[fe]->nodal_position(l,j) * dpsids(l,0);
  // 	  //   }

  // 	  // Add to list of tangents
  // 	  t.push_back(interpolated_dxds);
  // 	}

  //       // Calculate the angle between them (inverse cos of the dot product).
  //       return acos(t[0][0]*t[1][0] + t[0][1]*t[1][1]) / (2*Pi);
  //     }

  //   };

  //   //======================================================================
  //   /// 3D specialisation: calculate solid angles.
  //   //======================================================================
  //   template< class ELEMENT>
  //   class MicromagCornerAngleElement<ELEMENT,3> :
  //     public virtual FaceGeometry<FaceGeometry<ELEMENT> >,
  //     public virtual PointElement
  //   {

  //     //??ds calculating the angle is going to be harder here, do it later...
  //   };



}

#endif
