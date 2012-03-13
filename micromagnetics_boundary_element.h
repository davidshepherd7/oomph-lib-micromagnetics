#ifndef OOMPH_MICROMAGNETICS_BOUNDARY_ELEMENT_H
#define OOMPH_MICROMAGNETICS_BOUNDARY_ELEMENT_H


#include "generic.h"
#include "./micromagnetics_element.h"
#include "./micromagnetics_element.cc"
#include "./variable_order_quadrature.h"
#include <functional>

using namespace oomph;
using namespace MathematicalConstants;

void print_vec(std::ostream& out, const Vector<double> &vec)
{
  // std::copy(vec.begin(), vec.end(),
  // 	    ostream_iterator<double>(out, ","));
  for(unsigned i=0; i<vec.size(); i++)
    out << vec[i] << " ";
}


namespace oomph
{

  //===========================================================================
  ///
  //===========================================================================
  template<class ELEMENT,unsigned DIM>
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
				 VShape &psi, VShape &test) const
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
					 VShape &psi, VShape &test) const
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
    QVariableOrderGaussLegendre<el_dim>* v_int_pt =
      dynamic_cast<QVariableOrderGaussLegendre<el_dim>* >(integral_pt());
#else
    QVariableOrderGaussLegendre<el_dim>* v_int_pt =
      static_cast<QVariableOrderGaussLegendre<el_dim>*>(integral_pt());
#endif

    // Start of adaptive integration scheme
    //====================================

    // Setup error parameters for adaptive integration
    double reltol = 1e-4, reldiff;

    // Set up storage to keep track of orders used:
    std::vector<unsigned> order_list = {2,4,8,16,32,64,128,256,};
    unsigned n_order = order_list.size();
    Vector<unsigned> order_is_calculated(n_order,0);

    // Set up storage for data at each knot for each order
    Vector< Vector< Vector<double> > > x_kn(n_order), normal(n_order);
    Vector< Vector<double> > jw(n_order);
    Vector< Vector<VShape> > test(n_order), psi(n_order);

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
	    //??ds do the initilisation better?
	    temp_bm_prev = temp_bm;
	    for(unsigned l=0; l<n_element_node; l++)
	      temp_bm[l] = 0.0;

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
		    VShape psi_local(n_element_node);
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
    std::cout << "For source node at ";
    print_vec(std::cout,source_node_x);
    std::cout << std::endl;

    for(unsigned kn=0; kn<x_kn.size(); kn++)
      {
	print_vec(std::cout,s_kn[kn]);
	std::cout << gnd[kn] << std::endl;
      }

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
      = max_element(diff.begin(),diff.end());

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
    return -1/Pi * ndotr * std::pow((2*r),(dim()-1)) * 2;
  }




  //======================================================================
  /// Point element to sit at sharp corners and add the angle/solid angle of the
  /// corner to the boundary element matrix. DIM is the dimension of the entire
  /// problem not the dimension of the point element (which is always 0).
  //======================================================================
  template<class ELEMENT, unsigned DIM>
  class MicromagCornerAngleElement :
    public virtual FaceGeometry<FaceGeometry<ELEMENT> >,
    public virtual PointElement
  {};


  //======================================================================
  /// In 1D there are no sharp corners so no 1D version is needed.
  //======================================================================

  //======================================================================
  /// 2D specialisation: calculate angles.
  //======================================================================
  template< class ELEMENT>
  class MicromagCornerAngleElement<ELEMENT,2> :
    public virtual FaceGeometry<FaceGeometry<ELEMENT> >,
    public virtual PointElement
  {
  private:

    /// Pointer to the face elements to which this point element is
    /// attached.
    Vector<FaceGeometry<ELEMENT>* > Face_element_pt;

    /// The index of the node in face elements to which this point element is
    /// attached.
    Vector<unsigned> Face_node_number;

  public:

    MicromagCornerAngleElement() : Face_element_pt(2)
    {
      // Face_element_pt[0] = e1;
      // Face_element_pt[0] = e2;

      // Face_node_number.push_back(...);
    }

    /// Calculate the angle between the two attached face elements
    double calculate_corner_fractional_angle() const
    {
      Vector<Vector<double> > t;

      // For each attached face element (two of them) get the tangent vector
      for(unsigned fe=0; fe<2; fe++)
	{
	  // Find out the corresponding node in the face element
	  unsigned l = Face_node_number[fe];

	  // Find out the number of nodes in the face element
	  unsigned n_node_face = Face_element_pt[fe]->nnode();

	  // Get the local coordinate of this node in the face element
	  Vector<double> s_face(1,0.0);
	  Face_element_pt[fe]->local_coordinate_of_node(l,s_face);

	  // Get the value of the shape function derivatives at the node
	  Shape psi(n_node_face); // We have to calculate shape functions as well...
	  DShape dpsids(n_node_face,1);
	  Face_element_pt[fe]->dshape_local(s_face,psi,dpsids);

	  // Calculate all derivatives of the spatial coordinates wrt local
	  // coordinates
	  Vector<double> interpolated_dxds(2,0.0);
	  for(unsigned j=0;j<2;j++)
	    {
	      interpolated_dxds[j] +=
		Face_element_pt[fe]->nodal_position(l,j) * dpsids(l,0);
	    }

	  // Add to list of tangents
	  t.push_back(interpolated_dxds);
	}

      // Calculate the angle between them (inverse cos of the dot product).
      return acos(t[0][0]*t[1][0] + t[0][1]*t[1][1]);
    }

  };

  //======================================================================
  /// 3D specialisation: calculate solid angles.
  //======================================================================
  template< class ELEMENT>
  class MicromagCornerAngleElement<ELEMENT,3> :
    public virtual FaceGeometry<FaceGeometry<ELEMENT> >,
    public virtual PointElement
  {

    //??ds calculating the angle is going to be harder here, do it later...
  };

}

#endif
