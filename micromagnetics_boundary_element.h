#ifndef OOMPH_MICROMAGNETICS_BOUNDARY_ELEMENT_H
#define OOMPH_MICROMAGNETICS_BOUNDARY_ELEMENT_H


#include "generic.h"
#include "../micromagnetics_element.h"
#include "../micromagnetics_element.cc"

using namespace oomph;
using namespace MathematicalConstants;

namespace oomph
{

  //================================================================================
  ///
  //================================================================================
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

    /// Add the element's contribution to its residual vector (dummy function)
    inline void fill_in_contribution_to_residuals(Vector<double> &dummy)
    {
      // Do nothing - no residuals to add
    }

    /// \short Add the element's contribution to its residual vector and its
    /// Jacobian matrix
    inline void fill_in_contribution_to_jacobian(Vector<double> &dummy,
						 DenseMatrix<double> &boundary_matrix)
    {
      fill_in_boundary_element_contribution_micromag(boundary_matrix);
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
    void fill_in_boundary_element_contribution_micromag(DenseMatrix<double> &boundary_matrix)
      const;

    /// \short Pointer to the boundary mesh (needed to access nodes outside of this element
    /// for calculation of boundary matrix).
    Mesh* Mesh_pt;

    /// The dimension of the element surface/volume (i.e. one less than the dimension of the nodes.
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

    //??ds Set integral pointer to point to adaptive method
    //set_integration_scheme(....);

    // // Cast to the appropriate equation element so that we can
    // // find the index at which the variable is stored
    // //??ds this code seems horrible....
    // switch(Node_dim)
    //   {
    // 	//One dimensional problem
    //   case 1:
    // 	{
    // 	  MicromagEquations<1>* eqn_pt = dynamic_cast<MicromagEquations<1>*>(bulk_el_pt);
    // 	  //If the cast has failed die
    // 	  if(eqn_pt==0)
    // 	    {
    // 	      throw OomphLibError("Cannot cast the bulk element.",
    // 				  "MicromagFaceElement::MicromagFaceElement()",
    // 				  OOMPH_EXCEPTION_LOCATION);
    // 	    }
    // 	  else
    // 	    {
    // 	      // Read the indicies from the (cast) bulk element
    // 	      Phi_1_index_micromag = eqn_pt->phi_1_index_micromag();
    // 	      Phi_2_index_micromag = eqn_pt->phi_2_index_micromag();
    // 	    }
    // 	}
    // 	break;

    // 	//Two dimensional problem
    //   case 2:
    // 	{
    // 	  MicromagEquations<2>* eqn_pt = dynamic_cast<MicromagEquations<2>*>(bulk_el_pt);
    // 	  //If the cast has failed die
    // 	  if(eqn_pt==0)
    // 	    {
    // 	      throw OomphLibError("Cannot cast the bulk element.",
    // 				  "MicromagFaceElement::MicromagFaceElement()",
    // 				  OOMPH_EXCEPTION_LOCATION);
    // 	    }
    // 	  else
    // 	    {
    // 	      // Read the indicies from the (cast) bulk element
    // 	      Phi_1_index_micromag = eqn_pt->phi_1_index_micromag();
    // 	      Phi_2_index_micromag = eqn_pt->phi_2_index_micromag();
    // 	    }
    // 	}
    // 	break;

    // 	//Three dimensional problem
    //   case 3:
    // 	{
    // 	  MicromagEquations<3>* eqn_pt = dynamic_cast<MicromagEquations<3>*>(bulk_el_pt);
    // 	  //If the cast has failed die
    // 	  if(eqn_pt==0)
    // 	    {
    // 	      throw OomphLibError("Cannot cast the bulk element.",
    // 				  "MicromagFaceElement::MicromagFaceElement()",
    // 				  OOMPH_EXCEPTION_LOCATION);
    // 	    }
    // 	  else
    // 	    {
    // 	      // Read the indicies from the (cast) bulk element
    // 	      Phi_1_index_micromag = eqn_pt->phi_1_index_micromag();
    // 	      Phi_2_index_micromag = eqn_pt->phi_2_index_micromag();
    // 	    }
    // 	}
    // 	break;

    // 	//Any other case is an error
    //   default:
    // 	std::ostringstream error_stream;
    // 	error_stream <<  "Dimension of node is " << Node_dim
    // 		     << ". It should be 1,2, or 3!" << std::endl;

    // 	throw OomphLibError(error_stream.str(),
    // 			    "MicromagFaceElement::MicromagFaceElement()",
    // 			    OOMPH_EXCEPTION_LOCATION);
    // 	break;
    //   }
  }

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
  fill_in_boundary_element_contribution_micromag(DenseMatrix<double> &boundary_matrix)
    const
  {
    // Find out dimension of element
    const unsigned el_dim = Node_dim - 1;

    //Find out how many nodes there are
    const unsigned n_element_node = nnode();
    // std::cout << n_element_node << std::endl;

    //Set up memory for the shape and test functions
    Shape psi(n_element_node), test(n_element_node);

    // // Get current time
    // double time = time_pt()->time();

    //??ds adaptive quadrature: given acceptable error choose integration
    // method and return integral_pt

    //Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
      {
	// Get values of s (local coordinate)
	Vector<double> s(el_dim,0.0);
	for(unsigned j=0; j<el_dim; j++) {s[j] = integral_pt()->knot(ipt,j);}

	//Get the integral weight
	double w = integral_pt()->weight(ipt);

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
	for(unsigned i_source_node=0; i_source_node<n_boundary_node; i_source_node++)
	  {
	    // Get coordinates of source node
	    Vector<double> source_node_x(Node_dim,0.0);
	    mesh_pt()->node_pt(i_source_node)->position(source_node_x);

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
		boundary_matrix(l,i_source_node) -= dgreendn * test(l) * W;
		// boundary_matrix(l,i_source_node) += 1;
	      }

	  }
      }

    //??ds need to seperately add the contribution at each node from angles

    //??ds probably need to seperately calculate for the elements near the current one eventually

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
