#ifndef OOMPH_HYBRID_BOUNDARY_ELEMENT_DRIVER_H
#define OOMPH_HYBRID_BOUNDARY_ELEMENT_DRIVER_H

/*
  ??ds  description of file goes here
*/

#include "generic.h"
#include "../micromagnetics_element.h"

using namespace oomph;

namespace oomph
{

  // //============================================================
  // /// A class for the maths used in MagnetostaticBE.
  // //============================================================
  // template<unsigned DIM>
  // class MagnetostaticBEEquations : public virtual FiniteElement
  // {
  // private:

  // public:

  //   /// Constructor
  //   MagnetostaticBEEquations(){}

  //   /// Broken copy constructor
  //   MagnetostaticBEEquations(const MagnetostaticBEEquations& dummy)
  //   {
  //     BrokenCopy::broken_copy("MagnetostaticBEEquations");
  //   }

  //   /// Broken assignment operator
  //   void operator=(const MagnetostaticBEEquations&)
  //   {
  //     BrokenCopy::broken_assign("MagnetostaticBEEquations");
  //   }

  //   /// Self-test: Return 0 for OK.
  //   unsigned self_test()
  //   {return 0;} //??ds write a real test sometime

  //   /// Specify nodal index for phi_1
  //   unsigned phi_1_nodal_index_magnetostaticbe() const {return 0;}

  //   /// Specify nodal index for phi_2
  //   unsigned phi_2_nodal_index_magnetostaticbe() const {return 1;}

  //   /// Return the i-th value stored at local node n, DO take hanging nodes into account
  //   double nodal_value(const unsigned &n, const unsigned &i) const
  //   {return node_pt(n)->value(i);}

  //   /// Return the i-th value at time t stored at local node n, DO take hanging nodes into account
  //   double nodal_value(const unsigned &t, const unsigned &n, const unsigned &i) const
  //   {return node_pt(n)->value(t,i);}

  //   /// Output function at default number of plot points
  //   void output(std::ostream &outfile)
  //   {
  //     const unsigned n_plot = 5;
  //     output(outfile,n_plot);
  //   }

  //   /// Output function at n_plot^DIM points
  //   void output(std::ostream &outfile, const unsigned &n_plot);

  //   /// Get error by comparing with exact solution and get norm of exact solution.
  //   void compute_error(std::ostream &outfile,
  // 		       FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
  // 		       const double& time, double& error, double& norm);

  //   /// Add the element's contribution to its residual vector (wrapper)
  //   void fill_in_contribution_to_residuals(Vector<double> &residuals)
  //   {
  //     // Call the generic residuals function with flag set to 0 using a dummy matrix argument
  //     fill_in_generic_residual_contribution_micromag
  // 	(residuals, GeneralisedElement::Dummy_matrix, 0);
  //   }

  //   // ??ds re-activate this when you're confident that the Jacobian is correct!
  //   // /// Add the element's contribution to its residual vector and element Jacobian matrix (wrapper)
  //   // void fill_in_contribution_to_jacobian(Vector<double> &residuals,
  //   //                                       DenseMatrix<double> &jacobian)
  //   // {
  //   //   //Call the generic routine with the flag set to 1
  //   //   fill_in_generic_residual_contribution_micromag(residuals,jacobian,1);
  //   // }

  //   /// Fill in contribution to residuals and jacobian (if flag is set) from these equations (compatible with multiphysics)
  //   void fill_in_generic_residual_contribution_micromag(Vector<double> &residuals,
  // 							DenseMatrix<double> &jacobian,
  // 							const unsigned& flag);

  // protected:

  //   /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
  //   virtual double dshape_and_dtest_eulerian_micromag(const Vector<double> &s,
  // 						      Shape &psi,
  // 						      DShape &dpsidx,
  // 						      Shape &test,
  // 						      DShape &dtestdx) const=0;

  //   /// Shape, test functions & derivs. w.r.t. to global coords. at integration point ipt. Return Jacobian.
  //   virtual double dshape_and_dtest_eulerian_at_knot_micromag(const unsigned& ipt,
  // 							      Shape &psi,
  // 							      DShape &dpsidx,
  // 							      Shape &test,
  // 							      DShape &dtestdx) const=0;

  // }; // End of MagnetostaticBEEquations class

  // /// Fill in contribution to residuals and jacobian (if flag is set) from these equations (compatible with multiphysics)
  // template<unsigned DIM>
  // void MagnetostaticBEEquations<DIM>::
  // fill_in_generic_residual_contribution_micromag(Vector<double> &residuals,
  // 						 DenseMatrix<double> &jacobian,
  // 						 const unsigned& flag)
  // {
  //   //??ds fill in this bit
  // }




  //================================================================================
  /// A class combining the MicromagBE equations with a face element.
  //================================================================================
  template<class GEOMELEMENT>
  class MicromagFaceElement : public virtual FaceElement,
			      public virtual FaceGeometry<GEOMELEMENT>

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

    /// Add the element's contribution to its residual vector
    inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
      //Call the generic residuals function with flag set to 0
      //using a dummy matrix argument
      fill_in_generic_residual_contribution_micromag_face(residuals,
							  GeneralisedElement::Dummy_matrix,0);
    }


    /// \short Add the element's contribution to its residual vector and its
    /// Jacobian matrix
    inline void fill_in_contribution_to_jacobian(Vector<double> &residuals,
						 DenseMatrix<double> &jacobian)
    {
      //Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution_micromag_face(residuals,jacobian,1);
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
    void fill_in_generic_residual_contribution_micromag_face(Vector<double> &residuals,
							     DenseMatrix<double> &jacobian,
							     const unsigned& flag);

    /// The dimension of the element surface/volume (i.e. one less than the dimension of the nodes.
    unsigned Node_dim;

    /// The index at which phi_1 is stored
    unsigned Phi_1_index_micromag;

    /// The index at which phi_2 is stored
    unsigned Phi_2_index_micromag;

  };

  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, the
  /// index of the fixed local coordinate and its value represented
  /// by an integer (+/- 1), indicating that the face is located
  /// at the max. or min. value of the "fixed" local coordinate
  /// in the bulk element.
  //===========================================================================
  template<class GEOMELEMENT>
  MicromagFaceElement<GEOMELEMENT>::
  MicromagFaceElement(FiniteElement* const &bulk_el_pt,
		      const int &face_index) :
    FaceGeometry<GEOMELEMENT>(), FaceElement()
  {
    // Let the bulk element build the FaceElement, i.e. setup the pointers
    // to its nodes (by referring to the appropriate nodes in the bulk
    // element), etc.
    bulk_el_pt->build_face_element(face_index,this);

    // Extract the nodal dimension of the problem from the dimension of
    // the first (face) node.
    Node_dim = this->node_pt(0)->ndim();

    // Cast to the appropriate equation element so that we can
    // find the index at which the variable is stored
    //??ds this code seems horrible....
    switch(Node_dim)
      {
	//One dimensional problem
      case 1:
	{
	  MicromagEquations<1>* eqn_pt = dynamic_cast<MicromagEquations<1>*>(bulk_el_pt);
	  //If the cast has failed die
	  if(eqn_pt==0)
	    {
	      throw OomphLibError("Cannot cast the bulk element.",
				  "MicromagFaceElement::MicromagFaceElement()",
				  OOMPH_EXCEPTION_LOCATION);
	    }
	  else
	    {
	      // Read the indicies from the (cast) bulk element
	      Phi_1_index_micromag = eqn_pt->phi_1_index_micromag();
	      Phi_2_index_micromag = eqn_pt->phi_2_index_micromag();
	    }
	}
	break;

	//Two dimensional problem
      case 2:
	{
	  MicromagEquations<2>* eqn_pt = dynamic_cast<MicromagEquations<2>*>(bulk_el_pt);
	  //If the cast has failed die
	  if(eqn_pt==0)
	    {
	      throw OomphLibError("Cannot cast the bulk element.",
				  "MicromagFaceElement::MicromagFaceElement()",
				  OOMPH_EXCEPTION_LOCATION);
	    }
	  else
	    {
	      // Read the indicies from the (cast) bulk element
	      Phi_1_index_micromag = eqn_pt->phi_1_index_micromag();
	      Phi_2_index_micromag = eqn_pt->phi_2_index_micromag();
	    }
	}
	break;

	//Three dimensional problem
      case 3:
	{
	  MicromagEquations<3>* eqn_pt = dynamic_cast<MicromagEquations<3>*>(bulk_el_pt);
	  //If the cast has failed die
	  if(eqn_pt==0)
	    {
	      throw OomphLibError("Cannot cast the bulk element.",
				  "MicromagFaceElement::MicromagFaceElement()",
				  OOMPH_EXCEPTION_LOCATION);
	    }
	  else
	    {
	      // Read the indicies from the (cast) bulk element
	      Phi_1_index_micromag = eqn_pt->phi_1_index_micromag();
	      Phi_2_index_micromag = eqn_pt->phi_2_index_micromag();
	    }
	}
	break;

	//Any other case is an error
      default:
	std::ostringstream error_stream;
	error_stream <<  "Dimension of node is " << Node_dim
		     << ". It should be 1,2, or 3!" << std::endl;

	throw OomphLibError(error_stream.str(),
			    "MicromagFaceElement::MicromagFaceElement()",
			    OOMPH_EXCEPTION_LOCATION);
	break;
      }
  }

  //======================================================================
  /// Compute element residual Vector and/or element Jacobian matrix
  ///
  /// flag=1: compute both
  /// flag=0: compute only residual Vector
  ///
  /// Pure version without hanging nodes
  //======================================================================
  template<class ELEMENT>
  void MicromagFaceElement<ELEMENT>::
  fill_in_generic_residual_contribution_micromag_face(Vector<double> &residuals,
						      DenseMatrix<double> &jacobian,
						      const unsigned& flag)
  {
    // Find out dimension of element
    const unsigned el_dim = Node_dim - 1;

    //Find out how many nodes there are
    const unsigned n_node = nnode();

    //Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node,el_dim), dtestdx(n_node,el_dim);

    // Set up vector to store the local coordinates
    Vector<double> s(el_dim);

    // Get current time
    double time = time_pt()->time();

    //Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();
  }

} // End of oomph namespace

int main()
{
  return 0;
}

#endif
