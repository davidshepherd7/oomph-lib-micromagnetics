#ifndef OOMPH_HYBRID_BOUNDARY_ELEMENT_DRIVER_H
#define OOMPH_HYBRID_BOUNDARY_ELEMENT_DRIVER_H

/*
  ??ds  description of file goes here
*/

//??ds put in timestepper...

#include "generic.h"
#include "../micromagnetics_element.cc"
#include "meshes/rectangular_quadmesh.h"

using namespace oomph;
using namespace MathematicalConstants;

namespace oomph
{

  //================================================================================
  /// A class combining the MicromagBE equations with a face element.
  //================================================================================
  template<class ELEMENT>
  class MicromagFaceElement : public virtual FaceGeometry<ELEMENT>,
			      public virtual FaceElement
  {

  public:

    /// \short Constructor, takes the pointer to the bulk element and the
    /// index of the face to which the element is attached.
    MicromagFaceElement(FiniteElement* const &bulk_el_pt,
			const int& face_index,
			Mesh* mesh_pt);

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

    /// Calculate and add in the boundary element matrix contribution for this element
    void fill_in_elemental_contribution_to_boundary_element_matrix(DenseMatrix<double> &boundary_matrix) const;

    /// Add the element's contribution to its residual vector
    inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
      //Call the generic residuals function with flag set to 0
      //using a dummy matrix argument
      fill_in_generic_residual_contribution_micromag(residuals,
						     GeneralisedElement::Dummy_matrix,0);
    }

    /// \short Add the element's contribution to its residual vector and its
    /// Jacobian matrix
    inline void fill_in_contribution_to_jacobian(Vector<double> &residuals,
						 DenseMatrix<double> &jacobian)
    {
      //Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution_micromag(residuals,jacobian,1);
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
    void fill_in_generic_residual_contribution_micromag(Vector<double> &residuals,
							     DenseMatrix<double> &jacobian,
							     const unsigned& flag);

    /// \short Pointer to the mesh (needed to access nodes outside of this element
    /// for calculation of boundary matrix.
    Mesh* Mesh_pt;

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
  template<class ELEMENT>
  MicromagFaceElement<ELEMENT>::
  MicromagFaceElement(FiniteElement* const &bulk_el_pt, const int &face_index, Mesh* mesh_pt)
    : FaceGeometry<ELEMENT>(), FaceElement()
  {
    // Store the mesh pointer (needed to get access to positions of nodes not in element).
    Mesh_pt = mesh_pt;

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
  /// Dummy ??ds
  //======================================================================
  template<class ELEMENT>
  void MicromagFaceElement<ELEMENT>::
  fill_in_generic_residual_contribution_micromag(Vector<double> &residuals,
						      DenseMatrix<double> &dummy,
						      const unsigned& flag)
  {

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

    // Calculate the unit vector in direction r = x - y
    // and the distance from x to y.
    Vector<double> r_unit(Node_dim,0.0);
    double r(0.0);
    for(unsigned i=0; i<Node_dim; i++)
      {
	r += pow((x[i] - y[i]),2);
	r_unit[i] = (x[i] - y[i])/r;
      }

#ifdef PARANOID
    //??ds check that r is not ~machine error
#endif

    // Calculate dot product of n with the unit vector in direction r.
    double ndotr(0.0);
    for(unsigned i=0; i<Node_dim; i++)
      ndotr += r_unit[i]*n[i];

    std::cout << "ndotr = " << ndotr << std::endl;
    std::cout << "r = " << r << std::endl;
    std::cout << "greens/ndotr = " << -1/Pi * pow((2*r),-2) * 2 << std::endl;
    // dgreendn = -n dot r * 1/pi * (1/2)^(node_dim-1) * (1/r)^node_dim
    // See write up for details of calculation.
    return -1/Pi * ndotr * pow((2*r),-2) * 2;


    //??ds for testing
    //return -r;
  }

  //=======================================================================
  /// Compute the effect of phi_1 in this element on ALL nodes on the boundary.
  /// Note that this depends on all other elements but not on the values of phi_1.
  //??ds may have confused test and "psi" functions but it doesn't matter for now
  // since they are the same (Galerkin method)
  //=======================================================================
  template<class ELEMENT>
  void MicromagFaceElement<ELEMENT>::
  fill_in_elemental_contribution_to_boundary_element_matrix(DenseMatrix<double> &boundary_matrix)
    const
  {
    // Find out dimension of element
    const unsigned el_dim = Node_dim - 1;

    //Find out how many nodes there are
    const unsigned n_node = nnode();

    //Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);

    // Set up memory for the coordinate vectors and normal vector
    Vector<double> s(el_dim,0.0), interpolated_x(Node_dim,0.0),
      normal(Node_dim,0.0), target_node_x(Node_dim,0.0);

    // // Get current time
    // double time = time_pt()->time();

    //Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
      {
	//Get the integral weight
	double w = integral_pt()->weight(ipt);

	//Call the derivatives of the shape and test functions
	double J = shape_and_test(s,psi,test);

	//Premultiply the weights and the Jacobian
	double W = w*J;

	// Get values of s (local coordinate)
	for(unsigned j=0; j<el_dim; j++) {s[j] = integral_pt()->knot(ipt,j);}

	// Get values of x (global coordinate)
	for(unsigned l=0; l<n_node; l++)
	  {
	    for(unsigned i=0; i<Node_dim; i++)
	      interpolated_x[i] += nodal_position(l,i)*psi[l];
	  }

	// Compute the normal vector
	//??ds not sure hwo this works - might not work for curved boundaries?
	outer_unit_normal(ipt,normal);

	// Loop over ALL nodes on all boundaries (except the current ones?) (target_node)

	unsigned n_boundary = Mesh_pt->nboundary();
	for(unsigned i_boundary=2; i_boundary<3; i_boundary++)
	  {
	    unsigned n_boundary_nodes = Mesh_pt->nboundary_node(i_boundary);
	    for(unsigned i_target_node=0; i_target_node<n_boundary_nodes; i_target_node++)
	      {
		// Get coordinates of target node
		Mesh_pt->boundary_node_pt(i_boundary,i_target_node)->position(target_node_x);

		std::cout << "Gauss point at: (" << interpolated_x[0]
			  << ", " << interpolated_x[1] << ")"
			  << ". Target node at: (" << target_node_x[0]
			  << ", " << target_node_x[1] << ")"
			  << " on boundary " << i_boundary << std::endl;

		// Calculate dGreendn between target node and integration point
		double dgreendn = green_normal_derivative(interpolated_x,target_node_x,normal);

		// Loop over test functions i.e. (local nodes) adding contributions
		for(unsigned l=0; l<n_node; l++)
		  {
		    // Add contribution to integral (note dGreendn is negative in our definition)
		    // See write up for details
		    // Note that l is the LOCAL node number, need to do local to global mapping elsewhere, as seen in residual construction.
		    boundary_matrix(l,i_target_node) -= dgreendn * test(l) * W;
		  }
	      }
	  }
      }

    //??ds need to seperately add the contribution at each node from angles

    //??ds probably need to seperately calculate for the elements near the current one eventually

  }

  template<class ELEMENT>
  class TwoDBoundaryProblem : public Problem
  {

  public:

    /// Constructor: Pass number of elements and pointer to source function
    TwoDBoundaryProblem(const unsigned& n_x, const unsigned& n_y);

    /// Destructor (empty -- all the cleanup is done in the base class)
    ~TwoDBoundaryProblem(){};

    /// Doc the solution
    void doc_solution(DocInfo& doc_info);

    /// Get number of bulk elements
    unsigned get_n_bulk_element() {return N_bulk_element;}

  private:

    /// Pointer to control node at which the solution is documented
    Node* Control_node_pt;

    //  /// Doc info object
    //  DocInfo Doc_info;

    //  /// Trace file
    //  std::ofstream Trace_file;

    /// Update the problem specs before solve
    // Nothing to do here since no dirichlet boundaries.
    void actions_before_newton_solve(){};

    /// Update the problem specs after solve
    void actions_after_newton_solve(){};

    /// Update the problem specs after solve (empty)
    void actions_after_implicit_timestep() {}

    /// Update the problem specs before next timestep
    void actions_before_implicit_timestep(){};

    /// Set initial condition (incl previous timesteps) according to specified function.
    void set_initial_condition(){};

    /// Create face elements on the b-th boundary
    void create_face_elements(const unsigned& b);

    /// Number of bulk elements
    unsigned N_bulk_element;

  }; // end of problem class

  // /// Set number of values stored at each node
  // template<unsigned DIM, unsigned NNODE_1D>
  //   const unsigned MicromagFaceElement<NNODE_1D>::Initial_Nvalue = 2;

  //=====start_of_constructor===============================================
  ///
  //========================================================================
  template<class ELEMENT>
  TwoDBoundaryProblem<ELEMENT>::
  TwoDBoundaryProblem(const unsigned& n_x, const unsigned& n_y)
  {
    // Set domain size
    double l_x = 1.0;
    double l_y = 1.0;

    // Allocate steady state timestepper
    add_time_stepper_pt(new Steady<2>);

    // Build mesh
    Problem::mesh_pt() = new RectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y,time_stepper_pt());

    // Store number of bulk elements
    N_bulk_element = mesh_pt()->nelement();

    // Choose a control node at which the solution is documented
    unsigned control_el = unsigned(N_bulk_element/2); // Pick a control element in the middle
    Control_node_pt=mesh_pt()->finite_element_pt(control_el)->node_pt(0);  // Choose its first node as the control node
    std::cout << "Recording trace of the solution at: " << Control_node_pt->x(0) << std::endl;

    // Create face elements on all elements adjacent to all boundaries
    unsigned n_boundary = mesh_pt()->nboundary();
    for(unsigned i_boundary=0; i_boundary<n_boundary; i_boundary++)
      create_face_elements(i_boundary);

    //??ds boundary conditions?
    // No Dirichlet conditions so nothing to do here

    // Loop over elements to set pointers
    // No pointers to set so nothing to do here

    // Setup equation numbering scheme
    std::cout << "Number of equations: " << assign_eqn_numbers() << std::endl;

  } // end of constructor

  //============start_of_create_face_elements==============================
  /// Create Face Elements on the b-th boundary of the Mesh.
  //=======================================================================
  template<class ELEMENT>
  void TwoDBoundaryProblem<ELEMENT>::
  create_face_elements(const unsigned &b)
  {
    // How many bulk elements are on the boundary?
    unsigned n_boundary_element = mesh_pt()->nboundary_element(b);

    // Loop over these elements
    for(unsigned e=0; e<n_boundary_element; e++)
      {
	// Get pointer to bulk element
	ELEMENT* bulk_element_pt = dynamic_cast<ELEMENT*>
	  (mesh_pt()->boundary_element_pt(b,e));

	// What is the index of the face of the bulk element at the boundary
	int face_index = mesh_pt()->face_index_at_boundary(b,e);

	// Build the corresponding face element
	MicromagFaceElement<ELEMENT>* face_element_pt =
	  new MicromagFaceElement<ELEMENT>(bulk_element_pt,face_index,mesh_pt());

	//Add the prescribed-flux element to the mesh
	mesh_pt()->add_element_pt(face_element_pt);
      }

  }

  //=====================start_of_doc=======================================
  /// ??ds write this!
  //========================================================================
  template<class ELEMENT>
  void TwoDBoundaryProblem<ELEMENT>::
  doc_solution(DocInfo& doc_info)
  {

    std::ofstream some_file;
    char filename[100];

    // Number of plot points
    unsigned npts;
    npts=5;

    // Output solution
    sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
	    doc_info.number());
    some_file.open(filename);
    mesh_pt()->output(some_file,npts);
    some_file.close();

    // // Output exact solution
    // //----------------------
    // sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
    // 	    doc_info.number());
    // some_file.open(filename);
    // mesh_pt()->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u);
    // some_file.close();

  } // end of doc

} // End of oomph namespace

//==========start_of_main=================================================
/// ??ds
//========================================================================
int main()
{
  // Set number of elements in each direction
  unsigned n_x = 2;
  unsigned n_y = 1;


  //Set up the problem
  TwoDBoundaryProblem<QMicromagElement<2,2> >
    problem(n_x,n_y);

  // Setup doc info
  DocInfo doc_info;
  doc_info.set_directory("RESLT");

  // Check that we're ready to go
  std::cout << "\n\n\nProblem self-test ";
  if (problem.self_test()==0)
    {
      std::cout << "passed: Problem can be solved." << std::endl;
    }
  else
    {
      throw OomphLibError("Self test failed",
			  "main()",
			  OOMPH_EXCEPTION_LOCATION);
    }

  // // Solve the problem
  // problem.newton_solve();

  // //Output solution
  // problem.doc_solution(doc_info);

  // get boundary matrix contribution from first element

  // Eventually loop over all boundary elements?
  unsigned i_ele = problem.get_n_bulk_element();
  {
    // Get element pointer to ith element
    MicromagFaceElement<QMicromagElement<2,2> >* elem_pt =
      dynamic_cast<MicromagFaceElement<QMicromagElement<2,2> >*>
      (problem.mesh_pt()->element_pt(i_ele));

    // Get total number of boundary nodes
    unsigned n_boundary = problem.mesh_pt()->nboundary();
    unsigned n_boundary_nodes(0.0);
    for(unsigned i_boundary=0; i_boundary<n_boundary; i_boundary++)
      n_boundary_nodes += problem.mesh_pt()->nboundary_node(i_boundary);
    std::cout << "There are " << n_boundary_nodes << " nodes on the boundary" << std::endl;

    // Create matrix to store boundary data for this element
    DenseMatrix<double> boundary_matrix(2,n_boundary_nodes,0.0); //??ds not actually 2, jsut for now

    // Get boundary element matrix for this element
    elem_pt->fill_in_elemental_contribution_to_boundary_element_matrix(boundary_matrix);

    // dump for testing:
    std::ofstream matrix_file;
    char filename[100] = "matrix";
    matrix_file.open(filename);
    boundary_matrix.output(matrix_file);
    matrix_file.close();
  }

} //end of main

#endif
