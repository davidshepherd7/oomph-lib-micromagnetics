// Driver code for demag field finite element calculation using charges

// Generic oomph-lib routines
#include "generic.h"

// Include Poisson elements/equations
//#include "poisson.h"

// Include mesh
#include "meshes/one_d_mesh.h"

using namespace std;

using namespace oomph;

//==================================================
/// A class for the maths used in solving the Landau-Lifshitz-Gilbert equations.
//==================================================  
template<unsigned DIM>
class MicromagEquations : public virtual FiniteElement
{
private:

public:

  // CONSTRUCTORS ETC.
  /// Constructor (must initialise the Source_fct_pt to null), sets flag to not use ALE formulation of equations ?? change to use ALE once I understand it
  MicromagEquations() : Source_fct_pt(0), ALE_is_disabled(true) {}
 
  /// Broken copy constructor
  MicromagEquations(const MicromagEquations& dummy) 
  { 
    BrokenCopy::broken_copy("MicromagEquations");
  } 
 
  /// Broken assignment operator
  void operator=(const MicromagEquations&) 
  {
    BrokenCopy::broken_assign("MicromagEquations");
  }

  /// Self-test: Return 0 for OK
  unsigned self_test(){return 0;} //?? write a real test sometime


  // SET PARAMETERS
  // Specify the number of nodal values required at each node.
  // unsigned required_nvalue(const unsigned &n) const
  // {
  //   // Remember to change this when adding new fields.
  //   return 4;
  // } // commented out because defined in Qmicromagelement
  
  // Specify nodal index of value of phi
  unsigned phi_index_micromag() const
  {
    return 0;
  }

  // Specify nodal index of nth component of M
  unsigned M_index_micromag(const unsigned &n) const
  {
    return 1 + n;
  }

  // Get coefficients used in the LLG equation (at node n - could end up not constant)
  double get_LLG_damping_coeff(const unsigned &n=0) const
  { return 1;} //?? temporary ( gamma/(1+alpha^2) ) * alpha/saturisation_magnetisation;
 
  double get_LLG_precession_coeff(const unsigned &n=0) const
  { return 1;} //?? temporary  gamma/(1 + alpha^2); }


  /// Turn ALE on/off (needed if mesh is potentially moving - i.e. multiphysics but creates slightly more work)
  void disable_ALE(){ALE_is_disabled=true;}
  void enable_ALE(){ALE_is_disabled=false;}

  // GET VALUES
  /// Return the i-th value stored at local node n but do NOT take hanging nodes into account
  double raw_nodal_value(const unsigned &n, const unsigned &i)
  {return node_pt(n)->raw_value(i);}
  

  // SOURCE FUNCTION
  /// \short Function pointer to source function fct(x,f(x)) -- 
  /// x is a Vector! 
  typedef void (*PoissonSourceFctPt)(const Vector<double>& x, double& f);

  /// Access function: Pointer to source function
  PoissonSourceFctPt& source_fct_pt() {return Source_fct_pt;}

  /// Access function: Pointer to source function. Const version
  PoissonSourceFctPt source_fct_pt() const {return Source_fct_pt;}

  /// Get source term at (Eulerian) position x. This function is
  /// virtual to allow overloading in multi-physics problems where
  /// the strength of the source function might be determined by
  /// another system of equations.
  inline virtual void get_source_poisson(const unsigned& ipt,
					 const Vector<double>& x,
					 double& source) const
  {
    //If no source function has been set, return zero
    if(Source_fct_pt==0) {source = 0.0;}
    else
      {
	// Get source strength
	(*Source_fct_pt)(x,source);
      }
  }


  // OPERATIONS ON RESULTS
  /// Get demagnetising field = - grad(phi(s)) = - phi(s) * d_psi/dx (s)
  void get_demag(const Vector<double>& s, Vector<double>& flux) const
  {
    //Find out how many nodes there are in the element
    const unsigned n_node = nnode();

    //Get the index at which the unknown is stored
    const unsigned phi_nodal_index = phi_index_micromag();

    //Set up memory for the shape and test functions
    Shape psi(n_node);
    DShape dpsidx(n_node,3);
 
    //Call the derivatives of the shape and test functions
    dshape_eulerian(s,psi,dpsidx);
     
    //Initialise to zero
    for(unsigned j=0;j<3;j++)
      {
	flux[j] = 0.0;
      }
   
    // Loop over nodes
    for(unsigned l=0;l<n_node;l++) 
      {
	//Loop over derivative directions
	for(unsigned j=0;j<3;j++)
	  {                               
	    flux[j] += -1*this->nodal_value(l,phi_nodal_index)*dpsidx(l,j);
	  }
      }
  }

  /// Get divergence: divergence(M[i]) = dM_i/dx_i at local coordinate s
  double divergence_M(const Vector<double>& s) const
  {

    //Find out how many nodes there are in the element
    const unsigned n_node = nnode();

    // Create variable for div_M
    double div_M = 0;

    //Set up memory for the shape and test functions
    Shape psi(n_node);
    DShape dpsidx(n_node,DIM);
 
    //Call the derivatives of the shape and test functions
    dshape_eulerian(s,psi,dpsidx);
   
    // Loop over nodes (take sum over nodes)
    for(unsigned l=0;l<n_node;l++) 
      {
	// Loop over M directions and derivative directions (take sum over directions
	// ?? unsure how to handle dimensions here - M is always 3D but spatial variations only occur over DIM dimensions. For now will only calculate for DIM dimensions
	for(unsigned j=0;j<DIM;j++)
	  {                               
	    div_M += this->nodal_value(l,M_index_micromag(j))*dpsidx(l,j);
	  }
      }

    return div_M;
  }

  // Calculate the cross product of vectors A and B, store the result in vector output. NOTE: the cross product is only valid for 3-dimensional vectors
  void cross(Vector<double>& A, Vector<double>& B, Vector<double>& output) const
  {
    //?? #IFDEF PARANOID - check that vectors have length 3
    output[0] = A[1]*B[2] - A[2]*B[1];
    output[1] = A[2]*B[0] - A[0]*B[2];
    output[2] = A[0]*B[1] - A[1]*B[0];
  }
  
  /// Return FE representation of function value phi(s) at local coordinate s
  inline double interpolated_phi_micromag(const Vector<double> &s) const
  {
    //Find number of nodes
    const unsigned n_node = nnode();

    //Get the index at which the poisson unknown is stored
    const unsigned phi_nodal_index = phi_index_micromag();
   
    //Local shape function
    Shape psi(n_node);

    //Find values of shape function
    shape(s,psi);

    //Initialise value of u
    double interpolated_phi = 0.0;

    //Loop over the local nodes and sum
    for(unsigned l=0;l<n_node;l++) 
      {
	interpolated_phi += this->nodal_value(l,phi_nodal_index)*psi[l];
      }

    return(interpolated_phi);
  }


  // OUTPUT FUNCTIONS
  /// Output with default number of plot points
  void output(std::ostream &outfile) 
  {
    const unsigned n_plot=5;
    output(outfile,n_plot);
  }

  /// Output FE representation of soln: x,y,u or x,y,z,u at n_plot^DIM plot points
  void output(std::ostream &outfile, const unsigned &n_plot);

  /// C_style output with default number of plot points
  void output(FILE* file_pt)
  {
    const unsigned n_plot=5;
    output(file_pt,n_plot);
  }

  /// C-style output FE representation of soln: x,y,u or x,y,z,u at n_plot^DIM plot points
  void output(FILE* file_pt, const unsigned &n_plot);



  //?? MISSING COMPUTE ERROR FUNCTIONS - not sure how this can be done with no exact solution



  // RESIDUALS + JACOBIAN
  /// Add the element's contribution to its residual vector (wrapper)
  void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
    //Call the generic residuals function with flag set to 0 using a dummy matrix argument
    fill_in_generic_residual_contribution_micromag(residuals,GeneralisedElement::Dummy_matrix,0);
  }

  /// Add the element's contribution to its residual vector and element Jacobian matrix (wrapper)
  void fill_in_contribution_to_jacobian(Vector<double> &residuals, DenseMatrix<double> &jacobian)
  {
    //Call the generic routine with the flag set to 1
    fill_in_generic_residual_contribution_micromag(residuals,jacobian,1);
  }

  // Fill in contribution to residuals and jacobian (if flag is set) from these equations (compatible with multiphysics)
  void fill_in_generic_residual_contribution_micromag(Vector<double> &residuals, DenseMatrix<double> &jacobian, const unsigned& flag) ;

  /// Compute derivatives of elemental residual vector with respect to nodal coordinates. Overwrites default implementation in FiniteElement base class. dresidual_dnodal_coordinates(l,i,j) = d res(l) / dX_{ij}
  //  virtual void get_dresidual_dnodal_coordinates(RankThreeTensor<double>&dresidual_dnodal_coordinates);

protected:

  // Pointer to source function - only for testing purposes since div(M) is our source function in calculation of the demagnetising potential
  PoissonSourceFctPt Source_fct_pt;

  /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
  //  virtual double dshape_and_dtest_eulerian_micromag(const Vector<double> &s, Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx) const=0;


  /// Shape, test functions & derivs. w.r.t. to global coords. at integration point ipt. Return Jacobian.
  virtual double dshape_and_dtest_eulerian_at_knot_micromag(const unsigned& ipt, Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx) const=0;

  /// Boolean flag to indicate if ALE formulation is disabled when time-derivatives are computed. Only set to true if you're sure that the mesh is stationary.
  // ?? not actually implemented and ALE stuff yet...
  bool ALE_is_disabled;


}; // End of MicromagElements class






//======================================================================
/// Compute element residual Vector and/or element Jacobian matrix 
/// 
/// flag=1: compute both
/// flag=0: compute only residual Vector
///
/// Pure version without hanging nodes
//======================================================================
template<unsigned DIM>
void MicromagEquations<DIM>::fill_in_generic_residual_contribution_micromag(Vector<double> &residuals, DenseMatrix<double> &jacobian, const unsigned& flag) 
{

  //Find out how many nodes there are
  const unsigned n_node = nnode();

  //Set up memory for the shape and test functions
  Shape psi(n_node), test(n_node);
  DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);

  // Set up vector to store the local coordinates
  Vector<double> s(DIM);

  //Index at which the poisson unknown is stored
  const unsigned phi_nodal_index = phi_index_micromag();

  // Index at which the first direction of the magnetisation is stored
  const unsigned M_nodal_index = M_index_micromag(0);
 
  //Set the value of n_intpt
  const unsigned n_intpt = integral_pt()->nweight();

  //Integers to store the local equation and unknown numbers
  int phi_local_eqn=0, M_local_eqn=0, local_unknown=0;

  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
      //Get the integral weight
      double w = integral_pt()->weight(ipt);

      //Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_micromag(ipt,psi,dpsidx,test,dtestdx);
       
      //Premultiply the weights and the Jacobian
      double W = w*J;

      // Get values of s (local coordinate)
      for(unsigned j=0; j<DIM; j++) {s[j] = integral_pt()->knot(ipt,j);}

      //Allocate memory for local quantities and initialise to zero
      double interpolated_phi=0.0,  LLG_damping_coeff=0.0, LLG_precession_coeff=0.0;
      Vector<double> interpolated_x(DIM,0.0);
      Vector<double> interpolated_dphidx(DIM,0.0);
      Vector<double> dMdt(3,0.0); // M must be 3D for LLG equation to make sense
      Vector<double> H_eff(3,0.0), M_cross_H(3,0.0), M_cross_M_cross_H(3,0.0);
   
      //Calculate function value and derivatives:
      //-----------------------------------------
      // Loop over nodes
      for(unsigned l=0;l<n_node;l++) 
	{
	  //Get the nodal value of the poisson unknown (phi)
	  double phi_value = raw_nodal_value(l,phi_nodal_index);
	  interpolated_phi += phi_value*psi(l);
	  dMdt += dM_dt_micromag(l)*psi(l);

	  // Loop over spatial directions
	  for(unsigned j=0;j<DIM;j++)
	    {
	      interpolated_x[j] += raw_nodal_position(l,j)*psi(l);
	      interpolated_dphidx[j] += phi_value*dpsidx(l,j);
	    }

	} // end of calculating function values + derivatives

      // Get mesh velocity if needed
      if (!ALE_is_disabled)
	{
	  // ?? ALE section
	}

      // Get source function      
      double source;
      get_source_poisson(ipt,interpolated_x,source);

      //Get divergence of M
      double div_M  =  divergence_M(s);

      // Assemble residuals and Jacobian
      //--------------------------------
       
      // Loop over the test functions
      for(unsigned l=0;l<n_node;l++)
	{
	  
	  // Calculate residual for poisson equation (phi):
	  //----------------------------------------

	  // Get the local equation number for the poisson part
	  phi_local_eqn = nodal_local_eqn(l,phi_nodal_index);

	  if(phi_local_eqn >= 0)	  // If it's not a boundary condition:
	    {
	      // Add source term and 4*pi*divergence(M) 
	      residuals[phi_local_eqn] += (source + 4*MathematicalConstants::Pi*div_M)*test(l)*W;

	      // The Poisson bit
	      for(unsigned k=0;k<DIM;k++)
		{
		  residuals[phi_local_eqn] += interpolated_dphidx[k]*dtestdx(l,k)*W;
		}

	      // Calculate the jacobian
	      if(flag)
		{
		  //Loop over the velocity shape functions again
		  for(unsigned l2=0;l2<n_node;l2++)
		    { 
		      local_unknown = nodal_local_eqn(l2,phi_nodal_index);
		      //If at a non-zero degree of freedom add in the entry
		      if(local_unknown >= 0)
			{
			  //Add contribution to Elemental Matrix
			  for(unsigned i=0;i<DIM;i++)
			    {
			      jacobian(phi_local_eqn,local_unknown) 
				+= dpsidx(l2,i)*dtestdx(l,i)*W;
			    }
			}
		    }
		}
	    } // end of phi calculations



	  // Calculate residuals for the time evolution equations (Landau-Lifschitz-Gilbert):
	  //----------------------------------------
	  // dM/dt = -gamma/(1+alpha^2) [ (M x H) + (gamma/|M_s|)(M x (M x H)) ]

	  // Get the local equation number for the M_x part
	  M_local_eqn = nodal_local_eqn(l,M_nodal_index);
	  if(M_local_eqn >= 0)  // If it's not a boundary condition
	    {
	      // Get the demagnetising field (and put in H_eff)
	      get_demag(s,H_eff);
	      
	      // Get the coefficients for the LLG equation
	      LLG_damping_coeff = get_LLG_damping_coeff(l);
	      LLG_precession_coeff = get_LLG_precession_coeff(l);

	      // Get the cross products for the LLG equation
	      cross(M,H_eff,M_cross_H);
	      cross(M,M_cross_H,M_cross_M_cross_H);
	      
	      // Add the contributions to the residuals from the LLG equation
	      for(unsigned k=0; k<3; k++)
		{
		  residuals[M_local_eqn + k] += (dMdt[k] + LLG_precession_coeff*M_cross_H[k] + LLG_damping_coeff*M_cross_cross_H[k])*test(l)*W;
		}

	      // Calculate the jacobian
	      if(flag)
		{
		  //Loop over the velocity shape functions again
		  for(unsigned l2=0;l2<n_node;l2++)
		    { 
		      local_unknown = nodal_local_eqn(l2,phi_nodal_index);
		      //If it's not a boundary condition:
		      if(local_unknown >= 0)
		  	{
		  	  //Add contribution to Elemental Matrix
		  	  // for(unsigned i=0;i<DIM;i++)
		  	  //   {
		  	  //     jacobian(local_eqn,local_unknown) 
		  	  // 	+= dpsidx(l2,i)*dtestdx(l,i)*W;
		  	  //   }
		  	}
		  }
		}
	    }
	}

    } // End of loop over integration points

} // End of fill in residuals function 


//======================================================================
/// Output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
template <unsigned DIM>
void  MicromagEquations<DIM>::output(std::ostream &outfile, 
				     const unsigned &nplot)
{

  //Vector of local coordinates
  Vector<double> s(DIM);
 
  // Tecplot header info
  outfile << tecplot_zone_string(nplot);
 
  // Loop over plot points
  unsigned num_plot_points=nplot_points(nplot);
  for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
   
      // Get local coordinates of plot point
      get_s_plot(iplot,nplot,s);
   
      for(unsigned i=0;i<DIM;i++) 
	{
	  outfile << interpolated_x(s,i) << " ";
	}
      outfile << interpolated_phi_micromag(s) << std::endl;   
   
    }

  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(outfile,nplot);

}


//======================================================================
/// C-style output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
template <unsigned DIM>
void  MicromagEquations<DIM>::output(FILE* file_pt,
				     const unsigned &nplot)
{
  //Vector of local coordinates
  Vector<double> s(DIM);
 
  // Tecplot header info
  fprintf(file_pt,"%s",tecplot_zone_string(nplot).c_str());

  // Loop over plot points
  unsigned num_plot_points=nplot_points(nplot);
  for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot,nplot,s);
   
      for(unsigned i=0;i<DIM;i++) 
	{
	  fprintf(file_pt,"%g ",interpolated_x(s,i));
	}
      fprintf(file_pt,"%g \n",interpolated_phi_micromag(s));
    }

  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(file_pt,nplot);
}

//=================================================================
/// A class combining the micromag equations with a QElement geometry
//=================================================================
template < unsigned DIM, unsigned NNODE_1D>
class QMicromagElement : public QElement<DIM,NNODE_1D>, public MicromagEquations<DIM>
{

private:

  /// Static int that holds the number of variables at nodes: always the same
  static const unsigned Initial_Nvalue;

public:

  /// Constructor: Call constructors for QElement and Micromag equations
  QMicromagElement() : QElement<DIM,NNODE_1D>(), MicromagEquations<DIM>()
  {}
 
  /// Broken copy constructor
  QMicromagElement(const QMicromagElement<DIM,NNODE_1D>& dummy) 
  { 
    BrokenCopy::broken_copy("QMicromagElement");
  } 
 
  /// Broken assignment operator
  void operator=(const QMicromagElement<DIM,NNODE_1D>&) 
  {
    BrokenCopy::broken_assign("QMicromagElement");
  }

  /// Required  # of `values' (pinned or dofs) at node n.
  inline unsigned required_nvalue(const unsigned &n) const 
  {return Initial_Nvalue;}

  // OUTPUT FUNCTIONS (just call from MicromagEquations class)
  /// Output function: x,y,u or x,y,z,u
  void output(std::ostream &outfile)
  {MicromagEquations<DIM>::output(outfile);}


  /// Output function: x,y,u or x,y,z,u at n_plot^DIM plot points
  void output(std::ostream &outfile, const unsigned &n_plot) 
  {MicromagEquations<DIM>::output(outfile,n_plot);}


  /// C-style output function: x,y,u or x,y,z,u
  void output(FILE* file_pt)
  {MicromagEquations<DIM>::output(file_pt);}


  /// C-style output function: x,y,u or x,y,z,u at n_plot^DIM plot points
  void output(FILE* file_pt, const unsigned &n_plot)
  {MicromagEquations<DIM>::output(file_pt,n_plot);}


protected:

  /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
  inline double dshape_and_dtest_eulerian_micromag(const Vector<double> &s, Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx) const;

  /// Shape, test functions & derivs. w.r.t. to global coords. at integration point ipt. Return Jacobian.
  inline double dshape_and_dtest_eulerian_at_knot_micromag(const unsigned& ipt, Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx) const;

}; // end of QMicromagElement class declaration

//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
// Might be useful later so keep (but not used now)
template<unsigned DIM, unsigned NNODE_1D>
double QMicromagElement<DIM,NNODE_1D>::dshape_and_dtest_eulerian_micromag(const Vector<double> &s, Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx) const
{
  //Call the geometrical shape functions and derivatives  
  const double J = this->dshape_eulerian(s,psi,dpsidx);

  //Set the test functions equal to the shape functions
  test = psi;
  dtestdx= dpsidx;
 
  //Return the jacobian
  return J;
}


//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned DIM, unsigned NNODE_1D>
double QMicromagElement<DIM,NNODE_1D>::dshape_and_dtest_eulerian_at_knot_micromag(const unsigned &ipt, Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx) const
{
  //Call the geometrical shape functions and derivatives
  const double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);

  //Set the pointers of the test functions
  test = psi;
  dtestdx = dpsidx;

  //Return the jacobian
  return J;
}


//==start_of_namespace================================================
/// Namespace for calculation of scalar potential
//====================================================================
namespace OneDMicromagSetup
{
  // void source_function(const Vector<double>& x, double& source)
  // {
  //   // Initially try with M_x = x so source fn = -div(M) == -1
  //   // source = -1;

  //   // Try with M_x  like:  _/\_ (i.e. up to one then back down), then source fn is step fn with step at x = 0.5, source = -1 at x = 0, source = 0 at x = 1.
  //   // source = -1*tanh(500*(x[0] - 0.5));

  //   // unfortunately source_function is fn of x but divergence is function of s (local coord)
  //   source = MicromagEquations<DIM>::divergence_M(s);
  // }

  void get_initial_M(const Vector<double>& x, Vector<double>& M)
  {
    // Initialise x and y components of M to zero
    M[0] = 0;
    M[1] = 0;

    // Start with step change in M_z from 1 to -1 at x = 0.5
    // Should create a domain wall type structure
    M[2] = tanh(500*(x[0] - 0.5));

  }


  double get_boundary_phi(const Vector<double>& x)
  {
    // Set all boundaries to zero for now
    return 0;
  }
};



//==start_of_problem_class============================================
/// 1D Micromag problem in unit interval.
//====================================================================
template<class ELEMENT> 
class OneDMicromagProblem : public Problem
{

private:

  /// Pointer to source function
  MicromagEquations<1>::PoissonSourceFctPt Source_fct_pt;

  /// Pointer to control node at which the solution is documented ?? - not sure what this is
  Node* Control_node_pt;

  // Doc info object
  DocInfo Doc_info;
  
  // Trace file
  ofstream Trace_file;

public:

  /// Constructor: Pass number of elements and pointer to source function
  OneDMicromagProblem(const unsigned& n_element);

  /// Destructor (empty -- all the cleanup is done in the base class)
  ~OneDMicromagProblem(){};

  /// Update the problem specs before solve: Set boundary conditions
  void actions_before_newton_solve(){};

  /// Update the problem specs after solve (calculate demag field)
  void actions_after_newton_solve(const unsigned n_element);

  /// Doc the solution
  void doc_solution(DocInfo& doc_info, std::ofstream& trace_file);

  // TIME STEPPING FUNCTIONS
  /// Update the problem specs after solve (empty)
  void actions_after_implicit_timestep() {}

  /// Update the problem specs before next timestep
  void actions_before_implicit_timestep();

  /// Set initial condition (incl previous timesteps) according to specified function. 
  void set_initial_condition();

}; // end of problem class

template<unsigned DIM, unsigned NNODE_1D>
const unsigned QMicromagElement<DIM,NNODE_1D>::Initial_Nvalue = 4;



//=====start_of_constructor===============================================
/// \short Constructor for 1D Micromag problem in unit interval.
/// Discretise the 1D domain with n_element elements of type ELEMENT.
/// Specify function pointer to source function. 
//========================================================================
template<class ELEMENT>
OneDMicromagProblem<ELEMENT>::OneDMicromagProblem(const unsigned& n_element)
{  
  // Allocate the timestepper -- this constructs the Problem's time object with a sufficient amount of storage to store the previous timsteps. 
  add_time_stepper_pt(new BDF<2>);

  // Set domain length 
  double L=1.0;

  // Build mesh and store pointer in Problem
  Problem::mesh_pt() = new OneDMesh<ELEMENT>(n_element,L,time_stepper_pt());

  // Choose a control node at which the solution is documented
  unsigned control_el = unsigned(n_element/2); // Pick a control element in the middle
  Control_node_pt=mesh_pt()->finite_element_pt(control_el)->node_pt(0);  // Choose its first node as the control node
  cout << "Recording trace of the solution at: " << Control_node_pt->x(0) << std::endl;


  // Set up the boundary conditions for this problem: pin the nodes at either end
  mesh_pt()->boundary_node_pt(0,0)->pin(0);
  mesh_pt()->boundary_node_pt(1,0)->pin(0);


  // Loop over elements to set pointers to source function and time
  for(unsigned i=0;i<n_element;i++)
    {
      // Upcast from GeneralisedElement to the present element
      ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
      //Set the source function pointer
      elem_pt->source_fct_pt() = Source_fct_pt;

      // Set pointer to continous time
      elem_pt->time_pt() = time_pt();
    }

 
  // Setup equation numbering scheme
  assign_eqn_numbers();

} // end of constructor



//=========start of actions_before_implicit_timestep===============================
/// \short Actions before timestep: update the domain, then reset the 
/// boundary conditions for the current time.
//========================================================================
template<class ELEMENT>
void OneDMicromagProblem<ELEMENT>::actions_before_implicit_timestep()
{

  // Get pointer to (0th) element - exact element doesn't matter for this use, hopefully!
  //?? this will break if number of nodal data values changes in different elements, don't think it does change though
  ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0));

  // Get index at which phi is stored (in nodal data)
  unsigned phi_nodal_index = elem_pt->phi_index_micromag();

  // Set boundary conditions (more general method that works in higher dimensions)
  // Loop over all boundaries
  unsigned num_bound = mesh_pt()->nboundary();
  for(unsigned ibound=0;ibound<num_bound;ibound++)
    {
      // Loop over the nodes on this boundary
      unsigned num_nod=mesh_pt()->nboundary_node(ibound);
      for (unsigned inod=0;inod<num_nod;inod++)
	{
	  // Set boundary conditions at this node
	  Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
	  Vector<double> x(1,nod_pt->x(0));
	  double phi_value = OneDMicromagSetup::get_boundary_phi(x);
	  nod_pt->set_value(phi_nodal_index,phi_value);
	}
    }

}


//======================start_of_set_initial_condition====================
/// \short Set initial condition: Assign previous and current values
/// from exact solution.
//========================================================================
template<class ELEMENT>
void OneDMicromagProblem<ELEMENT>::set_initial_condition()
{ 
  // Backup time in global Time object
  double backed_up_time=time_pt()->time();
         
  // Past history needs to be established for t=time0-deltat, ...
  // Then provide current values (at t=time0) which will also form
  // the initial guess for the first solve at t=time0+deltat
 
  // Vector of exact solution value
  Vector<double> M(3);
  Vector<double> x(1);

  //Find number of nodes in mesh
  unsigned num_nod = mesh_pt()->nnode();

  // Get index at which M first component is stored (in node data)
  // ?? this wil break if M index changes in difference elements - don't think it does
  ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0)); // First get pointer to (0th) element
  unsigned M_nodal_index = elem_pt->M_index_micromag(0); // then use elem_pt to call the function to find the index

  // Set continuous times at previous timesteps:
  // How many previous timesteps does the timestepper use?
  int nprev_steps=time_stepper_pt()->nprev_values();
  Vector<double> prev_time(nprev_steps+1);
  for (int t=nprev_steps;t>=0;t--)
    {
      prev_time[t]=time_pt()->time(unsigned(t));
    } 

  // Loop over current & previous timesteps
  for (int t=nprev_steps;t>=0;t--)
    {
      // Continuous time
      double time = prev_time[t];
      cout << "setting IC at time =" << time << std::endl;
   
      // Loop over the nodes to set initial values everywhere
      for (unsigned n=0;n<num_nod;n++)
	{
	  // Get nodal coordinate
	  x[0]=mesh_pt()->node_pt(n)->x(0);

	  // Get initial value of M
	  OneDMicromagSetup::get_initial_M(x,M);
     
	  // Assign solution
	  for(unsigned i=0; i<3; i++)
	    {
	      // Set ith direction of M on node n at time t to be M[i]
	      mesh_pt()->node_pt(n)->set_value(t,M_nodal_index+i,M[i]);
	    }
     
	  // Loop over coordinate directions: Mesh doesn't move, so previous position = present position
	  // ?? presumably this is where the ALE formulation would/will/should come in
	  for (unsigned i=0;i<2;i++)
	    {
	      mesh_pt()->node_pt(n)->x(t,i)=x[i];
	    }
	} 
    }

  // Reset backed up time for global timestepper
  time_pt()->time()=backed_up_time;

} // end of set_initial_condition



//===start_of_doc=========================================================
/// Doc the solution in tecplot format. Label files with label.
//========================================================================
template<class ELEMENT>
void OneDMicromagProblem<ELEMENT>::doc_solution(DocInfo& doc_info, std::ofstream& trace_file)
{ 

  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts;
  npts=5; 

  cout << std::endl;
  cout << "=================================================" << std::endl;
  cout << "Doc'ing solution for t=" << time_pt()->time() << std::endl;
  cout << "=================================================" << std::endl;

  // Output solution 
  //-----------------
  sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
	  doc_info.number());
  some_file.open(filename);
  mesh_pt()->output(some_file,npts);

  // Write file as a tecplot text object
  some_file << "TEXT X=2.5,Y=93.6,F=HELV,HU=POINT,C=BLUE,H=26,T=\"time = " 
	    << time_pt()->time() << "\"";
  // ...and draw a horizontal line whose length is proportional
  // to the elapsed time
  some_file << "GEOMETRY X=2.5,Y=98,T=LINE,C=BLUE,LT=0.4" << std::endl;
  some_file << "1" << std::endl;
  some_file << "2" << std::endl;
  some_file << " 0 0" << std::endl;
  some_file << time_pt()->time()*20.0 << " 0" << std::endl;
  some_file.close();
} // end of doc

 

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////



  //======start_of_main==================================================
  /// Driver for 1D Micromag problem
  //=====================================================================
int main()
{

  // Set up the problem:
  unsigned n_element=20; //Number of elements
  OneDMicromagProblem<QMicromagElement<1,2> >
    problem(n_element);


  // SET UP OUTPUT
  // Setup labels for output
  DocInfo doc_info;

  // Output directory
  doc_info.set_directory("RESLT");
 
  // Output number
  doc_info.number()=0;
 
  // Open a trace file
  ofstream trace_file;
  char filename[100];   
  sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
  trace_file.open(filename);
  trace_file << "VARIABLES=\"time\",\"u<SUB>FE</SUB>\","
	     << "\"u<SUB>exact</SUB>\",\"norm of error\",\"norm of solution\""
	     << std::endl;
  

  // SET UP TIME STEPPING
  // Choose simulation interval and timestep
  double t_max=0.5;
  double dt=0.01;

  // Initialise timestep -- also sets the weights for all timesteppers
  // in the problem.
  problem.initialise_dt(dt);
 
  // Set initial condition (on M)
  problem.set_initial_condition();

  // Output initial conditions
  problem.doc_solution(doc_info,trace_file);   //Output initial condition
  doc_info.number()++; //Increment counter for solutions 

  // Check whether the problem can be solved
  cout << "\n\n\nProblem self-test ";
  if (problem.self_test()==0)  
    {
      cout << "passed: Problem can be solved." << std::endl;
    }
  else 
    {
      throw OomphLibError("failed!",
			  "main()",
			  OOMPH_EXCEPTION_LOCATION);
    }

  // SOLVE THE PROBLEM
  // Find number of steps
  unsigned nstep = unsigned(t_max/dt);

  // Timestepping loop
  for (unsigned istep=0;istep<nstep;istep++)
    {
      cout << " Timestep " << istep << std::endl;
   
      // Take timestep
      problem.unsteady_newton_solve(dt);
   
      //Output solution
      problem.doc_solution(doc_info,trace_file);
   
      //Increment counter for solutions 
      doc_info.number()++;
    }
 
  // Close trace file
  trace_file.close();
  
} // end of main
