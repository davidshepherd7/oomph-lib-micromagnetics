// Driver code for demag field finite element calculation using charges

// ??ds replace all capital M magnetisations with small m (for compliance with oomph-lib rules - capital means member data)

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
  /// Constructor (must initialise the Source_fct_pt to null), sets flag to not use ALE formulation of equations ??ds change to use ALE once I understand it
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
  unsigned self_test(){return 0;} //??ds write a real test sometime


  // SET PARAMETERS  
  // Specify nodal index of value of phi
  unsigned phi_index_micromag() const {return 0;}

  // Specify nodal index of nth component of M
  unsigned M_index_micromag(const unsigned &n) const {return 1 + n;}

  // Get coefficient of the precession term (M x H) of the LLG equation
  double get_llg_precession_coeff() const
  { 
    double alpha = 0.7;
    double gamma = 0.01;
    return gamma/(alpha*alpha);
  } //??ds temporary  gamma/(1 + alpha^2); }

  // Get coefficient of the damping term (M x (M x H)) of the LLG equation (at position x - could end up not constant if saturisation magnetisation changes)
  double get_llg_damping_coeff(const Vector<double>& x=0) const
  { 
    double alpha = 0.7;
    double gamma = 0.01;
    double M_s = 1.0;
    return (gamma/(alpha*alpha) )* alpha/M_s;
  } //??ds temporary ( gamma/(1+alpha^2) ) * alpha/saturisation_magnetisation;

  void get_applied_field(const Vector<double> &x, Vector<double> &H_applied) const
  {
    H_applied[0] = -20; //??ds temporary - later use a pointer as in poisson source
  }

  /// Turn ALE on/off (needed if mesh is potentially moving - i.e. multiphysics but creates slightly more work)
  void disable_ALE(){ALE_is_disabled=true;}
  void enable_ALE(){ALE_is_disabled=false;}


  // // GET VALUES ??ds don't think I need these functions at all?
  // /// Return the i-th value stored at local node n, DO take hanging nodes into account
  // double nodal_value(const unsigned &n, const unsigned &i) const
  // {return node_pt(n)->value(i);}

  // /// Return the i-th value at time t stored at local node n, DO take hanging nodes into account
  // double nodal_value(const unsigned &t, const unsigned &n, const unsigned &i) const
  // {return node_pt(n)->value(t,i);}
  

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
  void get_demag(const Vector<double>& s, Vector<double>& demag_field) const
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
	demag_field[j] = 0.0;
      }
    
    // Loop over nodes
    for(unsigned l=0;l<n_node;l++) 
      {
	//Loop over derivative directions
	for(unsigned j=0;j<3;j++)
	  {                               
	    demag_field[j] -= this->nodal_value(l,phi_nodal_index)*dpsidx(l,j);
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

  /// Return FE representation of M[j] at local coordinate s
  inline double interpolated_m_micromag(const Vector<double> &s,const unsigned &j) const
  {
    //Find number of nodes
    const unsigned n_node = nnode();
   
    //Local shape function
    Shape psi(n_node);

    //Find values of shape function
    shape(s,psi);
    
    // Initialise m
    double interpolated_m = 0;
	
    //Loop over the local nodes and sum
    for(unsigned l=0;l<n_node;l++) 
      {
	interpolated_m += this->nodal_value(l,M_index_micromag(j))*psi[l];
      }

    return interpolated_m;
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

  // /// C_style output with default number of plot points
  // void output(FILE* file_pt)
  // {
  //   const unsigned n_plot=5;
  //   output(file_pt,n_plot);
  // }

  // /// C-style output FE representation of soln: x,y,u or x,y,z,u at n_plot^DIM plot points
  // void output(FILE* file_pt, const unsigned &n_plot);



  //??ds MISSING COMPUTE ERROR FUNCTIONS - not sure how this can be done with no exact solution



  // RESIDUALS + JACOBIAN
  /// Add the element's contribution to its residual vector (wrapper)
  void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
    //Call the generic residuals function with flag set to 0 using a dummy matrix argument
    fill_in_generic_residual_contribution_micromag(residuals,GeneralisedElement::Dummy_matrix,0);
  }

  // ??ds re-activate this when you're confident that the Jacobian is correct!
  // /// Add the element's contribution to its residual vector and element Jacobian matrix (wrapper)
  // void fill_in_contribution_to_jacobian(Vector<double> &residuals, DenseMatrix<double> &jacobian)
  // {
  //   //Call the generic routine with the flag set to 1
  //   fill_in_generic_residual_contribution_micromag(residuals,jacobian,1);
  // }

  // Fill in contribution to residuals and jacobian (if flag is set) from these equations (compatible with multiphysics)
  void fill_in_generic_residual_contribution_micromag(Vector<double> &residuals, DenseMatrix<double> &jacobian, const unsigned& flag) ;

  /// Get dM/dt at local node n (using timestepper and history values).
  void dM_dt_micromag(const unsigned &n, Vector<double>& dMdt) const
  {
    // Get the data's timestepper
    TimeStepper* time_stepper_pt= this->node_pt(n)->time_stepper_pt();

    //Initialise dMdt to zero
    for (unsigned j=0; j<3; j++) {dMdt[j] = 0.0;}
   
    //Loop over the timesteps, if there is a non Steady timestepper
    if (!time_stepper_pt->is_steady())
      {
	// Loop over the directions
	for (unsigned j=0; j<3; j++)
	  {
	    // Get number of timsteps to use (past & present)
	    const unsigned n_time = time_stepper_pt->ntstorage();
     
	    // Loop over past and present times and add the contributions to the time derivative
	    for(unsigned t=0;t<n_time;t++)
	      {
		// ??ds The one in weights is the derrivative order (i.e. 1: dM/dt, 2: d^2M/dt^2) I think...
		dMdt[j] += time_stepper_pt->weight(1,t)*nodal_value(t,n,M_index_micromag(j));
	      }
	  }

      }

  } // end of dM_dt_micromag

protected:

  // Pointer to source function - only for testing purposes since div(M) is our source function in calculation of the demagnetising potential
  PoissonSourceFctPt Source_fct_pt;

  /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
  virtual double dshape_and_dtest_eulerian_micromag(const Vector<double> &s, Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx) const=0;

  /// Shape, test functions & derivs. w.r.t. to global coords. at integration point ipt. Return Jacobian.
  virtual double dshape_and_dtest_eulerian_at_knot_micromag(const unsigned& ipt, Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx) const=0;

  /// Boolean flag to indicate if ALE formulation is disabled when time-derivatives are computed. Only set to true if you're sure that the mesh is stationary.
  // ??ds not actually implemented and ALE stuff yet...
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

  //Set the value of n_intpt
  const unsigned n_intpt = integral_pt()->nweight();

  //Integers to store the local equation and unknown numbers
  int phi_local_eqn=0, M_local_eqn=0;

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
      double interpolated_phi=0.0,  llg_damping_coeff=0.0, llg_precession_coeff=0.0, div_m=0.0;
      Vector<double> interpolated_x(DIM,0.0), interpolated_dphidx(DIM,0.0);
      Vector<double> H_total(3,0.0), H_demag(3,0.0), H_cryst_anis(3,0.0), H_applied(3,0.0);
      Vector<double> interpolated_m(3,0.0), interpolated_mxH(3,0.0), interpolated_mxmxH(3,0.0);
      Vector<double> dmdt(3,0.0), interpolated_dmdt(3,0.0);
   
      //Calculate function value and derivatives:
      //-----------------------------------------
      // Loop over nodes
      for(unsigned l=0;l<n_node;l++) 
	{
	  //Get the nodal value of phi (the poisson unknown)
	  double phi_value = nodal_value(l,phi_nodal_index);
	  interpolated_phi += phi_value*psi(l);

	  // Get the nodal values of dM/dt
	  dM_dt_micromag(l,dmdt);
	  for(unsigned j=0; j<3; j++)
	    {
	      interpolated_dmdt[j] += dmdt[j]*psi(l);
	      interpolated_m[j] += nodal_value(l,M_index_micromag(j))*psi(l);
	    }

	  // Loop over spatial directions
	  for(unsigned j=0;j<DIM;j++)
	    {
	      interpolated_x[j] += nodal_position(l,j)*psi(l);
	      interpolated_dphidx[j] += phi_value*dpsidx(l,j);
	      div_m += nodal_value(l,M_index_micromag(j))*dpsidx(l,j);
	    }

	} // end of calculating values + derivatives

      // // Get source function  ??ds causes a seg fault - not entirely sure why but don't need it really anyway    
      // double source;
      // get_source_poisson(ipt,interpolated_x,source);

      // Get the demagnetising field (and put in H_total)
      //get_demag(s,H_demag);

      // Get applied field at this position
      get_applied_field(interpolated_x, H_applied);

      // Get crystalline anisotropy effective field ??ds for now just set to [0.5,0,0]
      // use dot products in real implementation
      if( interpolated_m[0] > 0) H_cryst_anis[0] = 0.5;
      else H_cryst_anis[0] = -0.5;
      
      // Take total of all fields used ??ds pass this entire section out to a function eventuall if possible?
      // ??ds add 0.1 to push off maximum (i.e. thermal-ish...)
      for(unsigned j=0; j<3; j++)
	{
	  H_total[j] = H_cryst_anis[j] + H_applied[j] - 1;
	}
            
      // Get the coefficients for the LLG equation (damping could be a function of position if saturation magnetisation varies)
      llg_damping_coeff = get_llg_damping_coeff(interpolated_x);
      llg_precession_coeff = get_llg_precession_coeff();
      
      // Get the cross products for the LLG equation
      cross(interpolated_m, H_total, interpolated_mxH);
      cross(interpolated_m, interpolated_mxH, interpolated_mxmxH);    

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
	      residuals[phi_local_eqn] += (4.0*MathematicalConstants::Pi*div_m)*test(l)*W;

	      // The Poisson bit
	      for(unsigned k=0;k<DIM;k++)
		{
		  residuals[phi_local_eqn] += interpolated_dphidx[k]*dtestdx(l,k)*W;
		}
	      // ??ds add in residuals calculation eventually
	    }

	  // Calculate residuals for the time evolution equations (Landau-Lifschitz-Gilbert):
	  //----------------------------------------
	  // dM/dt = -gamma/(1+alpha^2) [ (M x H) + (gamma/|M_s|)(M x (M x H)) ]

	  // Add the contributions to the residuals from the LLG equation      
	  // loop over M directions
	  for(unsigned k=0; k<3; k++)
	    {
	      // Get the local equation number for the M_x part
	      M_local_eqn = nodal_local_eqn(l,M_index_micromag(k));

	      if(M_local_eqn >= 0)  // If it's not a boundary condition
		{
		  residuals[M_local_eqn] += (interpolated_dmdt[k] + llg_precession_coeff*interpolated_mxH[k] + llg_damping_coeff*interpolated_mxmxH[k])*test(l)*W;
		}

	    }

	}

    }// End of loop over integration points

} // End of fill in residuals function 


//======================================================================
/// Output function:
///
///   x,phi,M_x,M_y,M_z
///
/// nplot points
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
   
      // Output position
      for(unsigned i=0;i<DIM;i++) 
	{
	  outfile << interpolated_x(s,i) << " ";
	}

      // Output phi value at position
      outfile << interpolated_phi_micromag(s) << " ";   

      // Output all M values at position
      for(unsigned i=0; i<3; i++)
	{
	  outfile << interpolated_m_micromag(s,i) << " ";
	}

      // Output div(M) as position
      outfile << divergence_M(s) << " ";
   
      // End the line ready for next point
      outfile << std::endl;
    }

  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(outfile,nplot);

}


// //======================================================================
// /// C-style output function:
// ///
// ///   x,y,u   or    x,y,z,u
// ///
// /// nplot points in each coordinate direction
// //======================================================================
// template <unsigned DIM>
// void  MicromagEquations<DIM>::output(FILE* file_pt,
// 				     const unsigned &nplot)
// {
//   //Vector of local coordinates
//   Vector<double> s(DIM);
 
//   // Tecplot header info
//   fprintf(file_pt,"%s",tecplot_zone_string(nplot).c_str());

//   // Loop over plot points
//   unsigned num_plot_points=nplot_points(nplot);
//   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
//     {
//       // Get local coordinates of plot point
//       get_s_plot(iplot,nplot,s);
   
//       for(unsigned i=0;i<DIM;i++) 
// 	{
// 	  fprintf(file_pt,"%g ",interpolated_x(s,i));
// 	}
//       fprintf(file_pt,"%g \n",interpolated_phi_micromag(s));
//     }

//   // Write tecplot footer (e.g. FE connectivity lists)
//   write_tecplot_zone_footer(file_pt,nplot);
// }

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


  // /// C-style output function: x,y,u or x,y,z,u
  // void output(FILE* file_pt)
  // {MicromagEquations<DIM>::output(file_pt);}


  // /// C-style output function: x,y,u or x,y,z,u at n_plot^DIM plot points
  // void output(FILE* file_pt, const unsigned &n_plot)
  // {MicromagEquations<DIM>::output(file_pt,n_plot);}


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
    // Start with step change in M_x from 1 to -1 at x = 0.5
    // Should create a domain wall type structure
    M[0] = 1; //tanh(5*(x[0] - 0.5));

    // Initialise y and z components of M to zero
    M[1] = 0;
    M[2] = 0;



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

  /// Pointer to control node at which the solution is documented ??ds - not sure what this is
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
  //??ds this will break if number of nodal data values changes in different elements, don't think it does change though
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

  // Get pointer to an element (any will do so take 0th)
  ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0)); 

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
	      mesh_pt()->node_pt(n)->set_value(t,elem_pt->M_index_micromag(i),M[i]);
	    }
     
	  // Loop over coordinate directions: Mesh doesn't move, so previous position = present position
	  // ??ds presumably this is where the ALE formulation would/will/should come in
	  for (unsigned i=0;i<1;i++)
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

  // // Write file as a tecplot text object
  // some_file << "TEXT X=2.5,Y=93.6,F=HELV,HU=POINT,C=BLUE,H=26,T=\"time = " 
  // 	    << time_pt()->time() << "\"";
  // // ...and draw a horizontal line whose length is proportional
  // // to the elapsed time
  // some_file << "GEOMETRY X=2.5,Y=98,T=LINE,C=BLUE,LT=0.4" << std::endl;
  // some_file << "1" << std::endl;
  // some_file << "2" << std::endl;
  // some_file << " 0 0" << std::endl;
  // some_file << time_pt()->time()*20.0 << " 0" << std::endl;

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
  unsigned n_element=40; //Number of elements
  OneDMicromagProblem<QMicromagElement<1,2> > problem(n_element);


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
  double t_max=30;
  double dt=0.1;

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
      cout << "Timestep " << istep << std::endl;
   
      // Take timestep
      problem.unsteady_newton_solve(dt);
   
      //Output solution
      problem.doc_solution(doc_info,trace_file);
   
      //Increment counter for solutions 
      doc_info.number()++;
    }
 
  // Close trace file
  trace_file.close();

  return 0;
  
} // end of main
