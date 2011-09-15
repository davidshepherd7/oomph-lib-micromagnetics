// An element for micromagnetics calculations (also including a geometric Qelement version).

// ??ds replace all capital M magnetisations with small m (for compliance with oomph-lib rules - capital means member data)

// Generic oomph-lib routines
#include "generic.h"

// Include mesh
#include "meshes/one_d_mesh.h"

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
  /// Constructor (must initialise the various function pointers to null), sets flag to not use ALE formulation of equations ??ds change to use ALE once I understand it
  MicromagEquations() : Source_fct_pt(0), Applied_field_fct_pt(0), Cryst_anis_field_fct_pt(0), ALE_is_disabled(true) {}
 
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
    double gamma = 0.1;
    return gamma/(alpha*alpha);
  } //??ds temporary  gamma/(1 + alpha^2); }

  // Get coefficient of the damping term (M x (M x H)) of the LLG equation (at position x - could end up not constant if saturisation magnetisation changes)
  double get_llg_damping_coeff(const Vector<double>& x=0) const
  { 
    double alpha = 0.7;
    double gamma = 0.1;
    double M_s = 1.0;
    return (gamma/(alpha*alpha) )* alpha/M_s;
  } //??ds temporary ( gamma/(1+alpha^2) ) * alpha/saturisation_magnetisation;

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

  // APPLIED FIELD
  /// Function pointer to applied field function
  typedef void (*AppliedFieldFctPt)(const Vector<double>& x,Vector<double>& H_applied);

  /// Access function: Pointer to applied field function
  AppliedFieldFctPt& applied_field_fct_pt() {return Applied_field_fct_pt;}

  /// Access function: Pointer to applied field function. Const version
  AppliedFieldFctPt applied_field_fct_pt() const {return Applied_field_fct_pt;}

  /// Get the applied field at Eulerian position x.
  inline virtual void get_applied_field(const Vector<double> &x, Vector<double> &H_applied) const
  {
    //If no applied field has been set, return zero vector
    if(Applied_field_fct_pt==0) {for(unsigned j=0;j<3;j++) H_applied[j] = 0.0;}
    else
      {
	// Otherwise get applied field strength
	(*Applied_field_fct_pt)(x,H_applied);
      }
  }
  
  // EFFECTIVE ANISOTROPY FIELD
  /// Function pointer to crystalline anisotropy field function
  typedef void (*CrystAnisFieldFctPt)(const Vector<double>& x,Vector<double>& H_cryst_anis);

  /// Access function: Pointer to crystalline anisotropy field function
  CrystAnisFieldFctPt& cryst_anis_field_fct_pt() {return Cryst_anis_field_fct_pt;}

  /// Access function: Pointer to crystalline anisotropy field function. Const version
  CrystAnisFieldFctPt cryst_anis_field_fct_pt() const {return Cryst_anis_field_fct_pt;}

  /// Get the crystalline anisotropy field at Eulerian position x.
  inline virtual void get_H_cryst_anis_field(const Vector<double> &x, Vector<double> &H_cryst_anis) const
  {
    //If no crystalline anisotropy function has been set, return zero vector
    if(Cryst_anis_field_fct_pt==0) {for(unsigned j=0;j<3;j++) H_cryst_anis[j] = 0.0;}
    else
      {
	// Otherwise get crystalline anisotropy field strength
	(*Cryst_anis_field_fct_pt)(x,H_cryst_anis);
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

  /// Pointer to poisson source function - only for testing purposes since div(M) is our source function in calculation of the demagnetising potential.
  PoissonSourceFctPt Source_fct_pt;
  
  /// Pointer to function giving applied field.
  AppliedFieldFctPt Applied_field_fct_pt;

  /// Pointer to function giving effective field due to the crystalline anisotropy
  CrystAnisFieldFctPt Cryst_anis_field_fct_pt;


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
      // Interpolate x, m, phi, derrivatives:
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

	}

      // Get source function   
      double source;
      get_source_poisson(ipt,interpolated_x,source);

      // Get applied field at this position
      get_applied_field(interpolated_x, H_applied);
      
      // Get crystalline anisotropy effective field
      get_H_cryst_anis_field(interpolated_x, H_cryst_anis);
       
      // Take total of all fields used ??ds pass this entire section out to a function eventually if possible?
      // ??ds add 0.1 to push off maximum (i.e. thermal-ish...)
      for(unsigned j=0; j<3; j++)
	{
	  H_total[j] = H_cryst_anis[j] + H_applied[j] - 0.1;
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
	  
	  // Calculate residual for poisson equation (to find phi):
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
	      // ??ds add in jacobian calculation eventually
	    }

	  // Calculate residual for 

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
	      //??ds put in jacobian calculation eventually

	    } // End of micromagnetics section

	} // End of calculating residuals and jacobians

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
