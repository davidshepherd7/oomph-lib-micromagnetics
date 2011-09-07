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
  /// Constructor (must initialise the Source_fct_pt to null)
  MicromagEquations() : Source_fct_pt(0) {}
 
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
  unsigned required_nvalue(const unsigned &n) const
  {
    // Remember to change this when adding new fields.
    return 4;
  }
  
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


  // GET VALUES
  /// Return the i-th value stored at local node n but do NOT take hanging nodes into account
  double raw_nodal_value(const unsigned &n, const unsigned &i)
  {return node_pt(n)->raw_value(i);}
  

  // SOURCE FUNCTION
  // Function pointer to source function fct(x, f(x)), x is a vector
  typedef void (*MicromagSourceFctPt)(const Vector<double>& x, double& f);

  ///  Function pointer to gradient of source function  fct(x,g(x)), x is a Vector
  typedef void (*MicromagSourceFctGradientPt)(const Vector<double>& x, Vector<double>& gradient);

  /// Access function: Pointer to source function
  MicromagSourceFctPt& source_fct_pt() {return Source_fct_pt;}

  /// Access function: Pointer to source function. Const version
  MicromagSourceFctPt source_fct_pt() const {return Source_fct_pt;}

  /// Access function: Pointer to gradient of source function
  MicromagSourceFctGradientPt& source_fct_gradient_pt() {return Source_fct_gradient_pt;}

  /// Access function: Pointer to gradient source function. Const version
  MicromagSourceFctGradientPt source_fct_gradient_pt() const {return Source_fct_gradient_pt;}

  // Return the poisson source term at (Eulerian) position x. ??Replace with div(M) eventually??
  inline virtual void get_source_micromag(const unsigned& ipt, const Vector<double>& x, double& source) const
  {
    if(Source_fct_pt==0) {source = 0.0;}    // If no source function has been set, return zero
    else {(*Source_fct_pt)(x,source);}	// Else get source strength
  }

  /// Get gradient of source term at (Eulerian) position x. Computed via function pointer (if set) or by finite differencing (default)
  inline virtual void get_source_gradient_micromag(const unsigned& ipt, const Vector<double>& x, Vector<double>& gradient) const
  {
    //If no gradient function has been set, FD it
    if(Source_fct_gradient_pt==0)
      {
	// Reference value
	double source=0.0;
	get_source_micromag(ipt,x,source);

	// FD it
	double eps_fd=GeneralisedElement::Default_fd_jacobian_step;
	double source_pls=0.0;
	Vector<double> x_pls(x);
	for (unsigned i=0;i<DIM;i++)
	  {
	    x_pls[i]+=eps_fd;
	    get_source_micromag(ipt,x_pls,source_pls);
	    gradient[i]=(source_pls-source)/eps_fd;
	    x_pls[i]=x[i];
	  }
      }
    else
      {
	// Get gradient
	(*Source_fct_gradient_pt)(x,gradient);
      }
  }

  // OPERATIONS ON RESULTS
  /// Get flux: flux[i] = dphi/dx_i (H_demag = -1 * flux)
  void get_flux(const Vector<double>& s, Vector<double>& flux) const
  {
    //Find out how many nodes there are in the element
    const unsigned n_node = nnode();

    //Get the index at which the unknown is stored
    const unsigned phi_nodal_index = phi_index_micromag();

    //Set up memory for the shape and test functions
    Shape psi(n_node);
    DShape dpsidx(n_node,DIM);
 
    //Call the derivatives of the shape and test functions
    dshape_eulerian(s,psi,dpsidx);
     
    //Initialise to zero
    for(unsigned j=0;j<DIM;j++)
      {
	flux[j] = 0.0;
      }
   
    // Loop over nodes
    for(unsigned l=0;l<n_node;l++) 
      {
	//Loop over derivative directions
	for(unsigned j=0;j<DIM;j++)
	  {                               
	    flux[j] += this->nodal_value(l,phi_nodal_index)*dpsidx(l,j);
	  }
      }
  }

  /// Get divergence: divergence(M[i]) = dM_i/dx_i at local coordinate s (used as source function for poisson)
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
  virtual void get_dresidual_dnodal_coordinates(RankThreeTensor<double>&
						dresidual_dnodal_coordinates);

protected:

  // Pointer to source function, ?? should only be temporary
  MicromagSourceFctPt Source_fct_pt;

  /// Pointer to gradient of source function
  MicromagSourceFctGradientPt Source_fct_gradient_pt;

  /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
  virtual double dshape_and_dtest_eulerian_micromag(const Vector<double> &s, Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx) const=0;


  /// Shape, test functions & derivs. w.r.t. to global coords. at integration point ipt. Return Jacobian.
  virtual double dshape_and_dtest_eulerian_at_knot_micromag(const unsigned& ipt, Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx) const=0;


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

  //Index at which the poisson unknown is stored
  const unsigned phi_nodal_index = phi_index_micromag();
 
  //Set the value of n_intpt
  const unsigned n_intpt = integral_pt()->nweight();

  //Integers to store the local equation and unknown numbers
  int local_eqn=0, local_unknown=0;

  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
      //Get the integral weight
      double w = integral_pt()->weight(ipt);

      //Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_micromag(ipt,psi,dpsidx,
							    test,dtestdx);
       
      //Premultiply the weights and the Jacobian
      double W = w*J;

      //Calculate local values of unknown
      //Allocate and initialise to zero
      double interpolated_phi=0.0;
      Vector<double> interpolated_x(DIM,0.0);
      Vector<double> interpolated_dphidx(DIM,0.0);
   
      //Calculate function value and derivatives:
      //-----------------------------------------
      // Loop over nodes
      for(unsigned l=0;l<n_node;l++) 
	{
	  //Get the nodal value of the poisson unknown (phi)
	  double phi_value = raw_nodal_value(l,phi_nodal_index);
	  interpolated_phi += phi_value*psi(l);

	  // Loop over spatial directions
	  for(unsigned j=0;j<DIM;j++)
	    {
	      interpolated_x[j] += raw_nodal_position(l,j)*psi(l);
	      interpolated_dphidx[j] += phi_value*dpsidx(l,j);
	    }

	}


      //Get source function ?? in this case we use div(M) - hard coded, maybe bad?
      // div_m is calculated in terms of local coordinate - ?? not sure if this is a good choice?
      // ?? have to get s (local coordinate) just for this - might be a better way somehow?
      Vector<double> s(DIM,0);
      for(unsigned j=0; j<DIM; j++) 
	{
	  s[j] = integral_pt()->knot(ipt,j);
	}
      double source = divergence_M(s);
      cout << s[0] << source << "\t" << endl;
      //get_source_micromag(ipt,interpolated_x,source);

      // Assemble residuals and Jacobian
      //--------------------------------
       
      // Loop over the test functions
      for(unsigned l=0;l<n_node;l++)
	{
	  //Get the local equation
	  local_eqn = nodal_local_eqn(l,phi_nodal_index);
	  /*IF it's not a boundary condition*/
	  if(local_eqn >= 0)
	    {
	      // Add source term here 
	      residuals[local_eqn] += source*test(l)*W;
             
	      // The Poisson bit itself
	      for(unsigned k=0;k<DIM;k++)
		{
		  residuals[local_eqn] += interpolated_dphidx[k]*dtestdx(l,k)*W;
		}

	      // Calculate the jacobian
	      //-----------------------
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
			      jacobian(local_eqn,local_unknown) 
				+= dpsidx(l2,i)*dtestdx(l,i)*W;
			    }
			}
		    }
		}
	    }
	}

    } // End of loop over integration points

} // End of fill in residuals function 



  //======================================================================
  /// Compute derivatives of elemental residual vector with respect
  /// to nodal coordinates. 
  /// dresidual_dnodal_coordinates(l,i,j) = d res(l) / dX_{ij}
  /// Overloads the FD-based version in the FE base class.
  //======================================================================
template <unsigned DIM>
void  MicromagEquations<DIM>::get_dresidual_dnodal_coordinates(
							       RankThreeTensor<double>&
							       dresidual_dnodal_coordinates)
{

  //Find out how many nodes there are
  const unsigned n_node = nnode();

  //Set up memory for the shape and test functions
  Shape psi(n_node), test(n_node);
  DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);
  DShape dpsidx_pls(n_node,DIM), dtestdx_pls(n_node,DIM);

  // Deriatives of shape fct derivatives w.r.t. nodal coords
  RankFourTensor<double> d_dpsidx_dX(DIM,n_node,n_node,DIM);
  RankFourTensor<double> d_dtestdx_dX(DIM,n_node,n_node,DIM);

  // Derivative of Jacobian of mapping w.r.t. to nodal coords
  DenseMatrix<double> dJ_dX(DIM,n_node);

  // Derivatives of derivative of phi w.r.t. nodal coords
  RankThreeTensor<double> d_dphidx_dX(DIM,n_node,DIM);

  // Gradient of source fct
  Vector<double> d_source_dx(DIM);

  //Index at which the poisson unknown is stored
  const unsigned phi_nodal_index = phi_index_micromag();
 
  //Set the value of n_intpt
  const unsigned n_intpt = integral_pt()->nweight();

  //Integers to store the local equation number
  int local_eqn=0;

  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
      //Get the integral weight
      double w = integral_pt()->weight(ipt);

      //Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_micromag(ipt,psi,dpsidx,
							    test,dtestdx);
       
      //Calculate local values 
      //Allocate and initialise to zero
      Vector<double> interpolated_x(DIM,0.0);
      Vector<double> interpolated_dphidx(DIM,0.0);
   
      //Calculate function value and derivatives:
      //-----------------------------------------
      // Loop over nodes
      for(unsigned l=0;l<n_node;l++) 
	{
	  //Get the nodal value of the Poisson unknown
	  double phi_value = raw_nodal_value(l,phi_nodal_index);
	  // Loop over directions
	  for(unsigned j=0;j<DIM;j++)
	    {
	      interpolated_x[j] += raw_nodal_position(l,j)*psi(l);
	      interpolated_dphidx[j] += phi_value*dpsidx(l,j);
	    }
	}

      //Get source function
      //-------------------
      double source;
      get_source_micromag(ipt,interpolated_x,source);

      // FD step 
      double eps_fd=GeneralisedElement::Default_fd_jacobian_step;
   
      // Do FD loop
      for (unsigned jj=0;jj<n_node;jj++)
	{
	  // Get node
	  Node* nod_pt=node_pt(jj);
     
	  // Loop over coordinate directions
	  for (unsigned ii=0;ii<DIM;ii++)
	    {
	      // Make backup
	      double backup=nod_pt->x(ii);
       
	      // Do FD step. No node update required as we're
	      // attacking the coordinate directly...
	      nod_pt->x(ii)+=eps_fd;
       
	      //Call the derivatives of the shape and test functions
	      //at advanced level
	      double J_pls = 
		dshape_and_dtest_eulerian_at_knot_micromag(ipt,psi,dpsidx_pls,
							   test,dtestdx_pls);
       
	      // Assign
	      dJ_dX(ii,jj)=(J_pls-J)/eps_fd;
	      for (unsigned i=0;i<DIM;i++)
		{
		  for (unsigned j=0;j<n_node;j++)
		    {
		      d_dpsidx_dX(ii,jj,j,i)=(dpsidx_pls(j,i)-dpsidx(j,i))/eps_fd;
		      d_dtestdx_dX(ii,jj,j,i)=(dtestdx_pls(j,i)-dtestdx(j,i))/eps_fd;
		    }
		}

	      // Shape deriv of du/dx_i
	      for (unsigned i=0;i<DIM;i++)
		{
		  double aux=0.0;
		  for (unsigned j_nod=0;j_nod<n_node;j_nod++)
		    {
		      aux+=raw_nodal_value(j_nod,phi_nodal_index)*
			d_dpsidx_dX(ii,jj,j_nod,i);
		    }
		  d_dphidx_dX(ii,jj,i)=aux;
		}
  
	      // Reset coordinate. No node update required as we're
	      // attacking the coordinate directly...
	      nod_pt->x(ii)=backup;
	    }
	}

      // Get gradient of source function
      get_source_gradient_micromag(ipt,interpolated_x, d_source_dx);


      // Assemble shape derivatives
      //---------------------------
       
      // Loop over the test functions
      for(unsigned l=0;l<n_node;l++)
	{
	  //Get the local equation
	  local_eqn = nodal_local_eqn(l,phi_nodal_index);

	  /*IF it's not a boundary condition*/
	  if(local_eqn >= 0)
	    {
	      // Loop over coordinate directions
	      for (unsigned ii=0;ii<DIM;ii++)
		{              
		  // Loop over nodes
		  for (unsigned jj=0;jj<n_node;jj++)
		    {       
		      double sum=source*test(l)*dJ_dX(ii,jj)+
			d_source_dx[ii]*test(l)*psi(jj)*J;
         
		      for (unsigned k=0;k<DIM;k++)
			{
			  sum+=interpolated_dphidx[k]*(dtestdx(l,k)*dJ_dX(ii,jj)+
						       d_dtestdx_dX(ii,jj,l,k)*J)
			    + d_dphidx_dX(ii,jj,k)*dtestdx(l,k)*J;
			}

		      // Multiply through by integration weight
		      dresidual_dnodal_coordinates(local_eqn,ii,jj)+=sum*w;
		    }
		}
	    }
	}

    } // End of loop over integration points
}   



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
namespace ScalarPotential
{
  void source_function(const Vector<double>& M, double& source)
  {
    // Initially try with M_x = x so source fn = -div(M) == -1
    // source = -1;

    // Try with M_x  like:  _/\_ (i.e. up to one then back down), then source fn is step fn with step at x = 0.5, source = -1 at x = 0, source = 0 at x = 1.
    source = -1*tanh(500*(M[0] - 0.5));
  }

  void get_boundary_u(const Vector<double>& x, Vector<double>& u)
  {
    u[0] = 0;
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
  MicromagEquations<1>::MicromagSourceFctPt Source_fct_pt;

public:

  /// Constructor: Pass number of elements and pointer to source function
  OneDMicromagProblem(const unsigned& n_element, 
		      MicromagEquations<1>::MicromagSourceFctPt source_fct_pt);

  /// Destructor (empty -- all the cleanup is done in the base class)
  ~OneDMicromagProblem(){};

  /// Update the problem specs before solve: Set boundary conditions
  void actions_before_newton_solve();

  /// Update the problem specs after solve (calculate demag field)
  void actions_after_newton_solve(const unsigned n_element);

  /// Doc the solution
  void doc_solution(const unsigned& label);

  // // TIME STEPPING FUNCTIONS
  // /// Update the problem specs after solve (empty)
  // void actions_after_implicit_timestep() {}

  // /// Update the problem specs before next timestep
  // void actions_before_implicit_timestep();

  // /// Set initial condition (incl previous timesteps) according to specified function. 
  // void set_initial_condition();

}; // end of problem class





//=====start_of_constructor===============================================
/// \short Constructor for 1D Micromag problem in unit interval.
/// Discretise the 1D domain with n_element elements of type ELEMENT.
/// Specify function pointer to source function. 
//========================================================================
template<class ELEMENT>
OneDMicromagProblem<ELEMENT>::OneDMicromagProblem(const unsigned& n_element,
						  MicromagEquations<1>::MicromagSourceFctPt source_fct_pt) : 
  Source_fct_pt(source_fct_pt)
{ 
  // Set domain length 
  double L=1.0;

  // Build mesh and store pointer in Problem
  Problem::mesh_pt() = new OneDMesh<ELEMENT>(n_element,L);

  // Set the boundary conditions for this problem: By default, all nodal
  // values are free -- we only need to pin the ones that have 
  // Dirichlet conditions. 

  // Pin the single nodal value at the single node on mesh 
  // boundary 0 (= the left domain boundary at x=0)
  mesh_pt()->boundary_node_pt(0,0)->pin(0);
 
  // Pin the single nodal value at the single node on mesh 
  // boundary 1 (= the right domain boundary at x=1)
  mesh_pt()->boundary_node_pt(1,0)->pin(0);


  // Loop over elements to set pointers to source function and create element objects?
  for(unsigned i=0;i<n_element;i++)
    {
      // Upcast from GeneralisedElement to the present element
      ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
      //Set the source function pointer
      elem_pt->source_fct_pt() = Source_fct_pt;
    }

 
  // Setup equation numbering scheme
  assign_eqn_numbers();

} // end of constructor




  //===start_of_actions_before_newton_solve========================================
  /// \short Update the problem specs before solve: (Re)set boundary values
  /// from the exact solution. 
  //========================================================================
template<class ELEMENT>
void OneDMicromagProblem<ELEMENT>::actions_before_newton_solve()
{
 
  // Assign boundary values for this problem by reading them out
  // from the exact solution.

  // Left boundary is node 0 in the mesh:
  Node* left_node_pt=mesh_pt()->node_pt(0);

  // Determine the position of the boundary node (the exact solution
  // requires the coordinate in a 1D vector!)
  Vector<double> x(1);
  x[0]=left_node_pt->x(0);
 
  // Boundary value (read in from exact solution which returns
  // the solution in a 1D vector)
  Vector<double> u(1);
  ScalarPotential::get_boundary_u(x,u);
 
  // Assign the boundary condition to one (and only) nodal value
  left_node_pt->set_value(0,u[0]);


  // Right boundary is last node in the mesh:
  unsigned last_node=mesh_pt()->nnode()-1;
  Node* right_node_pt=mesh_pt()->node_pt(last_node);

  // Determine the position of the boundary node
  x[0]=right_node_pt->x(0);
 
  // Boundary value (read in from exact solution which returns
  // the solution in a 1D vector)
  ScalarPotential::get_boundary_u(x,u);
 
  // Assign the boundary condition to one (and only) nodal value
  right_node_pt->set_value(0,u[0]);

 
} // end of actions before solve


  //===start_of_doc=========================================================
  /// Doc the solution in tecplot format. Label files with label.
  //========================================================================
template<class ELEMENT>
void OneDMicromagProblem<ELEMENT>::doc_solution(const unsigned& label)
{ 

  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts;
  npts=5; 

  // Output solution with specified number of plot points per element
  sprintf(filename,"soln%i.dat",label);
  some_file.open(filename);
  mesh_pt()->output(some_file,npts);
  some_file.close();

  // // Output exact solution at much higher resolution (so we can
  // // see how well the solutions agree between nodal points)
  // sprintf(filename,"exact_soln%i.dat",label);
  // some_file.open(filename);
  // mesh_pt()->output_fct(some_file,20*npts,ScalarPotential::get_exact_u); 
  // some_file.close();

  // // Doc pointwise error and compute norm of error and of the solution
  // double error,norm;
  // sprintf(filename,"error%i.dat",label);
  // some_file.open(filename);
  // mesh_pt()->compute_error(some_file,ScalarPotential::get_exact_u,
  //                          error,norm); 
  // some_file.close();

  // // Doc error norm:
  // cout << "\nNorm of error    : " << sqrt(error) << std::endl; 
  // cout << "Norm of solution : " << sqrt(norm) << std::endl << std::endl;
  // cout << std::endl;

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
  // Solve a 1D Micromag problem using a source function that generates
  // a fish shaped exact solution
  unsigned n_element=20; //Number of elements
  OneDMicromagProblem<QMicromagElement<1,2> >problem(n_element,ScalarPotential::source_function);
  //Element type as template parameter (num dimension, num nodes as other parameters)

  // Set initial condition (on M) - not yet implemented ??
  problem.set_initial_condition();


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

  // Solve the problem
  problem.newton_solve(n_element);

  //Output solution for this case (label output files with "1")
  problem.doc_solution(1);
  
} // end of main
