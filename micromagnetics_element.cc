#ifndef OOMPH_MICROMAGNETICS_ELEMENT_H
#include "micromagnetics_element.h"
#define OOMPH_MICROMAGNETICS_ELEMENT_H

using namespace oomph;
using namespace MathematicalConstants;

namespace oomph
{

  //======================================================================
  /// Set the data for the number of variables at each node, it is the same in
  /// all element types regardless of dimension, element order and geometry.
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  const unsigned QMicromagElement<DIM,NNODE_1D>::Initial_Nvalue = 5;

  template<unsigned DIM, unsigned NNODE_1D>
  const unsigned TMicromagElement<DIM,NNODE_1D>::Initial_Nvalue = 5;

  //======================================================================
  /// Compute element residual Vector and/or element Jacobian matrix
  ///
  /// flag=1: compute both
  /// flag=0: compute only residual Vector
  ///
  /// Pure version without hanging nodes
  //======================================================================
  template<unsigned DIM>
  void MicromagEquations<DIM>::fill_in_generic_residual_contribution_micromag
  (Vector<double> &residuals, DenseMatrix<double> &jacobian,
   const unsigned& flag) const
  {

    //======================================================================
    /// Get some useful numbers and set up memory.
    //======================================================================

    // Find out how many nodes there are
    const unsigned n_node = nnode();
    // Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);
    // Get current time
    const double time = time_pt()->time();
    // Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();


    //======================================================================
    /// Begin loop over the knots (integration points)
    //======================================================================
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
      {
	//======================================================================
	/// Calculate/get/interpolate all values for the residual calculations
	//======================================================================
	Vector<double> s(DIM);
  	for(unsigned j=0; j<DIM; j++) {s[j] = integral_pt()->knot(ipt,j);}
  	double J = dshape_dtest(s,psi,dpsidx,test,dtestdx);
  	double W = integral_pt()->weight(ipt) * J;

  	// Allocate memory for local quantities and initialise to zero. dphidx
  	// is also H_demag so we need all 3 components.
  	double itp_phi(0.0), itp_phi_1(0.0);
  	Vector<double> itp_x(DIM,0.0), itp_dphidx(3,0.0),
	  itp_dphi_1dx(3,0.0),itp_m(3,0.0), itp_dmdt(3,0.0);
	DenseDoubleMatrix itp_dmdx(3,3,0.0);

  	// Interpolate values at knot by looping over nodes adding contributions
  	for(unsigned l=0;l<n_node;l++)
  	  {
  	    itp_phi += nodal_value(l,phi_index_micromag()) * psi(l);
	    itp_phi_1 += nodal_value(l,phi_1_index_micromag()) * psi(l);

	    Vector<double> dmdt(3,0.0);
	    dm_dt_micromag(l,dmdt); // get dmdt at node l

  	    for(unsigned j=0; j<3; j++)
  	      {
  		itp_dmdt[j] += dmdt[j]*psi(l);
		itp_m[j] += nodal_value(l,m_index_micromag(j))*psi(l);
  	      }

  	    // Interpolate spatial values/derivatives
  	    for(unsigned j=0; j<DIM; j++)
  	      {
  		itp_x[j] += nodal_position(l,j)*psi(l);
  		itp_dphidx[j] += nodal_value(l,phi_index_micromag())*dpsidx(l,j);
		itp_dphi_1dx[j] += nodal_value(l,phi_1_index_micromag())*dpsidx(l,j);
		for(unsigned k=0; k<3; k++)
		  itp_dmdx(k,j) += nodal_value(l,m_index_micromag(k))*dpsidx(l,j);
  	      }
  	  }

	// Calculate other things (these need to go outside the loop over
	// nodes so that we have finished calculating itp_dmdx etc.).
	double itp_divm = 0.0;
	for(unsigned j=0; j<DIM; j++) itp_divm += itp_dmdx(j,j);
	Vector<double> itp_mxdmdt(3,0.0);
	VectorOps::cross(itp_m, itp_dmdt, itp_mxdmdt);

	// Source functions (for debugging, normally zero)
  	double phi_source = get_phi_source(time,itp_x);
  	double phi_1_source = get_phi_1_source(time,itp_x);

  	// Coefficients
  	double exch_c = get_exchange_coeff(time,itp_x);
	double llg_precess_c = get_llg_precession_coeff();
	double llg_damp_c = get_llg_damping_coeff();
	double magstatic_c = get_magnetostatic_coeff(time,itp_x);

	// Fields
	Vector<double> h_applied(3,0.0);
	get_applied_field(time, itp_x, h_applied);
	Vector<double> h_cryst_anis(3,0.0);
	get_H_cryst_anis_field(time, itp_x, itp_m, h_cryst_anis);


	//======================================================================
	/// Use the above values to calculate the residuals
	//======================================================================

  	// Loop over the test functions/nodes adding contributions to residuals
  	for(unsigned l=0;l<n_node;l++)
  	  {

	    // Total potential (phi)
  	    int phi_eqn = nodal_local_eqn(l,phi_index_micromag());
  	    if(phi_eqn >= 0) // if value is not pinned
  	      {
  		residuals[phi_eqn] -= phi_source*test(l)*W; // source
		residuals[phi_eqn] -= itp_divm*test(l)*W;	  // div(m)
  		for(unsigned k=0;k<DIM;k++)			  // Poisson
		  residuals[phi_eqn] -= itp_dphidx[k]*dtestdx(l,k)*W;
	      }

	    // Reduced potential (phi_1), only difference is in b.c.s
	    int phi_1_eqn = nodal_local_eqn(l,phi_1_index_micromag());
	    if(phi_1_eqn >= 0)
	      {
		residuals[phi_1_eqn] -= phi_1_source*test(l)*W;
		residuals[phi_1_eqn] -= itp_divm*test(l)*W;
		for(unsigned k=0;k<DIM;k++)
		  residuals[phi_1_eqn] -= itp_dphi_1dx[k]*dtestdx(l,k)*W;
	      }

	    // LLG itself (m, time evolution)
	    //=================================================

	    // Exchange residual contribution after integration by parts:
	    Vector<double> gradtestdotdmidx(3,0.0); //?? possibly could optimise
						    //- this is calculated twice
	    for(unsigned j=0; j<DIM; j++)
	      for(unsigned i=0; i<3; i++)
		gradtestdotdmidx[j] = dtestdx(l,i) * itp_dmdx(j,i);

	    // Total of field residual contributions.
	    // Magnetostatic field is -1* dphi/dx, coeff for debugging.
	    Vector<double> h_total(3,0.0);
	    for(unsigned j=0; j<DIM; j++)
	      h_total[j] = h_applied[j] + h_cryst_anis[j]
		- magstatic_c * itp_dphidx[j]
		- exch_c * gradtestdotdmidx[j];

	    // Cross product for the LLG equation
	    Vector<double> itp_mxh(3,0.0);
	    VectorOps::cross(itp_m, h_total, itp_mxh);

	    // add to residual
	    for(unsigned i=0; i<3; i++)
	      {
		int m_eqn = nodal_local_eqn(l,m_index_micromag(i));
		if(m_eqn >= 0)  // If it's not a boundary condition
		  {
		    residuals[m_eqn] += ( itp_dmdt[i]
					  + llg_precess_c * itp_mxh[i]
					  + llg_damp_c *itp_mxdmdt[i] )*test(l)*W;
		  }
	      }

	  } // End of loop over test functions, end of residual calculations


	//======================================================================
	/// If we want the Jacobian as well then calculate it
	//======================================================================
	if(flag){
	  // Double loop over nodes for the jacobian
	  for(unsigned l=0; l<n_node; l++){
	    for(unsigned l2=0;l2<n_node;l2++){

	      // Pre-calculations
	      Vector<double> gradpsil2(3,0.0), gradtestl(3,0.0),
		gradtestl2(3,0.0), itpmxgradpsil2(3,0.0);
	      for(unsigned j=0; j<DIM; j++){
		gradpsil2[j] = dpsidx(l2,j);
		gradtestl[j] = dtestdx(l,j);
		// gradtestl2[j] = dtestdx(l2,j);
	      }
	      double gradtestldotgradpsil2 = VectorOps::dot(gradtestl,gradpsil2);
	      VectorOps::cross(itp_m,gradpsil2,itpmxgradpsil2);

	      Vector<double> gradtestldotgradmi(3,0.0);
	      for(unsigned i=0; i<3; i++)
		for(unsigned j=0; j<DIM; j++)
		  gradtestldotgradmi[i] +=gradtestl[j] * itp_dmdx(j,i);

	      //=========================================================
	      /// Actual Jacobian calculation
	      //=========================================================

	      // Total potential (phi)
	      int phi_eqn = nodal_local_eqn(l,phi_index_micromag());
	      if(phi_eqn >= 0){

		// w.r.t. phi
		int phi_unknown = nodal_local_eqn(l2,phi_index_micromag());
		if(phi_unknown >= 0)
		  jacobian(phi_eqn,phi_unknown) += -gradtestldotgradpsil2 * W;

		// w.r.t. m
		for(unsigned j=0; j<3; j++){
		  int m_unknown = nodal_local_eqn(l2,m_index_micromag(j));
		  if(m_unknown >= 0)
		    jacobian(phi_eqn,m_unknown) += - dpsidx(l2,j) * test(l) * W;
		}

		// nothing w.r.t. phi_1
	      }


	      // Reduced potential (phi_1), only difference is in b.c.s
	      int phi_1_eqn = nodal_local_eqn(l,phi_1_index_micromag());
	      if(phi_1_eqn >= 0){

		// w.r.t. phi_1
		int phi_1_unknown = nodal_local_eqn(l2,phi_1_index_micromag());
		if(phi_1_unknown >= 0)
		  jacobian(phi_1_eqn,phi_1_unknown) += - gradtestldotgradpsil2 * W;

		// w.r.t. m
		for(unsigned j=0; j<3; j++){
		  int m_unknown = nodal_local_eqn(l2,m_index_micromag(j));
		  if(m_unknown >= 0)
		    jacobian(phi_1_eqn,m_unknown) += - dpsidx(l2,j) * test(l) * W;
		}

		// nothing w.r.t. phi
	      }


	      // ??ds need llg derrivatives in here
	      // llg derivatives
	      Vector<int> m_eqn(3), m_unknown(3);
	      for(unsigned j=0; j<3; j++)
		{
		  m_eqn[j] = nodal_local_eqn(l,m_index_micromag(j));
		  m_unknown[j] = nodal_local_eqn(l2,m_index_micromag(j));
		}


	      // w.r.t. phi
	      int phi_unknown = nodal_local_eqn(l2,phi_index_micromag());
	      if(phi_unknown >= 0){
		for(unsigned j=0; j<3; j++){
		  if(m_eqn[j] >= 0)
		    jacobian(m_eqn[j],phi_unknown) += itpmxgradpsil2[j] * W;
		}
	      }

	      // nothing w.r.t. phi1

	      // w.r.t. m
	      for(unsigned j=0; j<3; j++)
		{
		  if(m_unknown[j] >= 0){

		    // Construct the derivtives (we assume that almost always
		    // there will be some m[i] not pinned so don't bother to
		    // check yet).

		    Vector<double> jhatpsil2(3,0.0), jhatgradtestldotgradpsil2(3,0.0);
		    jhatpsil2[j] = psi(l2);
		    jhatgradtestldotgradpsil2[j] = gradtestldotgradpsil2;

		    Vector<double> integrand1(3,0.0), integrand2(3,0.0);
		    VectorOps::cross(jhatpsil2,gradtestldotgradmi,integrand1);
		    VectorOps::cross(itp_m,jhatgradtestldotgradpsil2,integrand1);

		    // Add each component:
		    for(unsigned i=0; i<3; i++)
		      {
			if(m_eqn[i] >= 0){
			  jacobian(m_eqn[i],m_unknown[j])
			    += - exch_c * W * ( integrand1[j] + integrand2[j]);
			}
		      }
		  }
		}


	    }
	  }
	}// End of Jacobian calculations

      }// End of loop over integration points
  } // End of fill in residuals function


  //======================================================================
  /// Output function:
  ///
  /// format is x, (y, z), phi, Mx, My, Mz, exact phi, exact Mx, exact My, exact Mz
  ///
  ///
  //======================================================================
template <unsigned DIM>
void MicromagEquations<DIM>::output(std::ostream &outfile,
				    const unsigned &n_plot)
{
  //Vector of local coordinates
  Vector<double> s(DIM);

  // Get number of values in solution at node 0
  // const unsigned nvalue = required_nvalue(0);

  // Tecplot header info
  outfile << tecplot_zone_string(n_plot);

  // Loop over plot points
  unsigned num_plot_points=nplot_points(n_plot);
  for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {

      // Get local coordinates of plot point
      get_s_plot(iplot,n_plot,s);

      // Get and output eulerian coordinates of plot point and output
      Vector<double> x(DIM,0.0);
      for(unsigned i=0; i<DIM; i++)
	{
	  x[i] = interpolated_x(s,i);
	  outfile << x[i] << " ";
	}

      // Output the magnetostatic field (= - dphidx) at this point
      Vector<double> itp_dphidx(3,0.0);
      interpolated_dphidx_micromag(s,itp_dphidx);
      for(unsigned i=0; i<3; i++)
	outfile << -itp_dphidx[i] << " ";

      // Phi 1 at this point
      outfile << interpolated_phi_1_micromag(s) << " ";

      // Phi total at this point
      outfile << interpolated_phi_micromag(s) << " ";

      // Output m at this point
      Vector<double> itp_m(3,0.0);
      interpolated_m_micromag(s,itp_m);
      for(unsigned i=0; i<3; i++)
	outfile << itp_m[i] << " ";

      // End the line ready for next point
      outfile << std::endl;
    }

  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(outfile,n_plot);

}

/// Output a time-dependent exact solution over the element.
template <unsigned DIM> void MicromagEquations<DIM>::
output_fct(std::ostream &outfile, const unsigned &n_plot,
	   const double& time,
	   const FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
{

  //Vector of local coordinates
  Vector<double> s(DIM);

  // Get number of values in solution at node 0
  const unsigned nvalue = required_nvalue(0);

  // Tecplot header info
  outfile << tecplot_zone_string(n_plot);

  // Loop over plot points
  unsigned num_plot_points=nplot_points(n_plot);
  for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {

      // Get local coordinates of plot point
      get_s_plot(iplot,n_plot,s);

      // Get and output eulerian coordinates of plot point and output
      Vector<double> x(DIM,0.0);
      for(unsigned i=0; i<DIM; i++)
	{
	  x[i] = interpolated_x(s,i);
	  outfile << x[i] << " ";
	}

      // Calculate and output exact solution at point x and time t
      Vector<double> exact_solution(nvalue,0.0);
      (*exact_soln_pt)(time,x,exact_solution);
      for(unsigned i=0; i<nvalue; i++)
	{
	  outfile << exact_solution[i] << " ";
	}

      // End the line ready for next point
      outfile << std::endl;
    }

  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(outfile,n_plot);
}


//======================================================================
/// Validate computed M against exact solution.
///
/// Solution is provided via function pointer.
/// Plot error at integration points.
///
//======================================================================
template <unsigned DIM>
void MicromagEquations<DIM>::
compute_error(std::ostream &outfile,
	      FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
	      const double& time, double& error_norm, double& exact_norm)
{

  // Initialise error and norm
  error_norm = 0.0;
  exact_norm = 0.0;

  // Find out how many nodes there are in the element
  unsigned n_node = nnode();

  // Get the shape function
  Shape psi(n_node);

  //Set the value of n_intpt
  unsigned n_intpt = integral_pt()->nweight();

  // Set how many values are output at each node, any will do so use 0th node.
  unsigned nvalues = node_pt(0)->nvalue();

  // Tecplot header info
  outfile << tecplot_zone_string(2);
  //??ds can't just use 2 here but seeing if it works
  //??ds use n_intpt/DIM? - only good for simple shaped elements?
  //??ds causes small issuse with output - points for soln/exact are at the corners of elements but for errors they are at the int pts

  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
      // Get s (local coordinate)
      Vector<double> s(DIM,0.0);
      for(unsigned i=0; i<DIM; i++) {s[i] = integral_pt()->knot(ipt,i);}

      // Get x (global coordinate) and output
      Vector<double> x(DIM,0.0);
      interpolated_x(s,x);
      for(unsigned i=0; i<DIM; i++){outfile << x[i] << " ";}

      //Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J=J_eulerian(s);

      //Premultiply the weights and the Jacobian
      double W = w*J;

      // Get entire ipl solution at position s and current time
      Vector<double> itp_soln(nvalues,0.0);
      interpolated_solution_micromag(s,itp_soln);

      // Get entire exact solution at point x and time "time"
      Vector<double> exact_soln(nvalues,0.0);
      (*exact_soln_pt)(time,x,exact_soln);

      // Output the error (difference between exact and itp solutions)
      for(unsigned i=0; i<nvalues; i++){outfile << exact_soln[i]- itp_soln[i] << " ";}
      outfile << std::endl;

      // Add contributions to the norms of the error and exact soln from this integration point
      for(unsigned i=0; i<nvalues; i++)
	{
	  error_norm += (exact_soln[i] - itp_soln[i])
	    *(exact_soln[i] - itp_soln[i])*W;
	  exact_norm += exact_soln[i]*exact_soln[i]*W;
	}

    }

  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(outfile,2);
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
//       fprintf(file_pt,"%g \n",itp_phi_micromag(s));
//     }

//   // Write tecplot footer (e.g. FE connectivity lists)
//   write_tecplot_zone_footer(file_pt,nplot);
// }

}

#endif
