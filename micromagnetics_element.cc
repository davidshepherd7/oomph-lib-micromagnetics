#ifndef OOMPH_MICROMAGNETICS_ELEMENT_H
#include "micromagnetics_element.h"
#define OOMPH_MICROMAGNETICS_ELEMENT_H

using namespace oomph;
using namespace MathematicalConstants;

namespace oomph
{

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
    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);

    // Set up vector to store the local coordinates
    Vector<double> s(DIM);

    // Get current time
    double time = time_pt()->time();

    // Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    // Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
      {
  	// Get the integral weight
  	double w = integral_pt()->weight(ipt);

  	// Call the derivatives of the shape and test functions
  	double J = dshape_and_dtest_eulerian_at_knot_micromag(ipt,psi,dpsidx,test,dtestdx);

  	//Premultiply the weights and the Jacobian
  	double W = w*J;

  	// Get values of s (local coordinate)
  	for(unsigned j=0; j<DIM; j++) {s[j] = integral_pt()->knot(ipt,j);}

  	// Allocate memory for local quantities and initialise to zero. dphidx
  	// is also H_demag so we need all 3 components.
  	double interpolated_phi(0.0), interpolated_phi_1(0.0), interpolated_divm(0.0);
  	Vector<double> interpolated_x(DIM,0.0), interpolated_dphidx(3,0.0),
	  interpolated_dphi_1dx(3,0.0),interpolated_m(3,0.0), interpolated_dmdt(3,0.0),
	  interpolated_H_exchange(3,0.0);
	DenseDoubleMatrix interpolated_dmdx(3,3,0.0);

  	// Interpolate values at knot by looping over nodes adding contributions
  	for(unsigned l=0;l<n_node;l++)
  	  {
  	    // Interpolate values of the total phi and phi_1
  	    interpolated_phi += nodal_value(l,phi_index_micromag()) * psi(l);
	    interpolated_phi_1 += nodal_value(l,phi_1_index_micromag()) * psi(l);

	    // Interpolate values of fields and magnetisations
	    Vector<double> dmdt(3,0.0);
  	    dm_dt_micromag(l,dmdt);
  	    for(unsigned j=0; j<3; j++)
  	      {
  		interpolated_dmdt[j] += dmdt[j]*psi(l);
  		interpolated_m[j] += nodal_value(l,m_index_micromag(j))*psi(l);
  		interpolated_H_exchange[j] += nodal_value(l,exchange_index_micromag(j))*psi(l);
  	      }

  	    // Interpolate spatial values/derivatives
  	    for(unsigned j=0; j<DIM; j++)
  	      {
  		interpolated_x[j] += nodal_position(l,j)*psi(l);
  		interpolated_dphidx[j] += nodal_value(l,phi_index_micromag())*dpsidx(l,j);
		interpolated_dphi_1dx[j] += nodal_value(l,phi_1_index_micromag())*dpsidx(l,j);
		for(unsigned k=0; k<3; k++)
		  interpolated_dmdx(k,j) += nodal_value(l,m_index_micromag(k))*dpsidx(l,j);
  	      }

	    // Calculate divergence of m
	    for(unsigned j=0; j<DIM; j++)
	      interpolated_divm += interpolated_dmdx(j,j);
  	  }

  	// Total potential (magnetostatic field calculations)
  	//----------------------------------------------------
  	// // Get source function
  	// double poisson_source = 0;
  	// get_poisson_source(time,ipt,interpolated_x,poisson_source);

  	// Loop over the test functions/nodes adding contributions
  	for(unsigned l=0;l<n_node;l++)
  	  {
  	    // Get the local equation numbers, check if it's a boundary condition
  	    int phi_local_eqn = nodal_local_eqn(l,phi_index_micromag());
  	    if(phi_local_eqn >= 0)
  	      {
  		// // Add source term
  		// residuals[phi_local_eqn] -= poisson_source*test(l)*W; //??ds not sure on sign

		// The divergence of M source term
		residuals[phi_local_eqn] -= interpolated_divm*test(l)*W;

  		// The Poisson part (after reducing the order of differentiation
  		// on phi using integration by parts gives a dot product).
  		for(unsigned k=0;k<DIM;k++)
  		  {
  		    residuals[phi_local_eqn] -= interpolated_dphidx[k]*dtestdx(l,k)*W;
  		  }

  		// ??ds add in jacobian calculation eventually
  	      }
  	  }
	// End of total potential section



	// Reduced potential (to get boundary conditions on total potential)
	//----------------------------------------------------------------------
	// The only difference between this and the total potential section is in
	// the boundary conditions.

  	// Loop over the test functions/nodes adding contributions
  	for(unsigned l=0;l<n_node;l++)
  	  {
  	    // Get the local equation numbers, check if it's a boundary condition
  	    int phi_1_local_eqn = nodal_local_eqn(l,phi_1_index_micromag());
  	    if(phi_1_local_eqn >= 0)
  	      {
		// The divergence of M source term
		residuals[phi_1_local_eqn] -= interpolated_divm*test(l)*W;

  		// The Poisson part (after reducing the order of differentiation
  		// on phi_1 using integration by parts gives a dot product).
  		for(unsigned k=0;k<DIM;k++)
  		  {
  		    residuals[phi_1_local_eqn] -= interpolated_dphi_1dx[k]*dtestdx(l,k)*W;
  		  }

  		// ??ds add in jacobian calculation eventually
  	      }
  	  }



  	// Exchange field section
  	//----------------------------------------------------
  	// H_ex = coeff * laplacian(M)
  	// Weak form: int( H_ex_i - coeff* grad(M).grad(test) ) = 0
  	// ??ds only when grad(M).n = 0 at boundaries, otherwise need another term!

  	// Get exchange coeff
  	double exchange_coeff = get_exchange_coeff(time,interpolated_x);

  	// Loop over H_exchange components
  	for(unsigned i=0; i<3; i++)
  	  {
  	    // Loop over test functions/nodes, adding contributions from each
  	    for(unsigned l=0; l<n_node; l++)
  	      {
  		// Get local equation number for H_exchange
  		int exchange_local_eqn = nodal_local_eqn(l,exchange_index_micromag(i));

  		// If it's not a boundary condition
  		if(exchange_local_eqn >=0)
  		  {
  		    // Add exchange field component at integration pt
  		    residuals[exchange_local_eqn] += interpolated_H_exchange[i] * test(l) * W;

  		    // Get grad(M).grad(test) at integration pt
		    double gradMi_dot_gradpsi(0.0);
  		    for(unsigned j=0; j<DIM; j++)
  		      {
			gradMi_dot_gradpsi += interpolated_dmdx(i,j) * dtestdx(l,j);
  		      }
		    residuals[exchange_local_eqn]
		      += exchange_coeff * gradMi_dot_gradpsi * W;

  		    //?? jacobian calculation here
  		  }
  	      } // End of loop over test functions/nodes
  	  } // End of loop over H directions



  	// LLG section (time evolution of magnetisation)
  	//----------------------------------------------------

  	// Get applied field at this position
  	Vector<double> H_applied(3,0.0);
  	get_applied_field(time, interpolated_x, H_applied);

  	// Get crystalline anisotropy effective field
  	Vector<double> H_cryst_anis(3,0.0);
  	get_H_cryst_anis_field(time, interpolated_x, interpolated_m, H_cryst_anis);

  	// // Get LLG source function
  	// Vector<double> llg_source(3,0.0);
  	// get_source_llg(time, interpolated_x, llg_source);

  	// Take total of all fields used (-dphidx = magnetostatic field)
  	Vector<double> H_total(3,0.0);
  	for(unsigned j=0; j<3; j++)
  	  {
  	    H_total[j] = H_applied[j] - interpolated_dphidx[j]
  	      + interpolated_H_exchange[j] + H_cryst_anis[j];
  	  }

  	// Get the cross products for the LLG equation
	Vector<double> interpolated_mxH(3,0.0), interpolated_mxmxH(3,0.0);
  	cross(interpolated_m, H_total, interpolated_mxH);
  	cross(interpolated_m, interpolated_mxH, interpolated_mxmxH);

	// Loop over test functions/nodes
  	for (unsigned l=0; l<n_node; l++)
  	  {
  	    // Loop over components of m
  	    for(unsigned i=0; i<3; i++)
  	      {
  		// Get the local equation number for the ith component of m
  		int m_local_eqn = nodal_local_eqn(l,m_index_micromag(i));

  		if(m_local_eqn >= 0)  // If it's not a boundary condition
  		  {
  		    residuals[m_local_eqn] +=
  		      ( interpolated_dmdt[i]
  			+ get_llg_precess() * interpolated_mxH[i]
  			+ get_llg_damp() * interpolated_mxmxH[i]
  			// - llg_source[i]
  			)*test(l)*W;
  		  }
  		//??ds put in jacobian calculation eventually

  	      }
  	  } // End of loop over test functions
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

    // Get time
    //double t = time_pt()->time();

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

	// Output solution vector at local coordinate s
	Vector<double> interpolated_solution(7,0.0); //??ds generalise the length?
	interpolated_solution_micromag(s,interpolated_solution);
	for(unsigned i=0; i<7; i++)
	  {
	    outfile << interpolated_solution[i] << " ";
	  }

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
    	Vector<double> exact_solution(7,0.0);
    	(*exact_soln_pt)(time,x,exact_solution);
    	for(unsigned i=0; i<7; i++)
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

	// Get x (global coordinate) and output (??ds at current time)
	Vector<double> x(DIM,0.0);
	interpolated_x(s,x);
	for(unsigned i=0; i<DIM; i++){outfile << x[i] << " ";}

	//Get the integral weight
	double w = integral_pt()->weight(ipt);

	// Get jacobian of mapping
	double J=J_eulerian(s);

	//Premultiply the weights and the Jacobian
	double W = w*J;

	// Get entire interpolated solution at position s and current time
	Vector<double> interpolated_soln(nvalues,0.0);
	interpolated_solution_micromag(s,interpolated_soln);

	// Get entire exact solution at point x and time "time"
	Vector<double> exact_soln(nvalues,0.0);
	(*exact_soln_pt)(time,x,exact_soln);

	// Output the error (difference between exact and interpolated solutions)
	for(unsigned i=0; i<nvalues; i++){outfile << exact_soln[i]- interpolated_soln[i] << " ";}
	outfile << std::endl;

	// Add contributions to the norms of the error and exact soln from this integration point
	for(unsigned i=0; i<nvalues; i++)
	  {
	    error_norm += (exact_soln[i] - interpolated_soln[i])
	      *(exact_soln[i] - interpolated_soln[i])*W;
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
  //       fprintf(file_pt,"%g \n",interpolated_phi_micromag(s));
  //     }

  //   // Write tecplot footer (e.g. FE connectivity lists)
  //   write_tecplot_zone_footer(file_pt,nplot);
  // }

  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////


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

}

#endif
