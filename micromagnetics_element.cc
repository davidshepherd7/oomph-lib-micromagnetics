#ifndef OOMPH_MICROMAGNETICS_ELEMENT_H
#include "micromagnetics_element.h"
#define OOMPH_MICROMAGNETICS_ELEMENT_H

using namespace oomph;
using namespace MathematicalConstants;

namespace oomph
{

  //======================================================================
  /// Set the data for the number of Variables at each node.
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  const unsigned QMicromagElement<DIM,NNODE_1D>::Initial_Nvalue = 8;

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

    // Get current time
    const double time = time_pt()->time();

    // Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    // Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
      {
  	// Get values of s (local coordinate)
	Vector<double> s(DIM);
  	for(unsigned j=0; j<DIM; j++) {s[j] = integral_pt()->knot(ipt,j);}

  	// Call the derivatives of the shape and test functions
  	double J = dshape_dtest(s,psi,dpsidx,test,dtestdx);

  	// Premultiply the integration weight and the Jacobian of the coordinate transform
  	double W = integral_pt()->weight(ipt) * J;

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
	    dm_dt_micromag(l,dmdt); // get dmdt at node l
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
  	  }

	// Calculate divergences (this needs to go outside the loop over
	// nodes so that we have finished calculating interpolated_dmdx).
	for(unsigned j=0; j<DIM; j++)
	  {
	    interpolated_divm += interpolated_dmdx(j,j);
	  }

	// for(unsigned l=0; l<n_node; l++)
	//   for(unsigned i=0; i<3; i++)
	//     std::cout << nodal_value(l,m_index_micromag(i))*psi(l) << " ";

	//     std::cout <<  " " << interpolated_m << std::endl;

	// std::cout << interpolated_phi_1 << "\n"
	// 	  << interpolated_phi << "\n"
	// 	  << interpolated_divm << "\n"
	// 	  << interpolated_x << "\n"
	// 	  << interpolated_dphidx << "\n"
	// 	  << interpolated_dphi_1dx << "\n"
	// 	  << interpolated_m << "\n"
	// 	  << interpolated_dmdt << "\n"
	// 	  << interpolated_H_exchange << std::endl << std::endl;

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
  		// residuals[phi_local_eqn] -= poisson_source*test(l)*W;

		// The divergence of M source term
		residuals[phi_local_eqn] -= interpolated_divm*test(l)*W;

  		// The Poisson part (after reducing the order of differentiation
  		// on phi using integration by parts gives a dot product).
  		for(unsigned k=0;k<DIM;k++)
  		  {
  		    residuals[phi_local_eqn] -= interpolated_dphidx[k]*dtestdx(l,k)*W;
  		  }

		// // Calculate the Jacobian
		// if(flag)
		//   {
		//     // Loop over test functions/nodes again
		//     for(unsigned l2=0;l2<n_node;l2++)
		//       {
		// 	// First the dependence of phi on all nodal values of
		// 	// phi in the element.
		// 	int phi_local_unknown = nodal_local_eqn(l2,phi_index_micromag());
		// 	if(phi_local_unknown >= 0)
		// 	  {
		// 	    for(unsigned k=0; k<DIM; k++) // grad(psi_l) . grad(psi_l2)
		// 	      jacobian(phi_local_eqn,phi_local_unknown)
		// 		-= dtestdx(l,k)*dpsidx(l2,k)*W;
		// 	  }

		// 	// The dependence of phi on M
		// 	double div_psi = 0.0; // Get divergence of shape fn for node l2
		// 	for(unsigned k=0; k<DIM; k++)
		// 	  div_psi += dpsidx(l2,k);

		// 	for(unsigned j=0; j<3; j++)
		// 	  {
		// 	    int m_local_unknown = nodal_local_eqn(l2,m_index_micromag(j));
		// 	    if(m_local_unknown >= 0)
		// 	      jacobian(phi_local_eqn,m_local_unknown)
		// 		+= div_psi * test(l) * W;// ??ds add jacobian bit...
		// 	  }

		// 	// No dependence of phi on exchange field

		// 	// Phi depends on phi_1 via the boundary element
		// 	// matrix. It relates all nodes to each other rather
		// 	// than just the nodes within a single element. Hence
		// 	// this part of the Jacobian is added seperately in ??ds
		// 	// actions_before....
		//       }
		//   } // End of phi Jacobian calculation
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

		// // Calculate the Jacobian - pretty much exactly the same as for phi
		// if(flag)
		//   {
		//     // Loop over test functions/nodes again
		//     for(unsigned l2=0;l2<n_node;l2++)
		//       {
		// 	// First the dependence of phi_1 on all nodal values of
		// 	// phi_1 in the element.
		// 	int phi_1_local_unknown = nodal_local_eqn(l2,phi_1_index_micromag());
		// 	if(phi_1_local_unknown >= 0)
		// 	  {
		// 	    for(unsigned k=0; k<DIM; k++) // grad(psi_l) . grad(psi_l2)
		// 	      jacobian(phi_1_local_eqn,phi_1_local_unknown)
		// 		-= dtestdx(l,k)*dpsidx(l2,k)*W;
		// 	  }

		// 	// The dependence of phi_1 on M
		// 	double div_psi = 0.0; // Get divergence of shape fn for node l2
		// 	for(unsigned k=0; k<DIM; k++)
		// 	  div_psi += dpsidx(l2,k);

		// 	for(unsigned j=0; j<3; j++)
		// 	  {
		// 	    int m_local_unknown = nodal_local_eqn(l2,m_index_micromag(j));
		// 	    if(m_local_unknown >= 0)
		// 	      jacobian(phi_1_local_eqn,m_local_unknown)
		// 		+= div_psi * test(l) * W;// ??ds add jacobian bit...
		// 	  }

		// 	// No dependence of phi_1 on exchange field

		// 	// Phi_1 does not depend on phi
		//       }
		//   } // End of phi_1 Jacobian calculation
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

    // std::cout << "Residuals from this element" << std::endl;
    // std::cout << "phi residual " << residuals[phi_index_micromag()] << std::endl;
    // std::cout << "phi_1 residual " << residuals[phi_1_index_micromag()] << std::endl;
    // std::cout << "M residuals " << residuals[m_index_micromag(0)]<< " "
    // 	      << residuals[m_index_micromag(1)] << " "
    // 	      << residuals[m_index_micromag(2)] << std::endl;
    // std::cout << "exchange residuals " << residuals[exchange_index_micromag(0)]<< " "
    // 	      << residuals[exchange_index_micromag(1)] << " "
    // 	      << residuals[exchange_index_micromag(2)] << std::endl;
    // std::cout << std::endl << std::endl;

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

    	// Output solution vector at local coordinate s
    	Vector<double> interpolated_solution(nvalue,0.0); //??ds generalise the length?
    	interpolated_solution_micromag(s,interpolated_solution);
    	for(unsigned i=0; i<nvalue; i++)
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

}

#endif
