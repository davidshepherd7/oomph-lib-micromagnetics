
#include "./micromagnetics_element.h"

using namespace oomph;

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
  void MicromagEquations<DIM>::fill_in_generic_residual_contribution_micromag(Vector<double> &residuals, DenseMatrix<double> &jacobian, const unsigned& flag)
  {

    //Find out how many nodes there are
    const unsigned n_node = nnode();

    //Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);

    // Set up vector to store the local coordinates
    Vector<double> s(DIM);

    // Get current time
    double time = time_pt()->time();

    //Index at which the poisson unknown is stored
    const unsigned phi_nodal_index = phi_index_micromag();

    //Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    //Integers to store the local equation and unknown numbers
    int phi_local_eqn=0, m_local_eqn=0;

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
	// dphidx is also H_demag so we need all 3 components initialised - even if they are zero.
	double interpolated_phi=0.0,  llg_damping_coeff=0.0, llg_precession_coeff=0.0, div_m=0.0;
	Vector<double> interpolated_x(DIM,0.0), interpolated_dphidx(3,0.0);
	Vector<double> H_total(3,0.0), H_cryst_anis(3,0.0), H_applied(3,0.0);
	Vector<double> interpolated_m(3,0.0), interpolated_mxH(3,0.0), interpolated_mxmxH(3,0.0);
	Vector<double> dmdt(3,0.0), interpolated_dmdt(3,0.0);
	Vector<double> interpolated_H_exchange(3,0.0);


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
	    dm_dt_micromag(l,dmdt);
	    for(unsigned j=0; j<3; j++)
	      {
		interpolated_dmdt[j] += dmdt[j]*psi(l);
		interpolated_m[j] += nodal_value(l,M_index_micromag(j))*psi(l);
		interpolated_H_exchange[j] +=nodal_value(l,exchange_index_micromag(j))*psi(l);
	      }

	    // Loop over spatial directions
	    for(unsigned j=0;j<DIM;j++)
	      {
		interpolated_x[j] += nodal_position(l,j)*psi(l);
		interpolated_dphidx[j] += phi_value*dpsidx(l,j);
		div_m += nodal_value(l,M_index_micromag(j)) *dpsidx(l,j);
	      }

	  }


	// Poisson section (demagnetising field calculation)
	//----------------------------------------------------

	// Get source function
	double poisson_source = 0;
	get_poisson_source(time,ipt,interpolated_x,poisson_source);

	// Loop over the test functions
	for(unsigned l=0;l<n_node;l++)
	  {

	    // Get the local equation number for the poisson part
	    phi_local_eqn = nodal_local_eqn(l,phi_nodal_index);

	    if(phi_local_eqn >= 0)	  // If it's not a boundary condition:
	      {
		// Add source term
		residuals[phi_local_eqn] -= poisson_source*test(l)*W; //??ds not sure on sign

		// standard residuals term (not "integrated by parts")
		// Seems to be unstable
		//residuals[phi_local_eqn] -= 4.0*MathematicalConstants::Pi*div_m;

		// The Poisson and divergence of M bits (after reducing the order of differentiation using integration by parts)
		for(unsigned k=0;k<DIM;k++)
		  {
		    residuals[phi_local_eqn] -= interpolated_dphidx[k]*dtestdx(l,k)*W;

		    // this rearrangement using "integration by parts" switches the sign,
		    // requires M_x = 0 at 0 and 1.
		    residuals[phi_local_eqn] += 4.0*MathematicalConstants::Pi
		      *interpolated_m[k]*dtestdx(l,k)*W;

		  }


		// ??ds add in jacobian calculation eventually

	      }

	  } // End of Poisson section]


	// Exchange field section
	//----------------------------------------------------
	// H_ex = coeff * laplacian(M)
	// Weak form: int( H_ex_i - coeff* grad(M).grad(test) ) = 0
	// ??ds only when grad(M).n = 0 at boundaries, otherwise need another term!
	double exchange_coeff = get_exchange_coeff(time, interpolated_x);

	for (unsigned l=0; l<n_node; l++)
	  {
	    // loop over field/magnetisation directions
	    for(unsigned k=0; k<3; k++)
	      {
		// Get the local equation number for the kth componenet of H_ex on node l
		int exchange_local_eqn = nodal_local_eqn(l,exchange_index_micromag(k));

		if(exchange_local_eqn >= 0)  // If it's not a boundary condition
		  {
		    residuals[exchange_local_eqn] += interpolated_H_exchange[k]*test(l)*W;

		    // add coeff * grad(M_k) * grad(test) (then *weight, as always)
		    // we only loop over DIM directions since derrivatives are automatically zero in non-spatial directions.
		    // This rearrangement requires dM/dn = 0 at boundaries.
		    for(unsigned i=0; i<DIM; i++)
		      {
			residuals[exchange_local_eqn] +=
			  exchange_coeff * (interpolated_m[k]*dpsidx(l,i)) * dtestdx(l,i) * W;
		      }
		    //??ds put in jacobian calculation eventually
		  }
	      }
	  } // End of loop over test functions


	// LLG section (time evolution of magnetisation)
	//----------------------------------------------------

	//??ds No boundary conditions on M, for now...

	// Get applied field at this position
	// get_applied_field(time, interpolated_x, H_applied);

	// Get crystalline anisotropy effective field
	// get_H_cryst_anis_field(time, interpolated_x, interpolated_m,  H_cryst_anis);

	// Get LLG source function
	Vector<double> llg_source(3,0.0);
	get_source_llg(time, interpolated_x, llg_source);

	//(-1*interpolated_dphidx is exactly the demagnetising/magnetostatic field: H_demag = - grad(phi))
	Vector<double> H_magnetostatic(3,0.0);
	for(unsigned j=0; j<DIM; j++)
	  {
	    H_magnetostatic[j] = -interpolated_dphidx[j];
	  }

	// Take total of all fields used
	// ??ds pass this entire section out to a function eventually if possible?
	// ??ds add 0.1 to push off maximum (i.e. thermal-ish...)
	for(unsigned j=0; j<3; j++)
	  {
	    H_total[j] = H_applied[j] + H_magnetostatic[j]
	      + interpolated_H_exchange[j] + H_cryst_anis[j];
	  }

	// Get the coefficients for the LLG equation (damping could be a function of position if saturation magnetisation varies)
	llg_damping_coeff = get_llg_damp(time, interpolated_x);
	llg_precession_coeff = get_llg_precess(time, interpolated_x);

	// Get the cross products for the LLG equation
	cross(interpolated_m, H_total, interpolated_mxH);
	cross(interpolated_m, interpolated_mxH, interpolated_mxmxH);

	for (unsigned l=0; l<n_node; l++)
	  {

	    // Calculate residuals for the time evolution equations (Landau-Lifschitz-Gilbert):
	    // dM/dt + gamma/(1+alpha^2) [ (M x H) + (gamma/|M_s|)(M x (M x H)) ] - llg_source = 0

	    // loop over M directions
	    for(unsigned k=0; k<3; k++)
	      {
		// Get the local equation number for the kth component of M part
		m_local_eqn = nodal_local_eqn(l,M_index_micromag(k));

		if(m_local_eqn >= 0)  // If it's not a boundary condition
		  {
		    residuals[m_local_eqn] +=
		      ( interpolated_dmdt[k]
			+ llg_precession_coeff*interpolated_mxH[k]
			+ llg_damping_coeff*interpolated_mxmxH[k]
			- llg_source[k]
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
    double t = time_pt()->time();

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

	// Output phi value at position
	outfile << interpolated_phi_micromag(s) << " ";

	// Output all M values at position
	Vector<double> interpolated_m(3,0.0);
	interpolated_m_micromag(s,interpolated_m);
	for(unsigned i=0; i<3; i++)
	  {
	    outfile << interpolated_m[i] << " ";
	  }

	// Output LLg source (just in case)
	Vector<double> llg_source(3,0.0);
	get_source_llg(t,x,llg_source);
	for(unsigned i=0; i<3; i++)
	  {
	    outfile << llg_source[i] << " ";
	  }

	// End the line ready for next point
	outfile << std::endl;
      }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile,n_plot);

  }

  /// Output a time-dependent exact solution over the element.
  template <unsigned DIM>
  void MicromagEquations<DIM>::output_fct(std::ostream &outfile, const unsigned &n_plot,
					  const double& time, FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt) const
  {

    //Vector of local coordinates
    Vector<double> s(DIM);

    // Get time
    double t = time_pt()->time();

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
	Vector<double> exact_solution(4,0.0);
	(*exact_soln_pt)(t,x,exact_solution);
	for(unsigned i=0; i<4; i++)
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
    outfile << "ZONE" << std::endl;

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

	// Get finite element interpolated magnetisation
	//??ds at the current time
	Vector<double> interpolated_m(3,0.0);
	interpolated_m_micromag(s,interpolated_m);

	// Combine with phi solution into one soln vector "interpolated_soln"
	Vector<double> interpolated_soln(nvalues,0.0);
	interpolated_soln[phi_index_micromag()] = interpolated_phi_micromag(s);
	for(unsigned i=0; i<3; i++)
	  { interpolated_soln[M_index_micromag(i)] = interpolated_m[i]; }
	//??ds include exchange field and anything else I add eventually

	// Get exact solution at point x and time "time"
	Vector<double> exact_soln(nvalues,0.0);
	(*exact_soln_pt)(time,x,exact_soln);

	// Output error then end the line
	// ============================================================
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
