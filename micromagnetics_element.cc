#ifndef OOMPH_MICROMAGNETICS_ELEMENT_CC
#define OOMPH_MICROMAGNETICS_ELEMENT_CC

#include "micromagnetics_element.h"

using namespace oomph;
using namespace MathematicalConstants;
using namespace VectorOps;

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

  template<unsigned DIM>
  const double MicromagEquations<DIM>::DummyBEMControlledEntry = 10;


  //======================================================================
  ///
  //======================================================================
  template<unsigned DIM>
  unsigned MicromagEquations<DIM>::self_test()
  {
    return 0;
  }

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

    // Get coefficients
    const double llg_precess_c = llg_precession_coeff();
    const double llg_damp_c = llg_damping_coeff();
    const double exch_c = exchange_coeff();
    const double magstatic_c = magnetostatic_coeff();

    // Cache time stepper weights needed in Jacobian (??ds write some
    // documentation on this? needs proper maths really...)
    double d_value_evaltime_by_dvalue_np1 = node_pt(0)->time_stepper_pt()->weight(0,0);
    double d_valuederivative_evaltime_by_dvalue_np1 = node_pt(0)->time_stepper_pt()->weight(1,0);

    //======================================================================
    /// Begin loop over the knots (integration points)
    //======================================================================
    for(unsigned ipt=0, nipt = integral_pt()->nweight(); ipt<nipt; ipt++)
      {
        //======================================================================
        /// Calculate/get/interpolate all values for the residual calculations
        //======================================================================
        Vector<double> s(DIM);
        for(unsigned j=0; j<DIM; j++) {s[j] = integral_pt()->knot(ipt,j);}

        // Create interpolator
        MMInterpolator<DIM> intp(this, s);

        double W = integral_pt()->weight(ipt) * intp.j();

        // {
        // // Allocate memory for local quantities and initialise to zero. dphidx
        // // is also H_demag so we need all 3 components.
        // double itp_phi(0.0), itp_phi_1(0.0);
        // Vector<double> itp_x(DIM,0.0), itp_dphidx(3,0.0),
        //   itp_dphi_1dx(3,0.0),itp_m(3,0.0), itp_dmdt(3,0.0);
        // DenseDoubleMatrix itp_dmdx(3,3,0.0);

        // // Interpolate values at knot by looping over nodes adding contributions
        // for(unsigned l=0;l<n_node;l++)
        //   {
        //     itp_phi += nodal_value(l,phi_index_micromag()) * intp.psi(l);
        //     itp_phi_1 += nodal_value(l,phi_1_index_micromag()) * intp.psi(l);

        //     Vector<double> dmdt(3,0.0);
        //     dm_dt_micromag(l,dmdt); // get dmdt at node l

        //     for(unsigned j=0; j<3; j++)
        //       {
        //         itp_dmdt[j] += dmdt[j]*intp.psi(l);
        //         itp_m[j] += nodal_value(l,m_index_micromag(j))*intp.psi(l);
        //       }

        //     // Interpolate spatial values/derivatives
        //     for(unsigned j=0; j<DIM; j++)
        //       {
        //         itp_x[j] += nodal_position(l,j)*intp.psi(l);
        //         itp_dphidx[j] += nodal_value(l,phi_index_micromag())*intp.dpsidx(l,j);
        //         itp_dphi_1dx[j] += nodal_value(l,phi_1_index_micromag())*intp.dpsidx(l,j);
        //         for(unsigned k=0; k<3; k++)
        //           itp_dmdx(k,j) += nodal_value(l,m_index_micromag(k))*intp.dpsidx(l,j);
        //       }
        //   }

        // if(two_norm_diff(intp.m(), itp_m) > 1e-12)
        //   {
        //     std::cout << itp_m << std::endl;
        //     std::cout << intp.m() << std::endl;
        //     std::cout << std::endl;
        //   }

        // if(two_norm_diff(intp.dmdt(), itp_dmdt) > 1e-12)
        //   {
        //     std::cout << itp_dmdt << std::endl;
        //     std::cout << intp.dmdt() << std::endl;
        //     std::cout << std::endl;
        //   }
        // }

        // Calculate other things:

        // Divergence of m
        const double itp_divm = intp.div_m();

        // Source functions (for debugging, normally zero)
        const double phi_source = get_phi_source(intp.time(),intp.x());
        const double phi_1_source = get_phi_1_source(intp.time(),intp.x());

        // Fields at integration point
        const Vector<double> h_applied = get_applied_field(intp.time(), intp.x(), s);
        const Vector<double> h_cryst_anis = get_H_cryst_anis_field(intp.time(),
                                                                   intp.x(), intp.m());

        // Cross product for magnetostatic contribution
        Vector<double> h_magnetostatic(3,0.0);
        for(unsigned j=0; j<3; j++) {h_magnetostatic[j] = -1 * magstatic_c * intp.dphidx()[j];}

        const Vector<double> itp_mxdmdt = cross(intp.m(), intp.dmdt());
        const Vector<double> itp_mxhapp = cross(intp.m(), h_applied);
        const Vector<double> itp_mxhca = cross(intp.m(), h_cryst_anis);
        const Vector<double> mxhms = cross(intp.m(), h_magnetostatic);

        if(two_norm(h_magnetostatic) != 0.0)
          {
            std::cout <<  "non-zero hms from built in phis!" << std::endl;
          }

        //======================================================================
        /// Use the above values to calculate the residuals
        //======================================================================

        // Loop over the test functions/nodes adding contributions to residuals
        for(unsigned l=0;l<n_node;l++)
          {
            // Total potential (phi)
            const int phi_eqn = nodal_local_eqn(l,phi_index_micromag());
            // If value is not pinned and not on the boundary.( We need to treat
            // boundary values of phi differently since they are determined by
            // the boundary element matrix and phi_1.)
            if((phi_eqn >= 0) && (!(node_pt(l)->is_on_boundary())))
              {
                std::cout <<  "unpinned phis!" << std::endl;
                residuals[phi_eqn] -= phi_source*intp.test(l)*W; // source
                residuals[phi_eqn] -= itp_divm*intp.test(l)*W;         // div(m)
                for(unsigned k=0;k<DIM;k++)                       // Poisson
                  residuals[phi_eqn] -= intp.dphidx()[k]*intp.dtestdx(l,k)*W;
              }

            // Reduced potential (phi_1), only difference is in b.c.s
            const int phi_1_eqn = nodal_local_eqn(l,phi_1_index_micromag());
            if(phi_1_eqn >= 0)
              {
                std::cout <<  "unpinned phis!" << std::endl;
                residuals[phi_1_eqn] -= phi_1_source*intp.test(l)*W;
                residuals[phi_1_eqn] -= itp_divm*intp.test(l)*W;
                for(unsigned k=0;k<DIM;k++)
                  residuals[phi_1_eqn] -= intp.dphi1dx()[k]*intp.dtestdx(l,k)*W;
              }

            // LLG itself (m, time evolution)
            //=================================================

            // Exchange residual contribution after integration by parts:
            //??ds possibly could optimise - this is calculated twice
            Vector<double> gradtestdotgradmi(3,0.0);
            for(unsigned i=0; i<3; i++)
              for(unsigned j=0; j<DIM; j++)
                gradtestdotgradmi[i] += intp.dtestdx(l,j) * intp.dmdx()[i][j];

            // Cross product for exchange contribution
            const Vector<double> mxexch = cross(intp.m(), gradtestdotgradmi);

            // std::cout << "mxhapp  = " << itp_mxhapp[0] << " "
            //           << itp_mxhapp[1] <<  " " << itp_mxhapp[2] << "\n";
            // std::cout << "mxhms  = " << mxhms[0] << " "
            //           << mxhms[1]<< " " << mxhms[2] << "\n";
            // std::cout << "mxexch  = " << mxexch[0] << " "
            //           << mxexch[1] << " " << mxexch[2] << std::endl;

            // add to residual
            for(unsigned i=0; i<3; i++)
              {
                const int m_eqn = nodal_local_eqn(l,m_index_micromag(i));
                if(m_eqn >= 0)  // If it's not a boundary condition
                  {
                    // dmdt, mxh_ap, mxh_ca, mxh_ms and mxdmdt terms
                    residuals[m_eqn] += ( intp.dmdt()[i]
                                          + llg_precess_c * itp_mxhapp[i]
                                          + llg_precess_c * itp_mxhca[i]
                                          + llg_precess_c * mxhms[i]
                                          - llg_damp_c *itp_mxdmdt[i] )*intp.test(l)*W;

                    // (m x exchange) term (seperate because it involves
                    // derivatives of the test function).
                    residuals[m_eqn] -= exch_c * llg_precess_c * mxexch[i] * W;
                  }
              }

          } // End of loop over test functions, end of residual calculations


        //======================================================================
        /// If we want the Jacobian as well then calculate it (otherwise
        /// continue to next integration point).
        //======================================================================
        if(!flag) continue;
        // Note that "continue" means go on to the next step of the loop, which
        // is perhaps not entirely clear. We use "if(not thing) then continue"
        // statements rather than "if(thing) then .." because the amount of
        // braces and indentation gets ridiculous with so many nested for loops
        // and if statements.

        // Double loop over nodes for the jacobian
        for(unsigned l=0; l<n_node; l++){

          for(unsigned l2=0;l2<n_node;l2++){

            // Timestepper weight for m at this node, at this
            // time. Assuming all nodes have same timestepper, seems pretty
            // safe bet.
            double mt0weight = d_valuederivative_evaltime_by_dvalue_np1;

            // Pre-calculations
            Vector<double> gradpsil2(3,0.0), gradtestl(3,0.0);
            for(unsigned j=0; j<DIM; j++){
              gradpsil2[j] = intp.dpsidx(l2,j);
              gradtestl[j] = intp.dtestdx(l,j);
            }
            double gradtestldotgradpsil2 = dot(gradtestl, gradpsil2);
            const Vector<double> itpmxgradpsil2 = cross(intp.m(), gradpsil2);

            Vector<double> gradtestdotgradmi(3,0.0);
            for(unsigned i=0; i<3; i++)
              for(unsigned j=0; j<DIM; j++) //??ds repeated calculation..
                gradtestdotgradmi[i] += intp.dtestdx(l,j) * intp.dmdx()[i][j];

            DenseMatrix<double> dhcadm(3,3,0.0);
            get_hca_derivative(intp.time(),intp.x(),intp.m(),intp.psi(l2),dhcadm);

            //=========================================================
            /// Actual Jacobian calculation
            //=========================================================

            // Total potential (phi)
            const int phi_eqn = nodal_local_eqn(l,phi_index_micromag());
            if(phi_eqn >= 0)
              {
                const int phi_unknown = nodal_local_eqn(l2,phi_index_micromag());

                // The diagonal of the phi phi block is non-zero but it is
                // determined by the BEM method. Because of the way FEM assembly
                // is handled we cannot put the values in here. So for now we
                // put in dummy values to reserve space in the sparse matrix.
                if(node_pt(l)->is_on_boundary())
                  {
                    if(phi_eqn == phi_unknown)
                      jacobian(phi_eqn,phi_unknown) = DummyBEMControlledEntry;
                    else
                      jacobian(phi_eqn,phi_unknown) = 0.0;
                  }
                else
                  {
                    // w.r.t. phi
                    if(phi_unknown >= 0)
                      jacobian(phi_eqn,phi_unknown) += -gradtestldotgradpsil2 * W;

                    // w.r.t. m
                    for(unsigned j=0; j<DIM; j++){
                      const int m_unknown = nodal_local_eqn(l2,m_index_micromag(j));
                      if(m_unknown >= 0)
                        jacobian(phi_eqn,m_unknown) += - intp.dpsidx(l2,j) * intp.test(l) * W;
                    }
                  }

                // nothing w.r.t. phi_1 (dependence is in BEM matrix)
              }


            // Reduced potential (phi_1), only difference between this and phi
            // is in b.c.s
            const int phi_1_eqn = nodal_local_eqn(l,phi_1_index_micromag());
            if(phi_1_eqn >= 0){

              // w.r.t. phi_1
              const int phi_1_unknown = nodal_local_eqn(l2,phi_1_index_micromag());
              if(phi_1_unknown >= 0)
                jacobian(phi_1_eqn,phi_1_unknown) += - gradtestldotgradpsil2 * W;

              // w.r.t. m
              for(unsigned j=0; j<DIM; j++){
                const int m_unknown = nodal_local_eqn(l2,m_index_micromag(j));
                if(m_unknown >= 0)
                  jacobian(phi_1_eqn,m_unknown) += - intp.dpsidx(l2,j) * intp.test(l) * W;
              }

              // nothing w.r.t. phi
            }



            //======================================================================
            /// LLg derivatives
            //======================================================================
            Vector<int> m_eqn(3), m_unknown(3);
            for(unsigned j=0; j<3; j++)
              {
                m_eqn[j] = nodal_local_eqn(l,m_index_micromag(j));
                m_unknown[j] = nodal_local_eqn(l2,m_index_micromag(j));
              }

            // w.r.t. phi
            const int phi_unknown = nodal_local_eqn(l2,phi_index_micromag());
            if(phi_unknown >= 0){
              for(unsigned j=0; j<3; j++){
                if(m_eqn[j] >= 0)
                  jacobian(m_eqn[j],phi_unknown) -= llg_precess_c
                    * magstatic_c * itpmxgradpsil2[j] * W * intp.test(l);
              }
            }

            // nothing w.r.t. phi1


            // w.r.t. m
            for(unsigned j=0; j<3; j++) // loop over the m we differentiate by
              {
                if((!(m_unknown[j] >= 0)) || (!(m_eqn[j] >= 0))) continue;

                Vector<double> jhat(3,0.0); jhat[j] = 1.0;

                // exchange pre-calcs
                const Vector<double> mxjhat = cross(intp.m(),jhat);
                const Vector<double> jhatxgradtestdotgradmi = cross(jhat,gradtestdotgradmi);

                // jhat x nondiffterms precalcs
                Vector<double> nondiffterms(3,0.0);
                for(unsigned i=0; i<3; i++)
                  nondiffterms[i] = - llg_damp_c * intp.dmdt()[i]
                    + llg_precess_c * h_magnetostatic[i]
                    + llg_precess_c * (h_cryst_anis[i] + h_applied[i]);
                const Vector<double> jhatxnondiffterms = cross(jhat,nondiffterms);

                // itp_m x terms precalcs
                Vector<double> diffterms(3,0.0);
                for(unsigned i=0; i<3; i++)
                  diffterms[i] += llg_precess_c * dhcadm(i,j);
                diffterms[j] -= llg_damp_c * intp.psi(l2) * mt0weight;
                const Vector<double> mxdiffterms = cross(intp.m(),diffterms);

                // mass matrix component due to time derivative
                jacobian(m_eqn[j],m_unknown[j])
                  += intp.test(l) * intp.psi(l2) * mt0weight * W;

                for(unsigned i=0; i<3; i++) // loop over the m we differentiate w.r.t.
                  {
                    if(!(m_eqn[i] >= 0)) continue;

                    // chain rule gives:
                    // d/dmj( m x (...)) = dmidmj x (....) + m x d/dmj(.....)

                    // dmidmj x (....)
                    jacobian(m_eqn[i],m_unknown[j]) +=
                      W * intp.test(l) * intp.psi(l2) * jhatxnondiffterms[i]
                      * d_value_evaltime_by_dvalue_np1;

                    // m x d/dmj(.....)
                    jacobian(m_eqn[i],m_unknown[j]) +=
                      W * intp.test(l) * mxdiffterms[i];

                    // Exchange contribution
                    jacobian(m_eqn[i],m_unknown[j])
                      -= llg_precess_c * exch_c * W *
                      ( intp.psi(l2) * jhatxgradtestdotgradmi[i]
                        + mxjhat[i] * gradtestldotgradpsil2)
                      * d_value_evaltime_by_dvalue_np1;

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

    // Tecplot header info
    outfile << tecplot_zone_string(n_plot);

    // Loop over plot points
    unsigned num_plot_points=nplot_points(n_plot);
    for (unsigned iplot=0;iplot<num_plot_points;iplot++)
      {
        get_s_plot(iplot,n_plot,s);

        MMInterpolator<DIM> intp(this, s);

        // output eulerian coordinates of plot point
        for(unsigned i=0; i<DIM; i++) outfile << intp.x(i) << " ";

        // std::cout << intp.x() << std::endl;
        // std::cout << intp.x(0) << " " << intp.x(1) << std::endl;

        // Output the magnetostatic field (= - dphidx) at this point
        Vector<double> intp_dphidx = intp.dphidx();
        for(unsigned i=0; i<3; i++) outfile << -intp_dphidx[i] << " ";

        // Phi 1 at this point
        outfile << intp.phi() << " ";

        // Phi total at this point
        outfile << intp.phi1() << " ";

        // Output m at this point
        for(unsigned i=0; i<3; i++) outfile << intp.m()[i] << " ";

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
  //                                  const unsigned &nplot)
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
  //     {
  //       fprintf(file_pt,"%g ",interpolated_x(s,i));
  //     }
  //       fprintf(file_pt,"%g \n",itp_phi_micromag(s));
  //     }

  //   // Write tecplot footer (e.g. FE connectivity lists)
  //   write_tecplot_zone_footer(file_pt,nplot);
  // }

}

#endif
