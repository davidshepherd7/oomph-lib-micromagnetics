
#include "micromagnetics_element.h"
#include "micromagnetics_flux_element.h"

// For output in
#include <iomanip>

using namespace oomph;
using namespace MathematicalConstants;
using namespace VectorOps;

namespace oomph
{

  /// ??ds Magic number, probably bad...
  const double MicromagEquations::DummyBEMControlledEntry = 10;

  /// \short Integrate a function given by func_pt over the element using
  /// the given integral_pt(). Because C++ sucks we have to do this with
  /// weird function objects.
  double MicromagEquations::integrate_over_element(const ElementalFunction* func_pt) const
  {
    double result = 0;

    // Loop over knots and sum
    for(unsigned ipt=0, nipt = integral_pt()->nweight(); ipt<nipt; ipt++)
      {
        // Get position in element
        Vector<double> s(this->dim());
        for(unsigned j=0; j<this->dim(); j++)
          {s[j] = integral_pt()->knot(ipt,j);}

        MMInterpolator intp(this, s);
        double J = intp.j();
        double w = integral_pt()->weight(ipt);

        // Add contribution
        result += func_pt->call(this, &intp) * w * J;
      }

    return result;
  }

  void MicromagEquations::check_interpolation(MMInterpolator& intp) const
  {
#ifdef PARANOID
    // Track whether any errors have been found
    bool failed = false;

    //??ds move to an option somewhere higher up?
    bool throw_on_fail = false;

    // Tolerence for differences in interpolation
    double tol = 5e-10;

    // Check that we are not using midpoint: old way doesn't work for midpoint
    if(node_pt(0)->time_stepper_pt()->nprev_values_for_value_at_evaluation_time() != 1)
      {
        return;
      }

    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Allocate memory for local quantities and initialise to zero. dphidx
    // is also H_demag so we need all 3 components.
    double itp_phi(0.0), itp_phi_1(0.0);
    Vector<double> itp_x(nodal_dimension(),0.0), itp_dphidx(3,0.0),
      itp_dphi_1dx(3,0.0),itp_m(3,0.0), itp_dmdt(3,0.0);
    DenseDoubleMatrix itp_dmdx(3,3,0.0);

    // Interpolate values at knot by looping over nodes adding contributions
    for(unsigned l=0;l<n_node;l++)
      {
        itp_phi += nodal_value(l,phi_index_micromag()) * intp.psi(l);
        itp_phi_1 += nodal_value(l,phi_1_index_micromag()) * intp.psi(l);

        Vector<double> dmdt(3,0.0);
        dm_dt_micromag(l,dmdt); // get dmdt at node l

        for(unsigned j=0; j<3; j++)
          {
            itp_dmdt[j] += dmdt[j]*intp.psi(l);
            itp_m[j] += nodal_value(l,m_index_micromag(j))*intp.psi(l);
          }

        // Interpolate spatial values/derivatives
        for(unsigned j=0; j<nodal_dimension(); j++)
          {
            itp_x[j] += nodal_position(l,j)*intp.psi(l);
            itp_dphidx[j] += nodal_value(l,phi_index_micromag())*intp.dpsidx(l,j);
            itp_dphi_1dx[j] += nodal_value(l,phi_1_index_micromag())*intp.dpsidx(l,j);
            for(unsigned k=0; k<3; k++)
              itp_dmdx(k,j) += nodal_value(l,m_index_micromag(k))*intp.dpsidx(l,j);
          }
      }

    std::cerr << std::setprecision(12);

    if(two_norm_diff(intp.x(), itp_x) > tol)
      {
        std::cerr << "Error in x interpolation (order is: twonorm, old, new)"
                  << two_norm_diff(intp.x(), itp_x) << std::endl;
        std::cerr << itp_x << std::endl;
        std::cerr << intp.x() << std::endl;
        std::cerr << std::endl;
        failed = true;
      }

    if(two_norm_diff(intp.m(), itp_m) > tol)
      {
        std::cerr << "Error in m interpolation (order is: twonorm, old, new)"
                  << two_norm_diff(intp.m(), itp_m) << std::endl;
        std::cerr << itp_m << std::endl;
        std::cerr << intp.m() << std::endl;
        std::cerr << std::endl;
        failed = true;
      }

    if(two_norm_diff(intp.dmdt(), itp_dmdt) > tol)
      {
        std::cerr << "Error in dmdt interpolation (order is: twonorm, old, new)"
                  << two_norm_diff(intp.dmdt(), itp_dmdt) << std::endl;
        std::cerr << itp_dmdt << std::endl;
        std::cerr << intp.dmdt() << std::endl;
        std::cerr << std::endl;
        failed = true;
      }

    //??ds dmdx?

    //??ds phis?

    if(failed && throw_on_fail)
      {
        std::string error_msg = "Error: difference between old interpolation and new";
        throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
  }


  void MicromagEquations::check_interpolation2(const Vector<double> &s) const
  {
#ifdef PARANOID
    // Track whether any errors have been found
    bool failed = false;

    //??ds move to an option somewhere higher up?
    bool throw_on_fail = false;

    // Tolerence for differences in interpolation
    double tol = 5e-10;

    MMInterpolator intp(this, s);
    MMArrayInterpolator<5> a_intp(this, s);

    std::cerr << std::setprecision(12);

    // Get arrays as vectors for easy comparison
    Vector<double> itp_x, itp_m, itp_dmdt;
    itp_x.assign(a_intp.x(), a_intp.x() + nodal_dimension());
    itp_m.assign(a_intp.m(), a_intp.m() + 3);
    itp_dmdt.assign(a_intp.dmdt(), a_intp.dmdt() + 3);

    if(two_norm_diff(intp.x(), itp_x) > tol)
      {
        std::cerr << "Error in x interpolation (order is: twonorm, old, new)"
                  << two_norm_diff(intp.x(), itp_x) << std::endl;
        std::cerr << itp_x << std::endl;
        std::cerr << intp.x() << std::endl;
        std::cerr << std::endl;
        failed = true;
      }

    if(two_norm_diff(intp.m(), itp_m) > tol)
      {
        std::cerr << "Error in m interpolation (order is: twonorm, old, new)"
                  << two_norm_diff(intp.m(), itp_m) << std::endl;
        std::cerr << itp_m << std::endl;
        std::cerr << intp.m() << std::endl;
        std::cerr << std::endl;
        failed = true;
      }

    if(two_norm_diff(intp.dmdt(), itp_dmdt) > tol)
      {
        std::cerr << "Error in dmdt interpolation (order is: twonorm, old, new)"
                  << two_norm_diff(intp.dmdt(), itp_dmdt) << std::endl;
        std::cerr << itp_dmdt << std::endl;
        std::cerr << intp.dmdt() << std::endl;
        std::cerr << std::endl;
        failed = true;
      }

    //??ds dmdx?

    //??ds phis?

    if(failed && throw_on_fail)
      {
        std::string error_msg = "Error: difference between old interpolation and new";
        throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
  }

  //======================================================================
  /// Compute element residual Vector and/or element Jacobian matrix
  ///
  /// flag=1: compute both
  /// flag=0: compute only residual Vector
  ///
  /// Pure version without hanging nodes
  //======================================================================
  void MicromagEquations::fill_in_generic_residual_contribution_micromag
  (Vector<double> &residuals, DenseMatrix<double> &jacobian,
   const unsigned& flag) const
  {

    //======================================================================
    /// Get some useful numbers and set up memory.
    //======================================================================

    // Find out how many nodes there are
    const unsigned n_node = nnode();
    const unsigned ndim = nodal_dimension();
    const unsigned eldim = dim();

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
        Vector<double> s(eldim);
        for(unsigned j=0; j<eldim; j++) {s[j] = integral_pt()->knot(ipt,j);}

        // Create interpolator //??ds maybe should move this out of loop,
        // add .new_point() function or something?
        MMArrayInterpolator<5> intp(this, s);

        double W = integral_pt()->weight(ipt) * intp.j();

        // Calculate other things:

        // Divergence of m
        const double itp_divm = intp.div_m();

        Vector<double> xvec(ndim, 0.0);
        for(unsigned i=0; i<ndim; i++) {xvec[i] = intp.x()[i];}

        Vector<double> mvec(3, 0.0);
        for(unsigned i=0; i<3; i++) {mvec[i] = intp.m()[i];}

        // Source functions (for debugging, normally zero)
        const double phi_source = get_phi_source(intp.time(),xvec);
        const double phi_1_source = get_phi_1_source(intp.time(),xvec);

        // Fields at integration point
        Vector<double> h_applied, h_cryst_anis, h_magnetostatic;
        get_applied_field(intp.time(), xvec, s, h_applied);
        get_H_cryst_anis_field(intp.time(), xvec, mvec, h_cryst_anis);
        get_magnetostatic_field(&intp, h_magnetostatic);


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
                std::cerr <<  "unpinned phis!" << std::endl;
                residuals[phi_eqn] -= phi_source*intp.test(l)*W; // source
                residuals[phi_eqn] -= itp_divm*intp.test(l)*W;         // div(m)
                for(unsigned k=0;k<ndim;k++)                       // Poisson
                  residuals[phi_eqn] -= intp.dphidx()[k]*intp.dtestdx(l,k)*W;
              }

            // Reduced potential (phi_1), only difference is in b.c.s
            const int phi_1_eqn = nodal_local_eqn(l,phi_1_index_micromag());
            if(phi_1_eqn >= 0)
              {
                std::cerr <<  "unpinned phis!" << std::endl;
                residuals[phi_1_eqn] -= phi_1_source*intp.test(l)*W;
                residuals[phi_1_eqn] -= itp_divm*intp.test(l)*W;
                for(unsigned k=0;k<ndim;k++)
                  residuals[phi_1_eqn] -= intp.dphi1dx()[k]*intp.dtestdx(l,k)*W;
              }

            // LLG itself (m, time evolution)
            //=================================================

            // Exchange residual contribution after integration by parts:
            //??ds possibly could optimise - this is calculated twice
            double gradtestdotgradmi[3] = {0,0,0};
            for(unsigned i=0; i<3; i++)
              for(unsigned j=0; j<ndim; j++)
                gradtestdotgradmi[i] += intp.dtestdx(l,j) * intp.dmdx(i)[j];

            // add to residual
            for(unsigned i=0; i<3; i++)
              {
                const int m_eqn = nodal_local_eqn(l, m_index_micromag(i));
                if(m_eqn >= 0)  // If it's not a boundary condition
                  {
                    // dmdt, mxh_ap, mxh_ca, mxh_ms and mxdmdt terms
                    residuals[m_eqn] +=
                      ( intp.dmdt()[i]
                        + llg_precess_c * opt_cross(i, intp.m(), h_applied)
                        + llg_precess_c * opt_cross(i, intp.m(), h_cryst_anis)
                        + llg_precess_c * opt_cross(i, intp.m(), h_magnetostatic)
                        - llg_damp_c * opt_cross(i, intp.m(), intp.dmdt())
                        )*intp.test(l)*W;

                    // (m x exchange) term (separate because it involves
                    // derivatives of the test function).
                    residuals[m_eqn] -= exch_c * llg_precess_c *
                      opt_cross(i, intp.m(), gradtestdotgradmi) * W;
                  }
                else
                  {
                    std::cerr << "PINNED magnetisation! shouldn't be any of these!"
                              <<std::endl;
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


        // nondiffterms precalcs
        double nondiffterms[3];
        for(unsigned i=0; i<3; i++)
          nondiffterms[i] = - llg_damp_c * intp.dmdt()[i]
            + llg_precess_c * h_magnetostatic[i]
            + llg_precess_c * (h_cryst_anis[i] + h_applied[i]);


        double gradpsil2[3];
        double gradtestl[3];
        double gradtestdotgradmi[3];
        double diffterms[3];
        int m_eqn[3], m_unknown[3];
        double dhcadm[3][3];



        // Double loop over nodes for the jacobian
        for(unsigned l=0; l<n_node; l++){

          // Pre calculate values with no l2 dependencies
          for(unsigned j=0; j<ndim; j++) {gradtestl[j] = intp.dtestdx(l,j);}
          for(unsigned i=0; i<3; i++)
            {
              m_eqn[i] = nodal_local_eqn(l,m_index_micromag(i));

              gradtestdotgradmi[i] = 0.0;
              for(unsigned j=0; j<ndim; j++) //??ds repeated calculation..
                {
                  gradtestdotgradmi[i] += intp.dtestdx(l,j) * intp.dmdx(i)[j];
                }
            }

          for(unsigned l2=0;l2<n_node;l2++){

            // Timestepper weight for m at this node, at this
            // time. Assuming all nodes have same timestepper, seems pretty
            // safe bet.
            double mt0weight = d_valuederivative_evaltime_by_dvalue_np1;

            for(unsigned j=0; j<ndim; j++) {gradpsil2[j] = intp.dpsidx(l2,j);}

            double gradtestldotgradpsil2 = 0.0;
            for(unsigned i=0; i < ndim; i++)
              {
                gradtestldotgradpsil2 += gradtestl[i] * gradpsil2[i];
              }

            get_hca_derivative(intp.time(),xvec,mvec,intp.psi(l2),dhcadm);

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
                    for(unsigned j=0; j<ndim; j++){
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
              for(unsigned j=0; j<ndim; j++){
                const int m_unknown = nodal_local_eqn(l2,m_index_micromag(j));
                if(m_unknown >= 0)
                  jacobian(phi_1_eqn,m_unknown) += - intp.dpsidx(l2,j) * intp.test(l) * W;
              }

              // nothing w.r.t. phi
            }



            //======================================================================
            /// LLg derivatives
            //======================================================================
            for(unsigned j=0; j<3; j++)
              {
                m_unknown[j] = nodal_local_eqn(l2,m_index_micromag(j));
              }

            // w.r.t. phi
            const int phi_unknown = nodal_local_eqn(l2,phi_index_micromag());
            if(phi_unknown >= 0){
              for(unsigned j=0; j<3; j++){
                if(m_eqn[j] >= 0)
                  jacobian(m_eqn[j],phi_unknown) -= llg_precess_c
                    * magstatic_c * opt_cross(j, intp.m(), gradpsil2) * W * intp.test(l);
              }
            }

            // nothing w.r.t. phi1


            // w.r.t. m
            for(unsigned j=0; j<3; j++) // loop over the m we differentiate by
              {
                if((!(m_unknown[j] >= 0)) || (!(m_eqn[j] >= 0))) continue;

                // Initialise jhat
                double jhat[3] = {0,0,0}; jhat[j] = 1.0;

                // itp_m x terms precalcs
                for(unsigned i=0; i<3; i++)
                  {
                    diffterms[i] = llg_precess_c * dhcadm[i][j];
                  }
                diffterms[j] -= llg_damp_c * intp.psi(l2) * mt0weight;

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
                      W * intp.test(l) * intp.psi(l2) * opt_cross(i, jhat, nondiffterms)
                      * d_value_evaltime_by_dvalue_np1;

                    // m x d/dmj(.....)
                    jacobian(m_eqn[i], m_unknown[j]) +=
                      W * intp.test(l) * opt_cross(i, intp.m(), diffterms);

                    // Exchange contribution
                    jacobian(m_eqn[i],m_unknown[j])
                      -= llg_precess_c * exch_c * W *
                      ( intp.psi(l2) * opt_cross(i, jhat, gradtestdotgradmi)

                        + opt_cross(i, intp.m(), jhat) * gradtestldotgradpsil2
                        )
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
  void MicromagEquations::output(std::ostream &outfile,
                                 const unsigned &n_plot)
  {
    //Vector of local coordinates
    Vector<double> s(nodal_dimension());

    // Tecplot header info
    outfile << tecplot_zone_string(n_plot);

    // Loop over plot points
    unsigned num_plot_points=nplot_points(n_plot);
    for (unsigned iplot=0;iplot<num_plot_points;iplot++)
      {
        get_s_plot(iplot,n_plot,s);

        MMInterpolator intp(this, s);

        // output eulerian coordinates of plot point
        for(unsigned i=0; i<nodal_dimension(); i++) outfile << intp.x(i) << " ";

        // Output the magnetostatic field (= - dphidx) at this point
        Vector<double> intp_dphidx = intp.dphidx();
        intp_dphidx.resize(3, 0.0); // make sure it has a 3rd entry
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
  void MicromagEquations::
  output_fct(std::ostream &outfile, const unsigned &n_plot,
             const double& time,
             const FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
  {

    //Vector of local coordinates
    Vector<double> s(nodal_dimension());

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
        Vector<double> x(nodal_dimension(),0.0);
        for(unsigned i=0; i<nodal_dimension(); i++)
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

  void MicromagEquations::
  get_magnetostatic_field(const Vector<double> &s,
                          Vector<double> &h_magnetostatic) const
  {
    // Construct an interpolator and call the underlying function.
    MMArrayInterpolator<5> intp(this, s);
    get_magnetostatic_field(&intp, h_magnetostatic);
  }

  void MicromagEquations::
  get_magnetostatic_field(MMArrayInterpolator<5>* intp_pt,
                          Vector<double> &h_magnetostatic) const
  {
    #ifdef PARANOID
    if(intp_pt == 0)
      {
        std::string error_msg = "Null interpolator!";
        throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

    const double magstatic_c = magnetostatic_coeff();

    h_magnetostatic.assign(3, 0.0);
    for(unsigned j=0; j<nodal_dimension(); j++)
      {
        h_magnetostatic[j] = -1 * magstatic_c * intp_pt->dphidx()[j];
      }
  }

  /// For micromagnetics the source function is the divergence of the
  /// magnetisation.
  void SemiImplicitMicromagEquations::
  get_magnetostatic_field(MMArrayInterpolator<5>* intp_pt,
                          Vector<double> &h_ms) const
  {
    // Lots of checks because this is stupid really...
#ifdef PARANOID
    if(this->nnode() != magnetostatic_field_element_pt()->nnode())
      {
        std::ostringstream error_msg;
        error_msg << "Elements must be the same geometry for this to "
                  << "work... sorry for the hackyness. Maybe you can fix it.";
        throw OomphLibError(error_msg.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

    if(this->dim() != magnetostatic_field_element_pt()->dim())
      {
        std::ostringstream error_msg;
        error_msg << "Elements must be the same geometry for this to "
                  << "work... sorry for the hackyness. Maybe you can fix it.";
        throw OomphLibError(error_msg.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

    if(this->integral_pt() != magnetostatic_field_element_pt()->integral_pt())
      {
        std::ostringstream error_msg;
        error_msg << "Elements must have the same integration scheme for this to"
                  << "work... sorry for the hackyness. Maybe you can fix it.";
        throw OomphLibError(error_msg.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

    // Get magnetostatic field from field element
    magnetostatic_field_element_pt()->magnetostatic_field(intp_pt->s(), h_ms);
  }

  /// For micromagnetics the source function is the divergence of the
  /// magnetisation.
  void SemiImplicitMicromagEquations::
  get_magnetostatic_field_time_derivative(const Vector<double> &s,
                                          Vector<double> &h_ms) const
  {
    // Lots of checks because this is stupid really...
#ifdef PARANOID
    if(this->nnode() != magnetostatic_field_element_pt()->nnode())
      {
        std::ostringstream error_msg;
        error_msg << "Elements must be the same geometry for this to "
                  << "work... sorry for the hackyness. Maybe you can fix it.";
        throw OomphLibError(error_msg.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

    if(this->dim() != magnetostatic_field_element_pt()->dim())
      {
        std::ostringstream error_msg;
        error_msg << "Elements must be the same geometry for this to "
                  << "work... sorry for the hackyness. Maybe you can fix it.";
        throw OomphLibError(error_msg.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

    if(this->integral_pt() != magnetostatic_field_element_pt()->integral_pt())
      {
        std::ostringstream error_msg;
        error_msg << "Elements must have the same integration scheme for this to"
                  << "work... sorry for the hackyness. Maybe you can fix it.";
        throw OomphLibError(error_msg.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

    // Get magnetostatic field from field element
    magnetostatic_field_element_pt()->magnetostatic_field_time_derivative(s, h_ms);
  }


  //======================================================================
  /// Validate computed M against exact solution.
  ///
  /// Solution is provided via function pointer.
  /// Plot error at integration points.
  ///
  //======================================================================
  void MicromagEquations::
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
    //??ds use n_intpt/nodal_dimension()? - only good for simple shaped elements?
    //??ds causes small issuse with output - points for soln/exact are at the corners of elements but for errors they are at the int pts

    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
      {
        // Get s (local coordinate)
        Vector<double> s(nodal_dimension(),0.0);
        for(unsigned i=0; i<nodal_dimension(); i++) {s[i] = integral_pt()->knot(ipt,i);}

        // Get x (global coordinate) and output
        Vector<double> x(nodal_dimension(),0.0);
        interpolated_x(s,x);
        for(unsigned i=0; i<nodal_dimension(); i++){outfile << x[i] << " ";}

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



  template <unsigned DIM, unsigned NNODE_1D>
  void QMicromagElement<DIM, NNODE_1D>::fill_in_face_element_contribution_to_jacobian
  (DenseMatrix<double> &jacobian) const
  {
    std::set<FiniteElement*>::iterator it;
    for(it=this->Face_element_pts.begin(); it!=this->Face_element_pts.end(); it++)
      {
        MicromagFluxElement<QMicromagElement<DIM,NNODE_1D> >* flux_ele_pt =
          dynamic_cast<MicromagFluxElement<QMicromagElement<DIM,NNODE_1D> >* >
          (*it);
        flux_ele_pt->fill_in_bulk_contribution_to_face_jacobian(jacobian);
      }
  }


  template < unsigned DIM, unsigned NNODE_1D>
  void TMicromagElement<DIM, NNODE_1D>::fill_in_face_element_contribution_to_jacobian
  (DenseMatrix<double> &jacobian) const
  {
    std::set<FiniteElement*>::iterator it;
    for(it=this->Face_element_pts.begin(); it!=this->Face_element_pts.end(); it++)
      {
        MicromagFluxElement<TMicromagElement<DIM,NNODE_1D> >* flux_ele_pt =
          dynamic_cast<MicromagFluxElement<TMicromagElement<DIM,NNODE_1D> >* >
          (*it);
        flux_ele_pt->fill_in_bulk_contribution_to_face_jacobian(jacobian);
      }
  }

  template < unsigned DIM, unsigned NNODE_1D>
  void QSemiImplicitMicromagElement<DIM, NNODE_1D>::fill_in_face_element_contribution_to_jacobian
  (DenseMatrix<double> &jacobian) const
  {
    std::set<FiniteElement*>::iterator it;
    for(it=this->Face_element_pts.begin(); it!=this->Face_element_pts.end(); it++)
      {
        MicromagFluxElement<QSemiImplicitMicromagElement<DIM,NNODE_1D> >* flux_ele_pt =
          dynamic_cast<MicromagFluxElement<QSemiImplicitMicromagElement<DIM,NNODE_1D> >* >
          (*it);
        flux_ele_pt->fill_in_bulk_contribution_to_face_jacobian(jacobian);
      }
  }

  template < unsigned DIM, unsigned NNODE_1D>
  void TSemiImplicitMicromagElement<DIM, NNODE_1D>::fill_in_face_element_contribution_to_jacobian
  (DenseMatrix<double> &jacobian) const
  {
    std::set<FiniteElement*>::iterator it;
    for(it=this->Face_element_pts.begin(); it!=this->Face_element_pts.end(); it++)
      {
        MicromagFluxElement<TSemiImplicitMicromagElement<DIM,NNODE_1D> >* flux_ele_pt =
          dynamic_cast<MicromagFluxElement<TSemiImplicitMicromagElement<DIM,NNODE_1D> >* >
          (*it);
        flux_ele_pt->fill_in_bulk_contribution_to_face_jacobian(jacobian);
      }
  }

  //====================================================================
  // Force building of templates
  //====================================================================
  template class QMicromagElement<2,2>;
  template class QMicromagElement<2,3>;
  template class QMicromagElement<3,2>;
  template class QMicromagElement<3,3>;

  template class TMicromagElement<2,2>;
  template class TMicromagElement<2,3>;
  template class TMicromagElement<3,2>;
  template class TMicromagElement<3,3>;

  template class QSemiImplicitMicromagElement<2,2>;
  template class QSemiImplicitMicromagElement<2,3>;
  template class QSemiImplicitMicromagElement<3,2>;
  template class QSemiImplicitMicromagElement<3,3>;

  template class TSemiImplicitMicromagElement<2,2>;
  template class TSemiImplicitMicromagElement<2,3>;
  template class TSemiImplicitMicromagElement<3,2>;
  template class TSemiImplicitMicromagElement<3,3>;
}
