#include "./residual_calculator.h"
#include "./micromagnetics_element.h"
#include "./vector_helpers.h"

using namespace oomph;
using namespace VectorOps;

namespace oomph
{

  void LLResidualCalculator::fill_in_generic_residual_contribution
  (const MicromagEquations* const e_pt,
   Vector<double> &residuals, DenseMatrix<double> &jacobian,
   const unsigned& flag) const
  {
    //======================================================================
    /// Get some useful numbers and set up memory.
    //======================================================================

    // Find out how many nodes there are
    const unsigned n_node = e_pt->nnode();
    const unsigned ndim = e_pt->nodal_dimension();
    const unsigned eldim = e_pt->dim();

    // Get coefficients
    const double llg_damp_c = e_pt->llg_damping_coeff();
#ifdef PARANOID
    if(llg_damp_c > 0.5)
      {
        throw OomphLibError("Can't handle high damping in Landau-Lifshitz yet",
                            OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif
    const double ll_damp_c = -(std::sqrt(-4*llg_damp_c*llg_damp_c + 1) - 1)
      /(2*llg_damp_c);

    // Cache time stepper weights needed in Jacobian (??ds write some
    // documentation on this? needs proper maths really...)
    double d_value_evaltime_by_dvalue_np1 =
      e_pt->node_pt(0)->time_stepper_pt()->weight(0,0);
    double d_valuederivative_evaltime_by_dvalue_np1 =
      e_pt->node_pt(0)->time_stepper_pt()->weight(1,0);


    //======================================================================
    /// Begin loop over the knots (integration points)
    //======================================================================
    for(unsigned ipt=0, nipt = e_pt->integral_pt()->nweight(); ipt<nipt; ipt++)
      {
        //======================================================================
        /// Calculate/get/interpolate all values for the residual calculations
        //======================================================================
        Vector<double> s(eldim);
        for(unsigned j=0; j<eldim; j++) {s[j] = e_pt->integral_pt()->knot(ipt,j);}

        // Create interpolator //??ds maybe should move this out of loop,
        // add .new_point() function or something?
        MMArrayInterpolator<5> intp(e_pt, s);

        double W = e_pt->integral_pt()->weight(ipt) * intp.j();

        // Calculate other things:

        // Divergence of m
        const double itp_divm = intp.div_m();

        Vector<double> xvec(ndim, 0.0);
        for(unsigned i=0; i<ndim; i++) {xvec[i] = intp.x()[i];}

        Vector<double> mvec(3, 0.0);
        for(unsigned i=0; i<3; i++) {mvec[i] = intp.m()[i];}

        // Source functions (for debugging, normally zero)
        const double phi_source = e_pt->get_phi_source(intp.time(),xvec);
        const double phi_1_source = e_pt->get_phi_1_source(intp.time(),xvec);

        // Fields at integration point
        Vector<double> h_applied, h_cryst_anis, h_magnetostatic;
        e_pt->get_applied_field(intp.time(), xvec, s, h_applied);
        e_pt->get_H_cryst_anis_field(intp.time(), xvec, mvec, h_cryst_anis);
        e_pt->get_magnetostatic_field(s, h_magnetostatic);

        Vector<double> h_simple(3, 0.0);
        h_simple[0] = h_applied[0] + h_cryst_anis[0] + h_magnetostatic[0];
        h_simple[1] = h_applied[1] + h_cryst_anis[1] + h_magnetostatic[1];
        h_simple[2] = h_applied[2] + h_cryst_anis[2] + h_magnetostatic[2];

        //======================================================================
        /// Use the above values to calculate the residuals
        //======================================================================

        // Loop over the test functions/nodes adding contributions to residuals
        for(unsigned l=0;l<n_node;l++)
          {

            // Cache test function derivative in a vector for easier access
            Vector<double> dtestdxl(ndim, 0.0);
            for(unsigned j=0; j<ndim; j++) dtestdxl[j] = intp.dtestdx(l, j);

            // Total potential (phi)
            const int phi_eqn = e_pt->nodal_local_eqn(l, e_pt->phi_index_micromag());
            // If value is not pinned and not on the boundary.( We need to treat
            // boundary values of phi differently since they are determined by
            // the boundary element matrix and phi_1.)
            if((phi_eqn >= 0) && (!(e_pt->node_pt(l)->is_on_boundary())))
              {
                std::cerr <<  "unpinned phis!" << std::endl;
                residuals[phi_eqn] -= phi_source*intp.test(l)*W; // source
                residuals[phi_eqn] -= itp_divm*intp.test(l)*W;         // div(m)
                for(unsigned k=0;k<ndim;k++)                       // Poisson
                  residuals[phi_eqn] -= intp.dphidx()[k]*intp.dtestdx(l,k)*W;
              }

            // Reduced potential (phi_1), only difference is in b.c.s
            const int phi_1_eqn = e_pt->nodal_local_eqn(l, e_pt->phi_1_index_micromag());
            if(phi_1_eqn >= 0)
              {
                std::cerr <<  "unpinned phis!" << std::endl;
                residuals[phi_1_eqn] -= phi_1_source*intp.test(l)*W;
                residuals[phi_1_eqn] -= itp_divm*intp.test(l)*W;
                for(unsigned k=0;k<ndim;k++)
                  residuals[phi_1_eqn] -= intp.dphi1dx()[k]*intp.dtestdx(l,k)*W;
              }

            // LL itself (m, time evolution)
            //=================================================

            Vector<double> gradmdotgradtest(3, 0.0);
            for(unsigned j=0; j<3; j++)
              {
                for(unsigned i=0; i<ndim; i++)
                  {
                    gradmdotgradtest[j] += intp.dmdx(j)[i] * intp.dtestdx(l, i);
                  }
              }


            // add to residual
            for(unsigned i=0; i<3; i++)
              {
                const int m_eqn = e_pt->nodal_local_eqn(l, e_pt->m_index_micromag(i));
                if(m_eqn >= 0)  // If it's not a boundary condition
                  {
                    // dmdt
                    residuals[m_eqn] += intp.dmdt()[i] * intp.test(l) * W;

                    // mxh for non-exchange fields
                    residuals[m_eqn] += opt_cross(i, intp.m(), h_simple)
                      * intp.test(l) * W;

                    // mxmxh for non-exchange fields
                    residuals[m_eqn] += ll_damp_c *
                      opt_double_cross(i, intp.m(), intp.m(), h_simple)
                      * intp.test(l) * W;

                    // mxex term
                    residuals[m_eqn] -= opt_cross(i, intp.m(), gradmdotgradtest)
                      * W;

                    // term 1 of mxmxex
                    residuals[m_eqn] += gradmdotgradtest[i] * W;

                    // term 2 of mxmxmex
                    double total = 0;
                    for(unsigned j=0; j<3; j++)
                      {
                        total += intp.m()[j] *
                          (intp.m()[i]*dot(dtestdxl, intp.dmdx(j), ndim)
                           + intp.test(l)*dot(intp.dmdx(i), intp.dmdx(j), ndim))

                          + intp.test(l) * intp.m()[i]
                          * dot(intp.dmdx(j), intp.dmdx(j), ndim);
                      }
                    residuals[m_eqn] -= ll_damp_c * total * W;
                  }
                else
                  {
                    std::cerr << "PINNED magnetisation! shouldn't be any of these!"
                              <<std::endl;
                  }
              }

          } // End of loop over test functions, end of residual calculations

        // Don't do Jacobian calculations yet
      }
  }


  void LLGResidualCalculator::fill_in_generic_residual_contribution
  (const MicromagEquations* const e_pt,
   Vector<double> &residuals, DenseMatrix<double> &jacobian,
   const unsigned& flag) const
  {

    //======================================================================
    /// Get some useful numbers and set up memory.
    //======================================================================

    // Find out how many nodes there are
    const unsigned n_node = e_pt->nnode();
    const unsigned ndim = e_pt->nodal_dimension();
    const unsigned eldim = e_pt->dim();

    // Get coefficients
    const double llg_precess_c = e_pt->llg_precession_coeff();
    const double llg_damp_c = e_pt->llg_damping_coeff();
    const double exch_c = e_pt->exchange_coeff();
    const double magstatic_c = e_pt->magnetostatic_coeff();

    // Cache time stepper weights needed in Jacobian (??ds write some
    // documentation on this? needs proper maths really...)
    double d_value_evaltime_by_dvalue_np1 =
      e_pt->node_pt(0)->time_stepper_pt()->weight(0,0);
    double d_valuederivative_evaltime_by_dvalue_np1 =
      e_pt->node_pt(0)->time_stepper_pt()->weight(1,0);

    //======================================================================
    /// Begin loop over the knots (integration points)
    //======================================================================
    for(unsigned ipt=0, nipt = e_pt->integral_pt()->nweight(); ipt<nipt; ipt++)
      {
        //======================================================================
        /// Calculate/get/interpolate all values for the residual calculations
        //======================================================================
        Vector<double> s(eldim);
        for(unsigned j=0; j<eldim; j++) {s[j] = e_pt->integral_pt()->knot(ipt,j);}

        // Create interpolator //??ds maybe should move this out of loop,
        // add .new_point() function or something?
        MMArrayInterpolator<5> intp(e_pt, s);

        double W = e_pt->integral_pt()->weight(ipt) * intp.j();

        // Calculate other things:

        // Divergence of m
        const double itp_divm = intp.div_m();

        Vector<double> xvec(ndim, 0.0);
        for(unsigned i=0; i<ndim; i++) {xvec[i] = intp.x()[i];}

        Vector<double> mvec(3, 0.0);
        for(unsigned i=0; i<3; i++) {mvec[i] = intp.m()[i];}

        // Source functions (for debugging, normally zero)
        const double phi_source = e_pt->get_phi_source(intp.time(),xvec);
        const double phi_1_source = e_pt->get_phi_1_source(intp.time(),xvec);

        // Fields at integration point
        Vector<double> h_applied, h_cryst_anis, h_magnetostatic;
        e_pt->get_applied_field(intp.time(), xvec, s, h_applied);
        e_pt->get_H_cryst_anis_field(intp.time(), xvec, mvec, h_cryst_anis);
        e_pt->get_magnetostatic_field(s, h_magnetostatic);


        //======================================================================
        /// Use the above values to calculate the residuals
        //======================================================================

        // Loop over the test functions/nodes adding contributions to residuals
        for(unsigned l=0;l<n_node;l++)
          {
            // Total potential (phi)
            const int phi_eqn = e_pt->nodal_local_eqn(l, e_pt->phi_index_micromag());
            // If value is not pinned and not on the boundary.( We need to treat
            // boundary values of phi differently since they are determined by
            // the boundary element matrix and phi_1.)
            if((phi_eqn >= 0) && (!(e_pt->node_pt(l)->is_on_boundary())))
              {
                std::cerr <<  "unpinned phis!" << std::endl;
                residuals[phi_eqn] -= phi_source*intp.test(l)*W; // source
                residuals[phi_eqn] -= itp_divm*intp.test(l)*W;         // div(m)
                for(unsigned k=0;k<ndim;k++)                       // Poisson
                  residuals[phi_eqn] -= intp.dphidx()[k]*intp.dtestdx(l,k)*W;
              }

            // Reduced potential (phi_1), only difference is in b.c.s
            const int phi_1_eqn = e_pt->nodal_local_eqn(l, e_pt->phi_1_index_micromag());
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
                const int m_eqn = e_pt->nodal_local_eqn(l, e_pt->m_index_micromag(i));
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
              m_eqn[i] = e_pt->nodal_local_eqn(l,e_pt->m_index_micromag(i));

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

            e_pt->get_hca_derivative(intp.time(),xvec,mvec,intp.psi(l2),dhcadm);

            //=========================================================
            /// Actual Jacobian calculation
            //=========================================================

            // Total potential (phi)
            const int phi_eqn = e_pt->nodal_local_eqn(l,e_pt->phi_index_micromag());
            if(phi_eqn >= 0)
              {
                const int phi_unknown = e_pt->nodal_local_eqn(l2,e_pt->phi_index_micromag());

                // The diagonal of the phi phi block is non-zero but it is
                // determined by the BEM method. Because of the way FEM assembly
                // is handled we cannot put the values in here. So for now we
                // put in dummy values to reserve space in the sparse matrix.
                if(e_pt->node_pt(l)->is_on_boundary())
                  {
                    if(phi_eqn == phi_unknown)
                      jacobian(phi_eqn,phi_unknown) = e_pt->DummyBEMControlledEntry;
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
                      const int m_unknown = e_pt->nodal_local_eqn(l2,e_pt->m_index_micromag(j));
                      if(m_unknown >= 0)
                        jacobian(phi_eqn,m_unknown) += - intp.dpsidx(l2,j) * intp.test(l) * W;
                    }
                  }

                // nothing w.r.t. phi_1 (dependence is in BEM matrix)
              }


            // Reduced potential (phi_1), only difference between this and phi
            // is in b.c.s
            const int phi_1_eqn = e_pt->nodal_local_eqn(l,e_pt->phi_1_index_micromag());
            if(phi_1_eqn >= 0){

              // w.r.t. phi_1
              const int phi_1_unknown = e_pt->nodal_local_eqn(l2,e_pt->phi_1_index_micromag());
              if(phi_1_unknown >= 0)
                jacobian(phi_1_eqn,phi_1_unknown) += - gradtestldotgradpsil2 * W;

              // w.r.t. m
              for(unsigned j=0; j<ndim; j++){
                const int m_unknown = e_pt->nodal_local_eqn(l2,e_pt->m_index_micromag(j));
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
                m_unknown[j] = e_pt->nodal_local_eqn(l2,e_pt->m_index_micromag(j));
              }

            // w.r.t. phi
            const int phi_unknown = e_pt->nodal_local_eqn(l2,e_pt->phi_index_micromag());
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
}
