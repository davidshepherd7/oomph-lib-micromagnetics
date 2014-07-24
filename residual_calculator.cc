#include "./residual_calculator.h"
#include "./micromagnetics_element.h"
#include "./vector_helpers.h"

#include "llg_factories.h"
#include "new_interpolators.h"

using namespace oomph;
using namespace VectorOps;

namespace oomph
{

  void LLGResidualCalculator::ll_residual
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

    // Factor to rescale time s.t. it matches with Gilbert form of llg
    const double ll_conversion_factor = (1+llg_damp_c*llg_damp_c);

    // Create interpolator
    std::unique_ptr<CachingMMArrayInterpolator>
      intp_pt(Factories::array_interpolator_factory(e_pt));


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

        // Set up interpolator for this point
        intp_pt->build(s);

        double W = e_pt->integral_pt()->weight(ipt) * intp_pt->j();

        // Calculate other things:

        // cache pointers from interpolator
        const double* intp_m = intp_pt->m();
        const double* intp_dmdt = intp_pt->dmdt();
        const double time = intp_pt->time();
        const double* intp_x = intp_pt->x();
        const double* intp_dmdx[3] = {intp_pt->dmdx(0),
                                      intp_pt->dmdx(1),
                                      intp_pt->dmdx(2)};

        // Copy some things into vectors ready for use in function calls
        // ??ds get rid of this?
        Vector<double> xvec(ndim, 0.0);
        for(unsigned i=0; i<ndim; i++) {xvec[i] = intp_x[i];}
        Vector<double> mvec(3, 0.0);
        for(unsigned i=0; i<3; i++) {mvec[i] = intp_m[i];}

        // Source functions (for debugging, normally zero)
        const double phi_source = e_pt->get_phi_source(time,xvec);
        const double phi_1_source = e_pt->get_phi_1_source(time,xvec);

        // Fields at integration point
        Vector<double> h_cryst_anis, h_magnetostatic;
        Vector<double> h_applied = e_pt->get_applied_field(time, xvec);
        e_pt->get_H_cryst_anis_field(time, xvec, mvec, h_cryst_anis);
        e_pt->get_magnetostatic_field(s, h_magnetostatic);
        # warning using slow get magnetostaic field call
        //??ds convert to use new interpolators

        Vector<double> h_simple(3, 0.0);
        h_simple[0] = h_applied[0] + h_cryst_anis[0] + h_magnetostatic[0];
        h_simple[1] = h_applied[1] + h_cryst_anis[1] + h_magnetostatic[1];
        h_simple[2] = h_applied[2] + h_cryst_anis[2] + h_magnetostatic[2];

        //======================================================================
        /// Use the above values to calculate the residuals
        //======================================================================

        // Loop over the test functions/nodes adding contributions to residuals
        for(unsigned l=0; l<n_node; l++)
          {

            // Cache test function derivative in a vector for easier access
            // ??ds remove?
            Vector<double> dtestdxl(ndim, 0.0);
            for(unsigned j=0; j<ndim; j++) dtestdxl[j] = intp_pt->dtestdx(l, j);

            // Total potential (phi)
            const int phi_eqn = e_pt->nodal_local_eqn(l, e_pt->phi_index_micromag());
            // If value is not pinned and not on the boundary.( We need to treat
            // boundary values of phi differently since they are determined by
            // the boundary element matrix and phi_1.)
            if((phi_eqn >= 0) && (!(e_pt->node_pt(l)->is_on_boundary())))
              {
                residuals[phi_eqn] -= phi_source*intp_pt->test(l)*W; // source
                residuals[phi_eqn] -= intp_pt->div_m()*intp_pt->test(l)*W;         // div(m)
                for(unsigned k=0;k<ndim;k++)                       // Poisson
                  residuals[phi_eqn] -= intp_pt->dphidx()[k]*intp_pt->dtestdx(l,k)*W;
              }

            // Reduced potential (phi_1), only difference is in b.c.s
            const int phi_1_eqn = e_pt->nodal_local_eqn(l, e_pt->phi_1_index_micromag());
            if(phi_1_eqn >= 0)
              {
                residuals[phi_1_eqn] -= phi_1_source*intp_pt->test(l)*W;
                residuals[phi_1_eqn] -= intp_pt->div_m()*intp_pt->test(l)*W;
                for(unsigned k=0;k<ndim;k++)
                  residuals[phi_1_eqn] -= intp_pt->dphi1dx()[k]*intp_pt->dtestdx(l,k)*W;
              }

            // LL itself (m, time evolution)
            //=================================================

            Vector<double> gradmdotgradtest(3, 0.0);
            for(unsigned j=0; j<3; j++)
              {
                gradmdotgradtest[j] = dot(intp_dmdx[j], dtestdxl, ndim);
              }


            // add to residual
            for(unsigned i=0; i<3; i++)
              {
                const int m_eqn = e_pt->nodal_local_eqn(l, e_pt->m_index_micromag(i));
                if(m_eqn >= 0)  // If it's not a boundary condition
                  {
                    // dmdt
                    residuals[m_eqn] += ll_conversion_factor *intp_dmdt[i]
                      * intp_pt->test(l) * W;

                    // mxh for non-exchange fields (precession)
                    residuals[m_eqn] += opt_cross(i, intp_m, h_simple)
                      * intp_pt->test(l) * W;

                    // mxmxh for non-exchange fields (damping)
                    residuals[m_eqn] += llg_damp_c *
                      opt_double_cross(i, intp_m, intp_m, h_simple)
                      * intp_pt->test(l) * W;

                    // mxex term (precession)
                    residuals[m_eqn] -= opt_cross(i, intp_m, gradmdotgradtest)
                      * W;

                    // term 1 of mxmxex (damping)
                    residuals[m_eqn] += llg_damp_c * gradmdotgradtest[i] * W;

                    // term 2 of mxmxex (damping)
                    double sum = 0;
                    for(unsigned j=0; j<3; j++)
                      {
                        sum += intp_m[j] *
                          (intp_m[i]*dot(dtestdxl, intp_dmdx[j], ndim)
                           + intp_pt->test(l)*dot(intp_dmdx[i], intp_dmdx[j], ndim))

                          + intp_pt->test(l) * intp_m[i]
                          * dot(intp_dmdx[j], intp_dmdx[j], ndim);
                      }
                    residuals[m_eqn] -= llg_damp_c * sum * W;
                  }
              }

          } // End of loop over test functions, end of residual calculations


        //======================================================================
        // If we want the Jacobian as well then calculate it (otherwise
        // continue to next integration point). Only for decoupled phi
        // terms so far! Doesn't handle BEM!
        // =====================================================================
        if(!flag) continue;


        // Double loop over nodes for the jacobian
        for(unsigned l=0; l<n_node; l++){

          for(unsigned l2=0;l2<n_node;l2++){

            double gradtestldotgradpsil2 = 0.0;
            for(unsigned i=0; i < ndim; i++)
              {
                gradtestldotgradpsil2 += intp_pt->dtestdx(l,i) * intp_pt->dpsidx(l2,i);
              }

            // Total potential (phi)
            const int phi_eqn = e_pt->nodal_local_eqn(l,e_pt->phi_index_micromag());
            if(phi_eqn >= 0)
              {
                const int phi_unknown = e_pt->nodal_local_eqn(l2,e_pt->phi_index_micromag());

                // w.r.t. phi
                if(phi_unknown >= 0)
                  {
                    jacobian(phi_eqn,phi_unknown) += -gradtestldotgradpsil2 * W;
                  }
#ifdef PARANOID
                const int phi_1_unknown =
                  e_pt->nodal_local_eqn(l2,e_pt->phi_1_index_micromag());
                if(phi_1_unknown >= 0)
                  {
                    throw OomphLibError("dphi/dphi1 not yet implemented",
                                        OOMPH_EXCEPTION_LOCATION,
                                        OOMPH_CURRENT_FUNCTION);
                  }
#endif
              }


            // Reduced potential (phi_1), only difference between this and phi
            // is in b.c.s
            const int phi_1_eqn = e_pt->nodal_local_eqn(l,e_pt->phi_1_index_micromag());
            if(phi_1_eqn >= 0)
              {
                // w.r.t. phi_1
                const int phi_1_unknown = e_pt->nodal_local_eqn(l2,e_pt->phi_1_index_micromag());
                if(phi_1_unknown >= 0)
                  {
                    jacobian(phi_1_eqn,phi_1_unknown) += - gradtestldotgradpsil2 * W;
                  }

#ifdef PARANOID
                // w.r.t. phi
                const int phi_unknown =
                  e_pt->nodal_local_eqn(l2,e_pt->phi_index_micromag());
                if(phi_unknown >= 0)
                  {
                    std::cout << phi_unknown << " " << phi_1_eqn << std::endl;
                    throw OomphLibError("dphi1/dphi not yet implemented",
                                        OOMPH_EXCEPTION_LOCATION,
                                        OOMPH_CURRENT_FUNCTION);
                  }
#endif
              }



#ifdef PARANOID
            int m_eqn[3], m_unknown[3];

            for(unsigned j=0; j<3; j++)
              {
                m_unknown[j] = e_pt->nodal_local_eqn(l2,e_pt->m_index_micromag(j));
                m_eqn[j] = e_pt->nodal_local_eqn(l,e_pt->m_index_micromag(j));
              }
            // w.r.t. m
            for(unsigned j=0; j<3; j++) // loop over the m we differentiate by
              {
                if((m_unknown[j] >= 0) || (m_eqn[j] >= 0))
                  {
                    throw OomphLibError("Jacobian calculation not yet implemented for magnetism part of ll (only for magnetostatics).",
                                        OOMPH_EXCEPTION_LOCATION,
                                        OOMPH_CURRENT_FUNCTION);
                  }

              }
#endif

          }
        }// End of Jacobian calculations

      }

  }


  void LLGResidualCalculator::llg_residual
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

    // Create interpolator
    std::unique_ptr<CachingMMArrayInterpolator>
      intp_pt(Factories::array_interpolator_factory(e_pt));

    // Allocate vectors outside intergration loop
    Vector<double> xvec(ndim, 0.0),  mvec(3, 0.0), h_cryst_anis,
      h_magnetostatic, s(eldim);

    //======================================================================
    /// Begin loop over the knots (integration points)
    //======================================================================
    for(unsigned ipt=0, nipt = e_pt->integral_pt()->nweight(); ipt<nipt; ipt++)
      {
        //======================================================================
        /// Calculate/get/interpolate all values for the residual calculations
        //======================================================================
        for(unsigned j=0; j<eldim; j++) {s[j] = e_pt->integral_pt()->knot(ipt,j);}

        // Set up interpolator for this point
        intp_pt->build(s);

        double W = e_pt->integral_pt()->weight(ipt) * intp_pt->j();

        // cache pointers from interpolator for speed of access
        const double* intp_m = intp_pt->m();
        const double* intp_dmdt = intp_pt->dmdt();
        const double time = intp_pt->time();
        const double* intp_x = intp_pt->x();
        const double* intp_dmdx[3] = {intp_pt->dmdx(0),
                                      intp_pt->dmdx(1),
                                      intp_pt->dmdx(2)};

        // Copy some things into vectors ready for use in function calls
        for(unsigned i=0; i<ndim; i++) {xvec[i] = intp_x[i];}
        for(unsigned i=0; i<3; i++) {mvec[i] = intp_m[i];}

        // Source functions (for debugging, normally zero)
        const double phi_source = e_pt->get_phi_source(time,xvec);
        const double phi_1_source = e_pt->get_phi_1_source(time,xvec);

        // Fields at integration point
        Vector<double> h_applied = e_pt->get_applied_field(time, xvec);
        e_pt->get_H_cryst_anis_field(time, xvec, mvec, h_cryst_anis);
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
                residuals[phi_eqn] -= phi_source*intp_pt->test(l)*W; // source
                residuals[phi_eqn] -= intp_pt->div_m()*intp_pt->test(l)*W;         // div(m)
                for(unsigned k=0;k<ndim;k++)                       // Poisson
                  residuals[phi_eqn] -= intp_pt->dphidx()[k]*intp_pt->dtestdx(l,k)*W;
              }

            // Reduced potential (phi_1), only difference is in b.c.s
            const int phi_1_eqn = e_pt->nodal_local_eqn(l, e_pt->phi_1_index_micromag());
            if(phi_1_eqn >= 0)
              {
                residuals[phi_1_eqn] -= phi_1_source*intp_pt->test(l)*W;
                residuals[phi_1_eqn] -= intp_pt->div_m()*intp_pt->test(l)*W;
                for(unsigned k=0;k<ndim;k++)
                  residuals[phi_1_eqn] -= intp_pt->dphi1dx()[k]*intp_pt->dtestdx(l,k)*W;
              }

            // LLG itself (m, time evolution)
            //=================================================

            // Exchange residual contribution after integration by parts:
            //??ds possibly could optimise - this is calculated twice
            double gradtestdotgradmi[3] = {0,0,0};
            for(unsigned i=0; i<3; i++)
              for(unsigned j=0; j<ndim; j++)
                gradtestdotgradmi[i] += intp_pt->dtestdx(l,j) * intp_dmdx[i][j];

            // add to residual
            for(unsigned i=0; i<3; i++)
              {
                const int m_eqn = e_pt->nodal_local_eqn(l, e_pt->m_index_micromag(i));
                if(m_eqn >= 0)  // If it's not a boundary condition
                  {
                    // dmdt, mxh_ap, mxh_ca, mxh_ms and mxdmdt terms
                    residuals[m_eqn] +=
                      ( intp_dmdt[i]
                        + llg_precess_c * opt_cross(i, intp_m, h_applied)
                        + llg_precess_c * opt_cross(i, intp_m, h_cryst_anis)
                        + llg_precess_c * opt_cross(i, intp_m, h_magnetostatic)
                        - llg_damp_c * opt_cross(i, intp_m, intp_dmdt)
                        )*intp_pt->test(l)*W;

                    // (m x exchange) term (separate because it involves
                    // derivatives of the test function).
                    residuals[m_eqn] -= exch_c * llg_precess_c *
                      opt_cross(i, intp_m, gradtestdotgradmi) * W;
                  }
              }

          } // End of loop over test functions, end of residual calculations


        //======================================================================
        /// If we want the Jacobian as well then calculate it (otherwise
        /// continue to next integration point).
        //======================================================================
        if(!flag) continue;
        // Note that "continue" means go on to the next step of the loop,
        // which is perhaps not entirely intuitive. We use "if(not thing)
        // then continue" statements rather than "if(thing) then { ... }"
        // because the amount of braces and indentation gets ridiculous
        // with so many nested for loops and if statements.


        // nondiffterms precalcs
        double nondiffterms[3];
        for(unsigned i=0; i<3; i++)
          nondiffterms[i] = - llg_damp_c * intp_dmdt[i]
            + llg_precess_c * h_magnetostatic[i]
            + llg_precess_c * (h_cryst_anis[i] + h_applied[i]);

        double gradtestdotgradmi[3];
        double diffterms[3];
        int m_eqn[3], m_unknown[3];
        double dhcadm[3][3];

        // Need a copy of grad psi here because the 3rd entry is used even
        // in 1d/2d, so we need to fill in a zero. Zero all entries then
        // copy whichever real entries we have.
        double gradpsil2[3] = {0.0, 0.0, 0.0};

        // Double loop over nodes for the jacobian
        for(unsigned l=0; l<n_node; l++){

          // Pre calculate values with no l2 dependencies
          for(unsigned i=0; i<3; i++)
            {
              m_eqn[i] = e_pt->nodal_local_eqn(l,e_pt->m_index_micromag(i));

              gradtestdotgradmi[i] = 0.0;
              for(unsigned j=0; j<ndim; j++) //??ds repeated calculation..
                {
                  gradtestdotgradmi[i] += intp_pt->dtestdx(l,j) * intp_dmdx[i][j];
                }
            }

          for(unsigned l2=0;l2<n_node;l2++){

            for(unsigned j=0; j<ndim; j++) {gradpsil2[j] = intp_pt->dpsidx(l2,j);}

            double gradtestldotgradpsil2 = 0.0;
            for(unsigned i=0; i < ndim; i++)
              {
                gradtestldotgradpsil2 += intp_pt->dtestdx(l,i) * intp_pt->dpsidx(l2,i);
              }

            e_pt->get_hca_derivative(time,xvec,mvec,intp_pt->psi(l2),dhcadm);

            //=========================================================
            /// Actual Jacobian calculation
            //=========================================================

            // Total potential (phi)
            const int phi_eqn = e_pt->nodal_local_eqn(l,e_pt->phi_index_micromag());
            if(phi_eqn >= 0)
              {
                const int phi_unknown = e_pt->nodal_local_eqn(l2,e_pt->phi_index_micromag());

                // The diagonal of the phi phi block is non-zero but it is
                // determined by the BEM method. The block is minus
                // identity matrix if the residuals are just r_nodal =
                // G[phi_1_b] - phi_b. Nodal values used because bem is
                // done by co-location rather than Galerkin to keep the bem
                // matrix simpler. Currently add this in the problem to avoid overlap issues.
                if(e_pt->node_pt(l)->is_on_boundary())
                  {
                    if(phi_eqn == phi_unknown)
                      jacobian(phi_eqn,phi_unknown) = 0.0;
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
                        jacobian(phi_eqn,m_unknown) += - intp_pt->dpsidx(l2,j) * intp_pt->test(l) * W;
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
                  jacobian(phi_1_eqn,m_unknown) += - intp_pt->dpsidx(l2,j) * intp_pt->test(l) * W;
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
                    * magstatic_c * opt_cross(j, intp_m, gradpsil2) * W * intp_pt->test(l);
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
                diffterms[j] -= llg_damp_c * intp_pt->psi(l2) *
                  d_valuederivative_evaltime_by_dvalue_np1;

                // mass matrix component due to time derivative
                jacobian(m_eqn[j],m_unknown[j])
                  += intp_pt->test(l) * intp_pt->psi(l2) *
                  d_valuederivative_evaltime_by_dvalue_np1 * W;

                for(unsigned i=0; i<3; i++) // loop over the m we differentiate w.r.t.
                  {
                    if(!(m_eqn[i] >= 0)) continue;

                    // chain rule gives:
                    // d/dmj( m x (...)) = dmidmj x (....) + m x d/dmj(.....)

                    // dmidmj x (....)
                    jacobian(m_eqn[i],m_unknown[j]) +=
                      W * intp_pt->test(l) * intp_pt->psi(l2) * opt_cross(i, jhat, nondiffterms)
                      * d_value_evaltime_by_dvalue_np1;

                    // m x d/dmj(.....)
                    jacobian(m_eqn[i], m_unknown[j]) +=
                      W * intp_pt->test(l) * opt_cross(i, intp_m, diffterms);

                    // Exchange contribution
                    jacobian(m_eqn[i],m_unknown[j])
                      -= llg_precess_c * exch_c * W *
                      ( intp_pt->psi(l2) * opt_cross(i, jhat, gradtestdotgradmi)

                        + opt_cross(i, intp_m, jhat) * gradtestldotgradpsil2
                        )
                      * d_value_evaltime_by_dvalue_np1;

                  }

              }

          }
        }// End of Jacobian calculations

      }// End of loop over integration points

  } // End of fill in residuals function
}
