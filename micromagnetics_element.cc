
#include "micromagnetics_element.h"
#include "micromagnetics_flux_element.h"

// For output in
#include <iomanip>

using namespace oomph;
using namespace MathematicalConstants;
using namespace VectorOps;

namespace oomph
{
  /// \short Integrate a function given by func_pt over the element using
  /// the given integral_pt(). Because C++ sucks we have to do this with
  /// weird function objects.
  double MicromagEquations::integrate_over_element(const ElementalFunction* func_pt) const
  {
     // Always do time derivative interpolation using BDF2 timestepper
     // to increase accuracy (midpoint is lower order for time
     // derivatives and implementing optional timestepper function
     // argument everywhere would be awful because C++ sucks).

     //??ds increase order of BDF? would require implementing it myself...

     // Create bdf time stepper, check some things and set it up
     BDF<2> bdf;
     TimeStepper* node_ts_pt = node_pt(0)->time_stepper_pt();
 #ifdef PARANOID
     // Check we have enough stored time steps
     if(bdf.ndt() > node_ts_pt->time_pt()->ndt())
       {
         std::string error_msg = "Not enough time steps in node's time stepper to use BDF2 for time derivative calculations.";
         throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                             OOMPH_EXCEPTION_LOCATION);
       }

     // Try to check if time steppers are compatible... ??ds
 #warning No check for compatability of time stepper with bdf2 used for integrals!
     if(0)
       {
         std::string error_msg = "";
         throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                             OOMPH_EXCEPTION_LOCATION);
       }
 #endif

     // Finish off time stepper construction
     bdf.time_pt() = node_ts_pt->time_pt();
     bdf.set_weights();


    double result = 0;

    // Loop over knots and sum
    for(unsigned ipt=0, nipt = integral_pt()->nweight(); ipt<nipt; ipt++)
      {
        // Get position in element
        Vector<double> s(this->dim());
        for(unsigned j=0; j<this->dim(); j++)
          {s[j] = integral_pt()->knot(ipt,j);}

        MMInterpolator intp(this, s, &bdf);
        double J = intp.j();
        double w = integral_pt()->weight(ipt);

        // Add contribution
        result += func_pt->call(this, &intp) * w * J;
      }

    return result;
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
   const unsigned& flag)
  {
    if(flag && Use_fd_jacobian)
      {
        // Generalised element version is by FD
        FiniteElement::fill_in_contribution_to_jacobian(residuals,
                                                        jacobian);
      }
    else
      {
        // Just run the function stored in the residual calculator slot
        Residual_calculator_pt->fill_in_generic_residual_contribution
          (this, residuals, jacobian, flag);
      }
  }


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

        // Output the magnetostatic at this point (for both decoupled and
        // fully coupled methods).
        Vector<double> h_magnetostatic;
        get_magnetostatic_field(s, h_magnetostatic);
        h_magnetostatic.resize(3, 0.0); // make sure it has a 3rd entry
        for(unsigned i=0; i<3; i++) outfile << h_magnetostatic[i] << " ";

        // Phi 1 at this point
        outfile << intp.phi1() << " ";

        // Phi total at this point
        outfile << intp.phi() << " ";

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
                          Vector<double> &hms) const
  {
#ifdef PARANOID
    if(intp_pt == 0)
      {
        std::string error_msg = "Null interpolator!";
        throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

    // Make sure the field has 3 dimensions (even if there are only two
    // spatial dimensions).
    hms.resize(3, 0.0);

    const double* hms_temp;
    // Copy the derivative elements into the field vector (it only has
    // [nodal dimension] entries).
    hms_temp = intp_pt->dvaluedx(this->phi_index_micromag());

    // Multiply by -1 and normalise
    for(unsigned j=0; j<3; j++)
      {
        hms[j] = -1 * magnetostatic_coeff() * hms_temp[j];
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

    // Get magnetostatic field from field element. Safe to assume that all
    // nodes have the same time stepper becuase otherwise our interpolators
    // don't work.
    magnetostatic_field_element_pt()->magnetostatic_field
      (intp_pt->s(), node_pt(0)->time_stepper_pt(), h_ms);
  }

  void MicromagEquations::
  get_magnetostatic_field_time_derivative(MMInterpolator* intp_pt,
                                          Vector<double> &dh_ms_dt) const
  {
#ifdef PARANOID
    if(intp_pt == 0)
      {
        std::string error_msg = "Null interpolator!";
        throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

    // Copy the derivative elements into the field vector (it only has
    // [nodal dimension] entries).
    dh_ms_dt = intp_pt->d2valuedxdt(this->phi_index_micromag());

    // Make sure the field has 3 dimensions (even if there are only two
    // spatial dimensions).
    dh_ms_dt.resize(3, 0.0);

    // Multiply by -1 and normalise
    for(unsigned j=0; j<3; j++)
      {
        dh_ms_dt[j] *= -1 * magnetostatic_coeff();
      }
  }

  void SemiImplicitMicromagEquations::
  get_magnetostatic_field_time_derivative(MMInterpolator* intp_pt,
                                          Vector<double> &dh_ms_dt) const
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

    // Get magnetostatic field derivative from field element
    magnetostatic_field_element_pt()->magnetostatic_field_time_derivative(intp_pt->s(),
                                                                          intp_pt->ts_pt(),
                                                                          dh_ms_dt);
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
