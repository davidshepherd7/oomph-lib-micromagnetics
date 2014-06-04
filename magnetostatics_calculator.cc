
#include "magnetostatics_calculator.h"

#include "../../src/generic/Vector.h"
#include "../../src/generic/oomph_utilities.h"
#include "micromagnetics_element.h"

namespace oomph
{
  void ImplicitMagnetostaticsCalculator::
  get_magnetostatic_field(MMArrayInterpolator* intp_pt,
                          Vector<double> &hms) const
  {
    // Make sure the field has 3 dimensions (even if there are only two
    // spatial dimensions).
    hms.resize(3, 0.0);

    const double* hms_temp;
    // Copy the derivative elements into the field vector (it only has
    // [nodal dimension] entries).
    hms_temp = intp_pt->dphidx();

    // Multiply by -1
    for(unsigned j=0; j<3; j++)
      {
        hms[j] = -1 * hms_temp[j];
      }
  }

  void ImplicitMagnetostaticsCalculator::get_magnetostatic_field_time_derivative
  (MMInterpolator* intp_pt, Vector<double> &dh_ms_dt) const
  {
    // Copy the derivative elements into the field vector (it only has
    // [nodal dimension] entries).
    dh_ms_dt = intp_pt->d2valuedxdt(phi_index_micromag());

    // Make sure the field has 3 dimensions (even if there are only two
    // spatial dimensions).
    dh_ms_dt.resize(3, 0.0);

    // Multiply by -1
    for(unsigned j=0; j<3; j++)
      {
        dh_ms_dt[j] *= -1;
      }
  }

  void SemiImplicitMagnetostaticsCalculator::
  get_magnetostatic_field(MMArrayInterpolator* intp_pt,
                          Vector<double> &h_ms) const
  {
    // Get magnetostatic field from field element. Safe to assume that all
    // nodes have the same time stepper becuase otherwise our interpolators
    // don't work.
    magnetostatic_field_element_pt()->magnetostatic_field
      (intp_pt->s(), intp_pt->ts_pt(), h_ms);
  }


  void SemiImplicitMagnetostaticsCalculator::get_magnetostatic_field_time_derivative
  (MMInterpolator* intp_pt, Vector<double> &dh_ms_dt) const
  {
    // Get magnetostatic field derivative from field element
    magnetostatic_field_element_pt()->
      magnetostatic_field_time_derivative(intp_pt->s(), intp_pt->ts_pt(),
                                          dh_ms_dt);
  }

  void AnalyticalMagnetostatics::
  get_magnetostatic_field(MMArrayInterpolator* intp_pt,
                          Vector<double> &h_ms) const
  {
#ifdef PARANOID
    if(Magnetostatic_field_fct_pt == 0)
      {
        std::string err = "Magnetostatic_field_fct_pt is null!";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif

    // Copy to vectors
    Vector<double> m(3, 0.0), x(intp_pt->dim(), 0.0);
    for(unsigned j=0; j<intp_pt->dim(); j++)
      {
        m[j] = intp_pt->m()[j];
        x[j] = intp_pt->x()[j];
      }

    h_ms = Magnetostatic_field_fct_pt(intp_pt->time(), x, m);

  }

}
