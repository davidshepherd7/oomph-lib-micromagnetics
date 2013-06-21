
#include "energy_functions.h"
#include "../../src/generic/oomph_utilities.h"
#include "./micromagnetics_element.h"

using namespace oomph;

namespace oomph
{


  /// \short Calculate the energy due to exchange at a point.
  double ExchangeEnergyFunction::call(const GeneralisedElement* ele_pt,
                                      const Vector<double> &s) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Create interpolator
    MMInterpolator intp(m_ele_pt, s);

    // Get parameters
    double A = m_ele_pt->magnetic_parameters_pt()->exchange_constant();

    // Calculate
    double energy = A * (VectorOps::two_norm(intp.dmdx(0)) +
                         VectorOps::two_norm(intp.dmdx(1)) +
                         VectorOps::two_norm(intp.dmdx(2)));
    return energy;
  }


  /// \short Calculate energy due to external applied field. ??ds
  /// Implementation details mean that this is WRONG for semi implicit
  /// methods! :(
  double ZeemanEnergyFunction::
  call(const GeneralisedElement* ele_pt, const Vector<double> &s) const
  {
#ifdef PARANOID
    if(dynamic_cast<const SemiImplicitMicromagEquations*>(ele_pt) != 0)
      {
        std::string error_msg = "Doesn't work with semi implicit!";
        throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Create interpolator
    MMInterpolator intp(m_ele_pt, s);

    // Get the field
    Vector<double> h_applied;
    m_ele_pt->get_applied_field(intp.time(), intp.x(), s, h_applied);

    // Get "re-normalisation" parameters
    double M_s = m_ele_pt->magnetic_parameters_pt()->saturation_magnetisation();
    double H_magnitude =   // ??ds slightly nasty..
      1/(m_ele_pt->magnetic_parameters_pt()->field_normalisation_factor());

    return - mag_parameters::mu0 * M_s * H_magnitude
      * VectorOps::dot(intp.m(), h_applied);
  }


  /// \short Calculate energy at a point s due to magnetocrystalline
  /// anisotropy. Source: my writeup (david).
  double CrystallineAnisotropyEnergyFunction::
  call(const GeneralisedElement* ele_pt, const Vector<double> &s) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Create interpolator
    MMInterpolator intp(m_ele_pt, s);

    // Get dimensional parameters
    double M_s = m_ele_pt->magnetic_parameters_pt()->saturation_magnetisation();
    double k1 = m_ele_pt->magnetic_parameters_pt()->k1();
    double mu0 = mag_parameters::mu0;
    Vector<double> e = m_ele_pt->magnetic_parameters_pt()->easy_axis(intp.m());

    // Return energy
    return k1 * (1 - (VectorOps::dot(intp.m(), e)/ ( mu0*mu0*M_s )));
  }


  /// \short Calculate energy due to external applied field. ??ds
  /// Implementation details mean that this is WRONG for semi implicit
  /// methods! :(
  double MagnetostaticEnergyFunction::call(const GeneralisedElement* ele_pt,
                                           const Vector<double> &s) const
  {
#ifdef PARANOID
    if(dynamic_cast<const SemiImplicitMicromagEquations*>(ele_pt) != 0)
      {
        std::string error_msg = "Doesn't work with semi implicit!";
        throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif


    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Create interpolator
    MMInterpolator intp(m_ele_pt, s);

    // Get the field
    Vector<double> h_ms(3, 0.0);
    for(unsigned j=0; j<m_ele_pt->nodal_dimension(); j++)
      {
        h_ms[j] = -1 * intp.dphidx()[j] *
          m_ele_pt->magnetic_parameters_pt()->magnetostatic_debug_coeff();
      }

    // Get "re-normalisation" parameters
    double M_s = m_ele_pt->magnetic_parameters_pt()->saturation_magnetisation();
    double H_magnitude =   // ??ds slightly nasty..
      1/(m_ele_pt->magnetic_parameters_pt()->field_normalisation_factor());

    return -0.5 * mag_parameters::mu0 * M_s * H_magnitude
      * VectorOps::dot(intp.m(), h_ms);
  }



}
