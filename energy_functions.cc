
#include "energy_functions.h"
#include "../../src/generic/oomph_utilities.h"
#include "./micromagnetics_element.h"

namespace oomph
{


  /// \short Calculate the energy due to exchange at a point at a single
  /// integration point.
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


  /// \short Calculate energy due to external applied field at a single
  /// integration point.
  double ZeemanEnergyFunction::
  call(const GeneralisedElement* ele_pt, const Vector<double> &s) const
  {
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


  /// \short Calculate energy due to magnetocrystalline
  /// anisotropy at a single integration point.
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
    Vector<double> e = m_ele_pt->magnetic_parameters_pt()->easy_axis(intp.m());

    // Return energy
    double temp = (VectorOps::dot(intp.m(), e) / M_s);
    return k1 * temp * temp;
  }


  /// \short Calculate energy due to external applied field at a single
  /// integration point.
  double MagnetostaticEnergyFunction::call(const GeneralisedElement* ele_pt,
                                           const Vector<double> &s) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Create interpolator
    MMInterpolator intp(m_ele_pt, s);

    // Get the field
    Vector<double> h_ms;
    m_ele_pt->get_magnetostatic_field(s, h_ms);

    // Get dimensional parameters
    double M_s = m_ele_pt->magnetic_parameters_pt()->saturation_magnetisation();
    double H_magnitude =   // ??ds slightly nasty..
      1/(m_ele_pt->magnetic_parameters_pt()->field_normalisation_factor());

    return -0.5 * mag_parameters::mu0 * M_s * H_magnitude
      * VectorOps::dot(intp.m(), h_ms);
  }


  /// \short Calculate (dm/dt)^2 for the previous time step at a single
  /// integration point.
  double DmdtSquaredFunction::call(const GeneralisedElement* ele_pt,
                                   const Vector<double> &s) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);
    MMInterpolator intp(m_ele_pt, s);
    return VectorOps::dot(intp.dmdt(), intp.dmdt());
  }

}
