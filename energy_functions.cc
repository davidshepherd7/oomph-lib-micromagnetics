
#include "../../src/generic/oomph_utilities.h"

#include "energy_functions.h"
#include "micromagnetics_element.h"
#include "../../src/generic/interpolator.h"


namespace oomph
{

  // double ElementalFunction::call(const GeneralisedElement* ele_pt,
  //                                const Vector<double> &s) const
  // {
  //   MMInterpolator intp(checked_dynamic_cast<const MicromagEquations*>(ele_pt),
  //                       s);
  //   return call(ele_pt, &intp);
  // }


  /// \short Calculate the energy due to exchange at a point at a single
  /// integration point.
  double ExchangeEnergyFunction::call(const GeneralisedElement* ele_pt,
                                      MMInterpolator* intp_pt) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Get parameters
    double A = m_ele_pt->magnetic_parameters_pt()->normalised_hex()/2;

    return A * (VectorOps::dot(intp_pt->dmdx(0), intp_pt->dmdx(0)) +
                VectorOps::dot(intp_pt->dmdx(1), intp_pt->dmdx(1)) +
                VectorOps::dot(intp_pt->dmdx(2), intp_pt->dmdx(2)));
  }


  /// \short Calculate the time derivative of energy due to exchange at a
  /// point at a single integration point.
  double dExchangeEnergydtFunction::call(const GeneralisedElement* ele_pt,
                                      MMInterpolator* intp_pt) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Get parameters
    double A = m_ele_pt->magnetic_parameters_pt()->normalised_hex()/2;

    return 2 *A * (VectorOps::dot(intp_pt->dmdx(0), intp_pt->d2mdxdt(0)) +
                   VectorOps::dot(intp_pt->dmdx(1), intp_pt->d2mdxdt(1)) +
                   VectorOps::dot(intp_pt->dmdx(2), intp_pt->d2mdxdt(2)));
  }


  /// \short Calculate energy due to external applied field at a single
  /// integration point.
  double ZeemanEnergyFunction::call(const GeneralisedElement* ele_pt,
                                    MMInterpolator* intp_pt) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Get the field
    Vector<double> h_applied;
    m_ele_pt->get_applied_field(intp_pt->time(), intp_pt->x(), intp_pt->s(), h_applied);

    return - VectorOps::dot(intp_pt->m(), h_applied);
  }

  /// \short Calculate time derivative of energy due to external applied
  /// field at a single integration point.
  double dZeemanEnergydtFunction::call(const GeneralisedElement* ele_pt,
                                       MMInterpolator* intp_pt) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Get the field
    Vector<double> h_applied;
    m_ele_pt->get_applied_field(intp_pt->time(), intp_pt->x(), intp_pt->s(), h_applied);

    return - VectorOps::dot(intp_pt->dmdt(), h_applied);
  }


  /// \short Calculate energy due to magnetocrystalline
  /// anisotropy at a single integration point.
  double CrystallineAnisotropyEnergyFunction::
  call(const GeneralisedElement* ele_pt, MMInterpolator* intp_pt) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Get parameters
    double k1 = m_ele_pt->magnetic_parameters_pt()->normalised_hk()/2;
    Vector<double> e = m_ele_pt->magnetic_parameters_pt()->easy_axis();

    return k1 * (1 - std::pow(VectorOps::dot(intp_pt->m(), e), 2));
  }


  /// \short Calculate time derivateive of energy due to magnetocrystalline
  /// anisotropy at a single integration point.
  double dCrystallineAnisotropydtEnergyFunction::
  call(const GeneralisedElement* ele_pt, MMInterpolator* intp_pt) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Get parameters
    double k1 = m_ele_pt->magnetic_parameters_pt()->normalised_hk()/2;
    Vector<double> e = m_ele_pt->magnetic_parameters_pt()->nearest_easy_axis(intp_pt->m());

    return - 2 * k1 * VectorOps::dot(intp_pt->m(), e) * VectorOps::dot(intp_pt->dmdt(), e);
  }

  /// \short Calculate energy due to external applied field at a single
  /// integration point.
  double MagnetostaticEnergyFunction::call(const GeneralisedElement* ele_pt,
                                           MMInterpolator* intp_pt) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Get the field
    Vector<double> h_ms;
    m_ele_pt->get_magnetostatic_field(intp_pt, h_ms);

    return -0.5 * VectorOps::dot(intp_pt->m(), h_ms);
  }


  /// \short Calculate time derivative of energy due to external applied
  /// field at a single integration point. Only for semi implicit version
  /// (but implicit is easy to implement..)
  double dMagnetostaticEnergydtFunction::call(const GeneralisedElement* ele_pt,
                                              MMInterpolator* intp_pt) const
  {
    const MicromagEquations* m_ele_pt
      = dynamic_cast<const MicromagEquations*>(ele_pt);

    // Get the field
    Vector<double> h_ms;
    m_ele_pt->get_magnetostatic_field(intp_pt, h_ms);

    // Get the time derivative of the field
    Vector<double> dh_ms_dt;
    m_ele_pt->get_magnetostatic_field_time_derivative(intp_pt, dh_ms_dt);

    return -0.5 * ( VectorOps::dot(intp_pt->dmdt(), h_ms) +
                    VectorOps::dot(intp_pt->m(), dh_ms_dt) );
  }

  /// \short Calculate (dm/dt)^2 for the previous time step at a single
  /// integration point.
  double DmdtSquaredFunction::call(const GeneralisedElement* ele_pt,
                                   MMInterpolator* intp_pt) const
  {
    return VectorOps::dot(intp_pt->dmdt(), intp_pt->dmdt());
  }

  /// \short Function for checking if integration is working ok, should
  /// give the area.
  double UnitFunction::call(const GeneralisedElement* ele_pt,
                            MMInterpolator* intp_pt) const
  {
    return 1.0;
  }

}
