
#include "../../src/generic/oomph_utilities.h"

#include "energy_functions.h"
#include "micromagnetics_element.h"
#include "../../src/generic/interpolator.h"


namespace oomph
{

  double ElementalFunction::call(const GeneralisedElement* ele_pt,
                                 const Vector<double> &s) const
  {
    MMInterpolator intp(checked_dynamic_cast<const MicromagEquations*>(ele_pt),
                        s);
    return call(ele_pt, &intp);
  }


  /// \short Calculate the energy due to exchange at a point at a single
  /// integration point.
  double ExchangeEnergyFunction::call(const GeneralisedElement* ele_pt,
                                      MMInterpolator* intp) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Get parameters
    double A = m_ele_pt->magnetic_parameters_pt()->normalised_hex()/2;

    return A * (VectorOps::dot(intp->dmdx(0), intp->dmdx(0)) +
                VectorOps::dot(intp->dmdx(1), intp->dmdx(1)) +
                VectorOps::dot(intp->dmdx(2), intp->dmdx(2)));
  }


  /// \short Calculate the time derivative of energy due to exchange at a
  /// point at a single integration point.
  double dExchangeEnergydtFunction::call(const GeneralisedElement* ele_pt,
                                      MMInterpolator* intp) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Get parameters
    double A = m_ele_pt->magnetic_parameters_pt()->normalised_hex()/2;

    return 2 *A * (VectorOps::dot(intp->dmdx(0), intp->d2mdxdt(0)) +
                   VectorOps::dot(intp->dmdx(1), intp->d2mdxdt(1)) +
                   VectorOps::dot(intp->dmdx(2), intp->d2mdxdt(2)));
  }


  /// \short Calculate energy due to external applied field at a single
  /// integration point.
  double ZeemanEnergyFunction::call(const GeneralisedElement* ele_pt,
                                    MMInterpolator* intp) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Get the field
    Vector<double> h_applied;
    m_ele_pt->get_applied_field(intp->time(), intp->x(), intp->s(), h_applied);

    return - VectorOps::dot(intp->m(), h_applied);
  }

  /// \short Calculate time derivative of energy due to external applied
  /// field at a single integration point.
  double dZeemanEnergydtFunction::call(const GeneralisedElement* ele_pt,
                                       MMInterpolator* intp) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Get the field
    Vector<double> h_applied;
    m_ele_pt->get_applied_field(intp->time(), intp->x(), intp->s(), h_applied);

    return - VectorOps::dot(intp->dmdt(), h_applied);
  }


  /// \short Calculate energy due to magnetocrystalline
  /// anisotropy at a single integration point.
  double CrystallineAnisotropyEnergyFunction::
  call(const GeneralisedElement* ele_pt, MMInterpolator* intp) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Get parameters
    double k1 = m_ele_pt->magnetic_parameters_pt()->normalised_hk()/2;
    Vector<double> e = m_ele_pt->magnetic_parameters_pt()->easy_axis();

    return k1 * (1 - std::pow(VectorOps::dot(intp->m(), e), 2));
  }


  /// \short Calculate time derivateive of energy due to magnetocrystalline
  /// anisotropy at a single integration point.
  double dCrystallineAnisotropydtEnergyFunction::
  call(const GeneralisedElement* ele_pt, MMInterpolator* intp) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Get parameters
    double k1 = m_ele_pt->magnetic_parameters_pt()->normalised_hk()/2;
    Vector<double> e = m_ele_pt->magnetic_parameters_pt()->nearest_easy_axis(intp->m());

    return - 2 * k1 * VectorOps::dot(intp->m(), e) * VectorOps::dot(intp->dmdt(), e);
  }

  /// \short Calculate energy due to external applied field at a single
  /// integration point.
  double MagnetostaticEnergyFunction::call(const GeneralisedElement* ele_pt,
                                           MMInterpolator* intp) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Get the field
    Vector<double> h_ms;
    m_ele_pt->get_magnetostatic_field(intp->s(), h_ms);

    return -0.5 * VectorOps::dot(intp->m(), h_ms);
  }


  /// \short Calculate (dm/dt)^2 for the previous time step at a single
  /// integration point.
  double DmdtSquaredFunction::call(const GeneralisedElement* ele_pt,
                                   MMInterpolator* intp) const
  {
    return VectorOps::dot(intp->dmdt(), intp->dmdt());
  }

  /// \short Function for checking if integration is working ok, should
  /// give the area.
  double UnitFunction::call(const GeneralisedElement* ele_pt,
                            MMInterpolator* intp) const
  {
    return 1.0;
  }

}
