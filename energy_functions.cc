
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
    return 0.5 * (VectorOps::dot(intp_pt->dmdx(0), intp_pt->dmdx(0)) +
                  VectorOps::dot(intp_pt->dmdx(1), intp_pt->dmdx(1)) +
                  VectorOps::dot(intp_pt->dmdx(2), intp_pt->dmdx(2)));
  }


  /// \short Calculate the time derivative of energy due to exchange at a
  /// point at a single integration point.
  double dExchangeEnergydtFunction::call(const GeneralisedElement* ele_pt,
                                      MMInterpolator* intp_pt) const
  {
    return (VectorOps::dot(intp_pt->dmdx(0), intp_pt->d2mdxdt(0)) +
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
    Vector<double> h_applied =
      m_ele_pt->get_applied_field(intp_pt->time(), intp_pt->x());

    return - VectorOps::dot(intp_pt->m(), h_applied);
  }

  /// \short Calculate time derivative of energy due to external applied
  /// field at a single integration point assuming constant external field.
  double dZeemanEnergydtFunction::call(const GeneralisedElement* ele_pt,
                                       MMInterpolator* intp_pt) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Get the field
    Vector<double> h_applied = m_ele_pt->get_applied_field(intp_pt->time(),
                                                           intp_pt->x());

    return - VectorOps::dot(intp_pt->dmdt(), h_applied);
  }


  /// \short Calculate energy due to magnetocrystalline anisotropy at a
  /// single integration point.
  double CrystallineAnisotropyEnergyFunction::
  call(const GeneralisedElement* ele_pt, MMInterpolator* intp_pt) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Get parameters
    double k1 = m_ele_pt->magnetic_parameters_pt()->normalised_hk();
    Vector<double> e = m_ele_pt->magnetic_parameters_pt()->easy_axis();

    return k1 * (1 - std::pow(VectorOps::dot(intp_pt->m(), e), 2));
  }


  /// \short Calculate time derivative of energy due to magnetocrystalline
  /// anisotropy at a single integration point.
  double dCrystallineAnisotropydtEnergyFunction::
  call(const GeneralisedElement* ele_pt, MMInterpolator* intp_pt) const
  {
    const MicromagEquations* m_ele_pt
      = checked_dynamic_cast<const MicromagEquations*>(ele_pt);

    // Get parameters
    double k1 = m_ele_pt->magnetic_parameters_pt()->normalised_hk();
    Vector<double> e = m_ele_pt->magnetic_parameters_pt()->
      nearest_easy_axis(intp_pt->m());

    return - k1 * VectorOps::dot(intp_pt->m(), e)
      * VectorOps::dot(intp_pt->dmdt(), e);
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
    m_ele_pt->get_magnetostatic_field(intp_pt->s(), h_ms);

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
    m_ele_pt->get_magnetostatic_field(intp_pt->s(), h_ms);

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

  double DmdtDotHeff::call(const GeneralisedElement* ele_pt,
                           MMInterpolator* intp_pt) const
  {
    const MicromagEquations* m_ele_pt
      = dynamic_cast<const MicromagEquations*>(ele_pt);

    // Some storage
    Vector<double> gradmi_dot_gradtest(3, 0.0);
    double val = 0.0;

    // Get fields, coeffs at this point
    Vector<double> h_ms, h_ca;
    m_ele_pt->get_magnetostatic_field(intp_pt->s(), h_ms);
    Vector<double> h_app = m_ele_pt->get_applied_field(intp_pt->time(),
                                                       intp_pt->x());
    m_ele_pt->get_H_cryst_anis_field(intp_pt->time(), intp_pt->x(),
                                     intp_pt->m(), h_ca);
    double exch_c = m_ele_pt->exchange_coeff();

    //??ds should probably have something with gamma in it in here?
    if(m_ele_pt->llg_precession_coeff() != 1)
      {
            std::string error_msg = "precession coeff != 1, don't know if this will work!";
            throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
      }

    // For each node in the element:
    for(unsigned l=0, nl=m_ele_pt->nnode(); l<nl; l++)
      {
        // Get the weird exchange term after integration by parts
        for(unsigned i=0; i<3; i++)
          for(unsigned j=0; j<m_ele_pt->nodal_dimension(); j++)
            gradmi_dot_gradtest[i] += intp_pt->dtestdx(l,j) * intp_pt->dmdx(i)[j];

        // Add contributions from this node
        val += (VectorOps::dot(intp_pt->dmdt(), h_ms)
                + VectorOps::dot(intp_pt->dmdt(), h_app)
                + VectorOps::dot(intp_pt->dmdt(), h_ca) ) * intp_pt->test(l)
          - exch_c * VectorOps::dot(intp_pt->dmdt(), gradmi_dot_gradtest);
      }

    return val;
  }

  InitialMFct Exact_fpt;

  double ExactFunctionDiffSquared::call(const GeneralisedElement* ele_pt,
                                        MMInterpolator* intp_pt) const
  {
    const MicromagEquations* m_ele_pt
      = dynamic_cast<const MicromagEquations*>(ele_pt);

    Vector<double> exact = (*Exact_pt)(intp_pt->time(), intp_pt->x());

    // Assume that m indices are contiguous and extract m from entire
    // solution.
    unsigned mi0 = m_ele_pt->m_index_micromag(0);
    Vector<double> exact_m; exact_m.assign(exact.begin()+mi0, exact.end());

    double diff = 0;
    const unsigned ni = exact_m.size();
    for(unsigned i=0; i<ni; i++)
      {
        diff += std::pow(intp_pt->m()[i] - exact_m[i], 2);
      }

    return diff;
  }



}
