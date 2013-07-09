#ifndef OOMPH_ENERGY_FUNCTIONS_H
#define OOMPH_ENERGY_FUNCTIONS_H

/*
  Functions that can be integrated over an element. Mostly for
  micromagnetics energy calculations so far.
*/

#include "../../src/generic/Vector.h"
#include "../../src/generic/elements.h"

namespace oomph
{

  class MMInterpolator;

  /// \short Function class for functions to integrate over elements.
  class ElementalFunction
  {
  public:

    /// \short Virtual destructor (to stop compiler complaining).
    virtual ~ElementalFunction() {}

    /// \short Function to integrate over the element.
    virtual double call(const GeneralisedElement* ele_pt,
                        MMInterpolator* intp) const = 0;

    /// \short Helper function to automatically create an interpolator
    /// object for call.
    double call(const GeneralisedElement* ele_pt, const Vector<double> &s) const;
  };

  // Boilerplate junk...

  class ExchangeEnergyFunction : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp) const;
  };

  class dExchangeEnergydtFunction : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp) const;
  };

  class ZeemanEnergyFunction : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp) const;
  };

  class dZeemanEnergydtFunction : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp) const;
  };

  class CrystallineAnisotropyEnergyFunction  : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp) const;
  };

  class dCrystallineAnisotropydtEnergyFunction  : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp) const;
  };

  class MagnetostaticEnergyFunction : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp) const;
  };

  class DmdtSquaredFunction : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp) const;
  };

  class UnitFunction : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp) const;
  };



} // End of oomph namespace

#endif
