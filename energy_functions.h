#ifndef OOMPH_ENERGY_FUNCTIONS_H
#define OOMPH_ENERGY_FUNCTIONS_H

/*
  Functions that can be integrated over an element. Mostly for
  micromagnetics energy calculations so far.

  Anything derived from ElementalFunction can be integrated by the function
  integrate_over_element(..) currently in MicromagEquations (hopefully will
  be moved to Element or Finite Element).
*/

#include "../../src/generic/Vector.h"
#include "../../src/generic/elements.h"

#include "micromag_types.h"

namespace oomph
{

  class MMInterpolator;

  // Boilerplate junk...

  class ExchangeEnergyFunction : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp_pt) const;
  };

  class dExchangeEnergydtFunction : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp_pt) const;
  };

  class ZeemanEnergyFunction : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp_pt) const;
  };

  class dZeemanEnergydtFunction : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp_pt) const;
  };

  class CrystallineAnisotropyEnergyFunction  : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp_pt) const;
  };

  class dCrystallineAnisotropydtEnergyFunction  : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp_pt) const;
  };

  class MagnetostaticEnergyFunction : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp_pt) const;
  };

  class dMagnetostaticEnergydtFunction : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp_pt) const;
  };

  class DmdtSquaredFunction : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp_pt) const;
  };

  class UnitFunction : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp_pt) const;
  };

  class DmdtDotHeff : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp_pt) const;
  };

  class ExactFunctionDiffSquared : public ElementalFunction
  {
    public:
    InitialMFct* Exact_pt;
    double call(const GeneralisedElement* ele_pt, MMInterpolator* intp_pt) const;
  };



} // End of oomph namespace

#endif
