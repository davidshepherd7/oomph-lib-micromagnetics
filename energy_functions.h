#ifndef OOMPH_ENERGY_FUNCTIONS_H
#define OOMPH_ENERGY_FUNCTIONS_H

/*
  Functions that can be integrated over an element. Mostly for
  micromagnetics energy calculations so far.
*/

#include "../../src/generic/Vector.h"
#include "../../src/generic/elements.h"

using namespace oomph;

namespace oomph
{

  /// \short Function class for functions to integrate over elements.
  class ElementalFunction
  {
  public:
    virtual double call(const GeneralisedElement* ele_pt,
                        const Vector<double> &s) const = 0;
  };


  class ExchangeEnergyFunction : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt,
                const Vector<double> &s) const;
  };


  class ZeemanEnergyFunction : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt,
                const Vector<double> &s) const;
  };


  class CrystallineAnisotropyEnergyFunction  : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt,
                const Vector<double> &s) const;
  };

  class MagnetostaticEnergyFunction : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt,
                const Vector<double> &s) const;
  };


  class DmdtSquaredFunction : public ElementalFunction
  {
    double call(const GeneralisedElement* ele_pt,
                const Vector<double> &s) const;
  };



} // End of oomph namespace

#endif
