#ifndef OOMPH_RESIDUAL_CALCULATOR_H
#define OOMPH_RESIDUAL_CALCULATOR_H

/*
  A class holding a single function which calculates the residual.
  This allows us to dynamically choose which residual function to use :)
*/

#include "../../src/generic/elements.h"
#include "../../src/generic/matrices.h"
#include "../../src/generic/Vector.h"


using namespace oomph;

namespace oomph
{

  // for micromag ones:
  class MicromagEquations;

  class ResidualCalculator
  {
  public:
    virtual ~ResidualCalculator(){};

    /// The residual function
    virtual void fill_in_generic_residual_contribution
    (const MicromagEquations* const ele_pt,
     Vector<double> &residuals, DenseMatrix<double> &jacobian,
     const unsigned& flag) const = 0;
  };

  class LLResidualCalculator : public ResidualCalculator
  {
  public:
    /// The residual function
    void fill_in_generic_residual_contribution
    (const MicromagEquations* const ele_pt,
     Vector<double> &residuals, DenseMatrix<double> &jacobian,
     const unsigned& flag) const;
  };

  class LLGResidualCalculator : public ResidualCalculator
  {
  public:
    /// The residual function
    void fill_in_generic_residual_contribution
    (const MicromagEquations* const ele_pt,
     Vector<double> &residuals, DenseMatrix<double> &jacobian,
     const unsigned& flag) const;
  };

} // End of oomph namespace

#endif
