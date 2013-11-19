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

  class LLGResidualCalculator : public ResidualCalculator
  {
  public:

    LLGResidualCalculator(bool _use_gilbert_form)
    {
      Use_gilbert_form = _use_gilbert_form;
    }

    /// The residual function. Pick which form to use based on the flag.
    void fill_in_generic_residual_contribution
    (const MicromagEquations* const ele_pt,
     Vector<double> &residuals, DenseMatrix<double> &jacobian,
     const unsigned& flag) const
      {
        if(use_gilbert_form())
          {
            llg_residual(ele_pt, residuals, jacobian, flag);
          }
        else
          {
            ll_residual(ele_pt, residuals, jacobian, flag);
          }
      }

    bool use_gilbert_form() const {return Use_gilbert_form;}

    void set_use_gilbert_form() {Use_gilbert_form = true;}
    void set_use_ll_form() {Use_gilbert_form = false;}

    private:

    /// Calculate residual using gilbert form
    void llg_residual(const MicromagEquations* const e_pt,
                 Vector<double> &residuals, DenseMatrix<double> &jacobian,
                 const unsigned& flag) const;

    /// Calculate residual using ll form
    void ll_residual(const MicromagEquations* const e_pt,
                     Vector<double> &residuals, DenseMatrix<double> &jacobian,
                     const unsigned& flag) const;

    bool Use_gilbert_form;

  };

} // End of oomph namespace

#endif
