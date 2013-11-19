#ifndef OOMPH_RESIDUAL_SWAPPING_EXPLICIT_TIMESTEPPER_H
#define OOMPH_RESIDUAL_SWAPPING_EXPLICIT_TIMESTEPPER_H

#include "../../src/generic/Vector.h"
#include "../../src/generic/mesh.h"
#include "../../src/generic/explicit_timesteppers.h"


using namespace oomph;

namespace oomph
{

  class ResidualSwappingExplicitTimestepper : public ExplicitTimeStepper
  {

  public:

    ResidualSwappingExplicitTimestepper()
    {
      // null pointers
      residual_pt = 0;
      underlying_time_stepper_pt = 0;
    }


    virtual ~ResidualSwappingExplicitTimestepper()
    {
      delete underlying_time_stepper_pt; underlying_time_stepper_pt = 0;
      // mesh pts will be handled by problem
    }


    LLGResidualCalculator* residual_pt;
    ExplicitTimeStepper* underlying_time_stepper_pt;


    /// Advance time step. Swap the residual, call the underlying explicit
    /// time stepper then swap the residuals back.
    void timestep(ExplicitTimeSteppableObject* const &object_pt,
                  const double &dt)
    {
      // Swap to explicit (ll) form
      residual_pt->set_use_ll_form();

      // timestep
      underlying_time_stepper_pt->timestep(object_pt, dt);

      // Swap to explicit (gilbert) form
      residual_pt->set_use_gilbert_form();
    }


  private:

    /// Inaccessible copy constructor
    ResidualSwappingExplicitTimestepper
    (const ResidualSwappingExplicitTimestepper&){};

    /// Inaccessible assignement operator
    void operator=(const ResidualSwappingExplicitTimestepper&){};
  };

} // End of oomph namespace

#endif
