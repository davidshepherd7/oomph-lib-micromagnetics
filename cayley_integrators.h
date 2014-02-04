#ifndef OOMPH_CAYLEY_INTEGRATORS_H
#define OOMPH_CAYLEY_INTEGRATORS_H

#include "../../src/generic/double_vector.h"
#include "../../src/generic/Vector.h"
#include "../../src/generic/matrices.h"
#include "../../src/generic/oomph_utilities.h"

#include "../../src/generic/explicit_timesteppers.h"


namespace oomph
{

  class CayleyEuler : public ExplicitTimeStepper
  {
  public:

    CayleyEuler() {}

    void timestep(ExplicitTimeSteppableObject* const &object_pt,
                  const double &dt);

  private:

    /// Disabled copy/assign
    CayleyEuler(const CayleyEuler&) {}
    void operator=(const CayleyEuler&) {}

  };


  class CayleyRK2 : public ExplicitTimeStepper
  {
  public:

    CayleyRK2() {}

    void timestep(ExplicitTimeSteppableObject* const &object_pt,
                  const double &dt);


  private:

    /// Disabled copy/assign
    CayleyRK2(const CayleyRK2&) {}
    void operator=(const CayleyRK2&) {}

  };

}

#endif
