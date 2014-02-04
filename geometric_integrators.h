#ifndef OOMPH_GEOMETRIC_INTEGRATORS_H
#define OOMPH_GEOMETRIC_INTEGRATORS_H

#include "../../src/generic/double_vector.h"
#include "../../src/generic/Vector.h"
#include "../../src/generic/matrices.h"
#include "../../src/generic/oomph_utilities.h"

#include "../../src/generic/explicit_timesteppers.h"


namespace oomph
{

  class GeomEuler : public ExplicitTimeStepper
  {
  public:

    GeomEuler() {}

    void timestep(ExplicitTimeSteppableObject* const &object_pt,
                  const double &dt);

  private:

    /// Disabled copy/assign
    GeomEuler(const GeomEuler&) {}
    void operator=(const GeomEuler&) {}

  };


  class GeomRK2 : public ExplicitTimeStepper
  {
  public:

    GeomRK2() {}

    void timestep(ExplicitTimeSteppableObject* const &object_pt,
                  const double &dt);


  private:

    /// Disabled copy/assign
    GeomRK2(const GeomRK2&) {}
    void operator=(const GeomRK2&) {}

  };

}

#endif
