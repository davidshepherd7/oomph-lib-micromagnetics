#ifndef OOMPH_RENORMALISATION_HANDLERS_H
#define OOMPH_RENORMALISATION_HANDLERS_H

#include "../../src/generic/oomph_utilities.h"

namespace oomph
{

  class LLGProblem;

  /// Pure virtual base
  class RenormalisationHandler
  {
  public:
    virtual void renormalise(LLGProblem* problem_pt) = 0;
  };


  class AlwaysRenormalise : public RenormalisationHandler
  {
  public:
    void renormalise(LLGProblem* problem_pt);
  };


  class NeverRenormalise : public RenormalisationHandler
  {
  public:
    void renormalise(LLGProblem* problem_pt)
    {
      // no-op
    }
  };


  class RenormaliseAfterTolerance : public RenormalisationHandler
  {
  public:

    void renormalise(LLGProblem* problem_pt)
    {
      throw OomphLibError("Not implemented (yet?).", OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
  };


  class RenomaliseOntoSphere : public RenormalisationHandler
  {
  public:
    void renormalise(LLGProblem* problem_pt)
    {
      throw OomphLibError("Not implemented (yet?).", OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
  };


  inline RenormalisationHandler*
  renormalisation_handler_factory(std::string name)
  {
    if(name == "always")
      {
        return new AlwaysRenormalise();
      }
    else if(name == "never")
      {
        return new NeverRenormalise();
      }
    else
      {
        std::string err = "Unrecognised renormaliation setting \"" + name +"\".";
        throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
  }


} // End of oomph namespace

#endif
