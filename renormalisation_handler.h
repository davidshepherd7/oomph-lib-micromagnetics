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

    RenormaliseAfterTolerance()
    {
      // ??ds hard code it for now
      handler_pt = new AlwaysRenormalise();
      tol = 1e-2; // This is the default value from magpar (see
                  // magpar-0.9/src/llg/checkiterationllg.c line 179).
    }

    RenormalisationHandler* handler_pt;
    double tol;

    void renormalise(LLGProblem* problem_pt);

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
    else if(name == "after-tol")
      {
        return new RenormaliseAfterTolerance();
      }
    else if(name == "onto-sphere")
      {
        return new RenomaliseOntoSphere();
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
