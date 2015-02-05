#include "renormalisation_handler.h"
#include "llg_problem.h"

namespace oomph {
  void AlwaysRenormalise::renormalise(LLGProblem* problem_pt)
  {
    problem_pt->force_renormalise_magnetisation();
  }

}
