#include "renormalisation_handler.h"
#include "llg_problem.h"
#include "vector_helpers.h"


namespace oomph {
  void AlwaysRenormalise::renormalise(LLGProblem* problem_pt)
  {
    oomph_info << "Renormalising nodal magnetisations." << std::endl;

    // Loop over meshes and renormalise m at each node
    for(unsigned nd=0; nd<problem_pt->mesh_pt()->nnode(); nd++)
      {
        Node* nd_pt = problem_pt->mesh_pt()->node_pt(nd);

        // Get m vector
        Vector<double> m_values(3,0.0);
        for(unsigned j=0; j<3; j++) m_values[j] = nd_pt->value(problem_pt->m_index(j));

        // Normalise
        VectorOps::normalise(m_values);

        // Write m vector
        for(unsigned j=0; j<3; j++) nd_pt->set_value(problem_pt->m_index(j),
                                                     m_values[j]);
      }
  }

}
