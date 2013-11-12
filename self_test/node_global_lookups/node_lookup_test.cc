
// Floating point debugging
#include <fenv.h>

#include "generic.h"
#include "../../poisson_test_problem.h"

#include "../../boundary_element_handler.h"

using namespace oomph;


// =================================================================
/// Node lookup test
// =================================================================
int node_lookup_test()
{
  // Build a test problem
  GenericPoissonForTests test_poisson;

  // Make a mesh of just the unpinned boundary nodes
  Mesh test_mesh;
  unsigned pdof = 0;
  for(unsigned b=0; b <test_poisson.bulk_mesh_pt()->nboundary(); b++)
    {
      for(unsigned nd=0, nnode=test_poisson.bulk_mesh_pt()->nboundary_node(b); nd<nnode; nd++)
        {
          Node* nd_pt = test_poisson.bulk_mesh_pt()->boundary_node_pt(b,nd);
          if(!(nd_pt->is_pinned(pdof)))
            {
              test_mesh.add_node_pt(nd_pt);
            }
        }
    }

  // Build lookup
  NodeGlobalNumbersLookup nd_lookup(&test_mesh, pdof);

  // std::cout << *(nd_lookup.node_to_global_mapping_pt()) << std::endl;
  // std::cout << std::endl;
  // std::cout << *(nd_lookup.global_to_node_mapping_pt()) << std::endl;

  // Check lookup works
  for(unsigned nd=0, nnode=test_mesh.nnode(); nd<nnode; nd++)
    {
      Node* nd_pt = test_mesh.node_pt(nd);

      unsigned global_from_node = nd_pt->eqn_number(pdof);
      unsigned global_from_lookup = nd_lookup.node_to_global(nd);
      if(global_from_lookup != global_from_node)
        {
          std::cout << "Forward (node to global) lookup failed"<< std::endl;
          std::cout << nd_lookup.node_to_global(nd) << " vs "
                    << global_from_node << std::endl;
          return 1;
        }

      int local_from_lookup = nd_lookup.global_to_node(global_from_node);
      if((local_from_lookup < 0) || (unsigned(local_from_lookup) != nd))
        {
          std::cout << "Reverse (global to node) lookup failed" << std::endl;
          std::cout << nd_lookup.global_to_node(global_from_node) << " vs "
                    << nd << std::endl;
          return 2;
        }
    }

  return 0;
}

int main()
{
  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  int result = node_lookup_test();

  if(result == 0)
    {
      std::cout << "***" <<std::endl;
      std::cout << "*** node_lookup_test() passed" << std::endl;
      std::cout << "***" <<std::endl;
    }

  return result;
}
