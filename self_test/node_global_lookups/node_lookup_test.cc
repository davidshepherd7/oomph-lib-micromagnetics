
// Floating point debugging
#include <fenv.h>

#include "generic.h"
#include "../../poisson_test_problem.h"

#include "../../boundary_element_handler.h"

using namespace oomph;


// =================================================================
/// Added lookup test
// =================================================================
int added_lookup_test()
{
  // Build a test problem
  GenericPoissonForTests test_poisson;

  // Make a mesh of just the unpinned boundary nodes
  Mesh test_mesh;
  unsigned pdof = 0;
  // Loop over all meshes in problem
  for(unsigned msh=0, nmsh=test_poisson.nsub_mesh(); msh<nmsh; msh++)
    {
      for(unsigned b=0; b <test_poisson.mesh_pt(msh)->nboundary(); b++)
        {
          for(unsigned nd=0, nnode=test_poisson.mesh_pt(msh)->nboundary_node(b); nd<nnode; nd++)
            {
              Node* nd_pt = test_poisson.mesh_pt(msh)->boundary_node_pt(b,nd);
              if(!(nd_pt->is_pinned(pdof)))
                {
                  test_mesh.add_node_pt(nd_pt);
                }
            }
        }
    }


  // Build lookup
  AddedMainNumberingLookup lookup(&test_mesh, pdof);

  // std::cout << *(lookup.added_to_main_mapping_pt()) << std::endl;
  // std::cout << std::endl;
  // std::cout << *(lookup.main_to_added_mapping_pt()) << std::endl;

  // Check lookup works
  for(unsigned nd=0, nnode=test_mesh.nnode(); nd<nnode; nd++)
    {
      Node* nd_pt = test_mesh.node_pt(nd);

      unsigned main_from_added = nd_pt->eqn_number(pdof);
      unsigned main_from_lookup = lookup.added_to_main(nd);
      if(main_from_lookup != main_from_added)
        {
          std::cout << "Forward (added to main) lookup failed"<< std::endl;
          std::cout << lookup.added_to_main(nd) << " vs "
                    << main_from_added << std::endl;
          return 1;
        }

      int local_from_lookup = lookup.main_to_added(main_from_added);
      if((local_from_lookup < 0) || (unsigned(local_from_lookup) != nd))
        {
          std::cout << "Reverse (main to added) lookup failed" << std::endl;
          std::cout << lookup.main_to_added(main_from_added) << " vs "
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

  int result = added_lookup_test();

  if(result == 0)
    {
      std::cout << "***" <<std::endl;
      std::cout << "*** added_lookup_test() passed" << std::endl;
      std::cout << "***" <<std::endl;
    }

  return result;
}
