
#include "generic.h"
#include "micromag.h"

#include "../../../../src/meshes/simple_cubic_mesh.h"
#include "../../../../src/meshes/simple_cubic_tet_mesh.h"


using namespace oomph;

namespace oomph
{


  void make_problem(const std::string& tag,
                    Mesh& mesh, GenericPoissonProblem& problem)
  {
    Steady<0>* ts_pt = new Steady<0>;
    problem.add_time_stepper_pt(ts_pt);

    problem.Doc_info.set_directory("Validation/" + tag);

    // Set the factory function used to create the surface mesh for
    // Neumman boundaries.
    if(dynamic_cast<QElementGeometricBase*>(mesh.element_pt(0)) != 0)
      {
        problem.set_flux_mesh_factory(&Factories::
                                      poisson_surface_mesh_factory
                                      <QTFPoissonFluxElement<3,2> >);
      }
    else if(dynamic_cast<TElementGeometricBase*>(mesh.element_pt(0)) != 0)
      {
        problem.set_flux_mesh_factory(&Factories::
                                      poisson_surface_mesh_factory
                                      <TTFPoissonFluxElement<3,2> >);
      }
    else
      {
        std::string err = "Unknown element geom.";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }

    // Assign b.c.s
    for(unsigned b=0; b<mesh.nboundary(); b++)
      {
        if(b != 1)
          {
            problem.add_dirichlet_boundary(&mesh, b,
                                           &TanhSolnForPoisson::get_exact_u);
          }
      }

    problem.add_neumann_boundary
      (&mesh, 1, &TanhSolnForPoisson::prescribed_flux_on_fixed_x_boundary);

    // Assign function pointers
    problem.set_source_fct_pt(&TanhSolnForPoisson::source_function);
    problem.exact_solution_fct_pt() = &TanhSolnForPoisson::get_exact_u;

    // Mesh vector for build function
    Vector<Mesh*> meshes(1);
    meshes[0] = &mesh;

    // build
    problem.build(meshes);
  }


}


int main()
{

  Steady<0> ts;

  // // Find out node positions for tets
  // SimpleCubicTetMesh<TTFPoissonElement<3, 2> > tm(5,5,5, 2,2,2, &ts);
  // const FiniteElement* ele_pt = tm.finite_element_pt(0);
  // for(unsigned nd=0, nnd=ele_pt->nnode(); nd<nnd; nd++)
  //   {
  //     Node* nd_pt = ele_pt->node_pt(nd);
  //     Vector<double> vec(3, 0.0);
  //     nd_pt->position(vec);
  //     std::cout << vec << std::endl;
  //   }
  // std::cout << "fin" << std::endl;


  // Build meshes
  SimpleCubicMesh<QTFPoissonElement<3, 2> > brick_mesh(10,10,10, 1.0,1.0,1.0, &ts);
  TetMeshBase* tet_mesh_pt = new TetMeshBase;
  MeshCreationHelpers::brick2tet(brick_mesh,
                                 MeshCreationHelpers::new_element<TTFPoissonElement<3, 2> >,
                                 *tet_mesh_pt);


  // Check boundaries ok
  if(brick_mesh.nboundary() != tet_mesh_pt->nboundary())
    {
      return 1;
    }
  for(unsigned b=0; b<brick_mesh.nboundary(); b++)
    {
      if(brick_mesh.nboundary_node(b) != tet_mesh_pt->nboundary_node(b))
        {
          return 2;
        }
    }


  // Set the orientation of the "step" to 45 degrees
  TanhSolnForPoisson::TanPhi=1.0;

  // Initial value for the steepness of the "step"
  TanhSolnForPoisson::Alpha=1.0;



  // Build + run with bricks
  GenericPoissonProblem brick_problem;
  make_problem("brick", brick_mesh, brick_problem);
  brick_problem.initial_doc();
  brick_problem.newton_solve();
  brick_problem.doc_solution();
  brick_problem.final_doc();
  std::cout << brick_problem.get_error_norm() << std::endl;

  // Build + run with tets
  GenericPoissonProblem tet_problem;
  make_problem("tet", *tet_mesh_pt, tet_problem);
  tet_problem.initial_doc();
  tet_problem.newton_solve();
  tet_problem.doc_solution();
  tet_problem.final_doc();
  std::cout << tet_problem.get_error_norm() << std::endl;

  // Dump to compare solution at nodes (hard to compare elsewhere since
  // elements not same).
  brick_mesh.dump("Validation/brick.dump", false);
  tet_mesh_pt->dump("Validation/tet.dump", false);

  return 0;
}
