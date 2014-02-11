
#include "generic.h"
#include "micromag.h"

#include "../../../../src/meshes/simple_cubic_mesh.h"
#include "../../../../src/meshes/simple_cubic_tet_mesh.h"


using namespace oomph;

namespace oomph
{
  /// The trivial factory function for an ELEMENT.
  template<class ELEMENT>
  inline  FiniteElement* new_element() { return new ELEMENT;}


  void brick2tet(const Mesh& brick_mesh,
                 ElementFactoryFctPt element_factory_fpt,
                 TetMeshBase& out_mesh)
  {
#ifdef PARANOID
    if(brick_mesh.finite_element_pt(0)->dim() != 3)
      {
        std::string err = "Only for bricks! (3D)";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
    if(brick_mesh.finite_element_pt(0)->nnode() != 8
       || brick_mesh.finite_element_pt(0)->nnode_1d() != 2)
      {
        std::string err = "Only for bricks w/ nnode1d = 2! (nnode = 8)";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
    if(brick_mesh.finite_element_pt(0)->nnode() != 8)
      {
        std::string err = "Only for bricks w/ nnode1d = 2! (nnode = 8)";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif

    // From "How to Subdivide Pyramids, Prisms and Hexahedra into
    // Tetrahedra" by Julien Dompierre Paul Labbé Marie-Gabrielle Vallet
    // Ricardo Camarero: Can subdivide by creating tets with these nodes
    // from the brick.

    // Note to go from their ordering to ours

    // 1) convert One-based -> Zero-based
    // {{0, 1, 2, 5},
    //  {0, 2, 7, 5},
    //  {0, 2, 3, 7},
    //  {0, 5, 7, 4},
    //  {3, 7, 5, 6}};

    // 2) 2 <-> 3 and 6 <->7 to convert between brick node orderings:
    // {{0, 1, 3, 5},
    //  {0, 3, 6, 5},
    //  {0, 3, 2, 6},
    //  {0, 5, 6, 4},
    //  {2, 6, 5, 7}};

    // 3) to convert between tet node orderings: swap nodes 1 and 2:
    unsigned brick2tet_map[][4] = {{0, 3, 1, 5},
                                   {0, 6, 3, 5},
                                   {0, 2, 3, 6},
                                   {0, 6, 5, 4},
                                   {2, 5, 6, 7}};


    // Copy nodes directly into mesh
    for(unsigned nd=0, nnd=brick_mesh.nnode(); nd<nnd; nd++)
      {
        Node* nd_pt = brick_mesh.node_pt(nd);
        out_mesh.add_node_pt(nd_pt);
      }


    // Create new elements
    for(unsigned ele=0, nele=brick_mesh.nelement(); ele<nele; ele++)
      {
        const FiniteElement* qele_pt = brick_mesh.finite_element_pt(ele);

        // for(unsigned nd=0, nnd=qele_pt->nnode(); nd<nnd; nd++)
        //   {
        //     Node* nd_pt = qele_pt->node_pt(nd);
        //     Vector<double> vec(3, 0.0);
        //     nd_pt->position(vec);
        //     std::cout << vec << std::endl;
        //   }
        // std::cout << std::endl;

#ifdef PARANOID
        if(dynamic_cast<const FaceElement*>(qele_pt) != 0)
          {
            throw OomphLibError("Function not implemented for face elements",
                                OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
          }
#endif

        // Create new elements and assign nodes
        for(unsigned j=0; j<5; j++)
          {
            FiniteElement* tele_pt = element_factory_fpt();

            tele_pt->node_pt(0) = qele_pt->node_pt(brick2tet_map[j][0]);
            tele_pt->node_pt(1) = qele_pt->node_pt(brick2tet_map[j][1]);
            tele_pt->node_pt(2) = qele_pt->node_pt(brick2tet_map[j][2]);
            tele_pt->node_pt(3) = qele_pt->node_pt(brick2tet_map[j][3]);

            out_mesh.add_element_pt(tele_pt);

#ifdef PARANOID
            // Check jacobian
            Shape psi(4);
            DShape dpsidx(4, 3);
            Vector<double> s(3, 0.3);
            const double J = tele_pt->dshape_eulerian(s,psi,dpsidx);
#endif
          }


      }


    // Use TetMeshBase to sort out boundary stuff
    out_mesh.set_nboundary(brick_mesh.nboundary());
    out_mesh.setup_boundary_element_info();

  }


  void make_problem(Mesh& mesh, GenericPoissonProblem& problem)
  {
    problem.Doc_info.set_directory("Validation");

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
  brick2tet(brick_mesh, new_element<TTFPoissonElement<3, 2> >, *tet_mesh_pt);



  // Set the orientation of the "step" to 45 degrees
  TanhSolnForPoisson::TanPhi=1.0;

  // Initial value for the steepness of the "step"
  TanhSolnForPoisson::Alpha=1.0;



  // Build + run with bricks
  GenericPoissonProblem brick_problem;
  make_problem(brick_mesh, brick_problem);
  brick_problem.newton_solve();
  brick_problem.final_doc();
  std::cout << brick_problem.get_error_norm() << std::endl;

  // Build + run with tets
  GenericPoissonProblem tet_problem;
  make_problem(*tet_mesh_pt, tet_problem);
  tet_problem.newton_solve();
  tet_problem.final_doc();
  std::cout << tet_problem.get_error_norm() << std::endl;


  // Dump to compare solution at nodes (hard to compare elsewhere since
  // elements not same).
  brick_mesh.dump("Validation/brick.dump", false);
  tet_mesh_pt->dump("Validation/tet.dump", false);

  return 0;
}
