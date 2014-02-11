#include "my_generic_problem.h"
#include "magnetics_helpers.h"
#include "energy_functions.h"

#include "../../src/meshes/tetgen_mesh.h"

namespace oomph
{

  namespace MManipulation
  {

    /// \short Compute the effective damping constant (alpha) for the
    /// previous time step (see Albuquerque2001).
    double alt_effective_damping_used(MyProblem* problem_pt,
                                      std::deque<double>& previous_energies)
    {
      // Integral over all space of (dm/dt)^2 used in last step
      double dmdt_squared = integral_of_dmdt_squared(problem_pt); //??ds

      // If no change then damping is undefined
      if(dmdt_squared  == 0) return nan("");

      // Forumla from Albuquerque2001 & dAquino2005
      double dEdt = alt_dEnergydt(problem_pt, previous_energies);
      double effective_alpha = - dEdt / dmdt_squared;

      return effective_alpha;
    }


    /// \short Compute the effective damping constant (alpha) for the
    /// previous time step (see Albuquerque2001).
    double effective_damping_used(MyProblem* problem_pt)
    {
      // Integral over all space of (dm/dt)^2 used in last step
      double dmdt_squared = integral_of_dmdt_squared(problem_pt);

      // If no change then damping is undefined
      if(dmdt_squared  == 0) return nan("");

      // Forumla from Albuquerque2001 & dAquino2005
      double dEdt = dEnergydt(problem_pt);
      double effective_alpha = - dEdt / dmdt_squared;

      return effective_alpha;
    }


    double exchange_energy(MyProblem* problem_pt)
    {
      ExchangeEnergyFunction f;
      return problem_pt->integrate_over_problem(&f);
    }


    double zeeman_energy(MyProblem* problem_pt)
    {
      ZeemanEnergyFunction f;
      return problem_pt->integrate_over_problem(&f);
    }

    double crystalline_anisotropy_energy(MyProblem* problem_pt)
    {
      CrystallineAnisotropyEnergyFunction f;
      return problem_pt->integrate_over_problem(&f);
    }


    double magnetostatic_energy(MyProblem* problem_pt)
    {
      MagnetostaticEnergyFunction f;
      return problem_pt->integrate_over_problem(&f);
    }

    double integral_of_dmdt_squared(MyProblem* problem_pt)
    {
      DmdtSquaredFunction f;
      return problem_pt->integrate_over_problem(&f);
    }

    double dEnergydt(MyProblem* problem_pt)
    {
      dExchangeEnergydtFunction de_exdt;
      double I_de_exdt = problem_pt->integrate_over_problem(&de_exdt);

      dZeemanEnergydtFunction de_zeedt;
      double I_de_zeedt = problem_pt->integrate_over_problem(&de_zeedt);

      dCrystallineAnisotropydtEnergyFunction de_cadt;
      double I_de_cadt = problem_pt->integrate_over_problem(&de_cadt);

      dMagnetostaticEnergydtFunction de_ms;
      double I_de_ms = problem_pt->integrate_over_problem(&de_ms);

      return I_de_exdt + I_de_zeedt + I_de_cadt + I_de_ms;
    }

    double alt_dEnergydt(MyProblem* problem_pt,
std::deque<double>& previous_energies)
    {
      // Make a BDF2 time stepper to look up weights from (because I'm
      // lazy...)
      BDF<2> bdf;
      TimeStepper* node_ts_pt = problem_pt->mesh_pt()->
        finite_element_pt(0)->node_pt(0)->time_stepper_pt();
      bdf.time_pt() = node_ts_pt->time_pt();
      bdf.set_weights();

      // Calculate first derivative
      double deriv = 0.0;
      for(unsigned t=0;t<bdf.ntstorage();t++)
        {
          deriv += bdf.weight(1,t) * previous_energies[t];
        }

      return deriv;
    }
  }

  namespace MeshCreationHelpers
  {
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
                                  OOMPH_EXCEPTION_LOCATION,
                                  OOMPH_CURRENT_FUNCTION);
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
              double J = tele_pt->dshape_eulerian(s,psi,dpsidx);
              J++; // suppress unused variable warning
#endif
            }


        }


      // Use TetMeshBase to sort out boundary stuff
      out_mesh.set_nboundary(brick_mesh.nboundary());
      out_mesh.setup_boundary_element_info();

    }

  }
}
