#include "my_generic_problem.h"
#include "magnetics_helpers.h"
#include "energy_functions.h"

#include "llg_problem.h"

#include "../../src/meshes/tetgen_mesh.h"

namespace oomph
{

  namespace InitialM
  {
    /// Wave-like exact solution from Jeong2014 and various other papers.
    Vector<double> LLGWaveSolution::operator()(const double& t,
                                               const Vector<double>& x) const
    {
      using namespace MathematicalConstants;
      using namespace std;

#ifdef PARANOID
      if(x.size() != dim)
        {
          std::string err = "x vector is the wrong length.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      // solution parameters
      double a = Pi/12;
      double k = 2 * Pi;

      // Rescale time because this is a solution to the LL equation
      double t_scaled = t / (1 + damping*damping);


      double sum_x = 0.0;
      for(unsigned j=0; j<dim; j++) {sum_x += x[j];}

      Vector<double> m(5, 0.0);
      if(damping == 0.0)
        {
          m[2] = sin(a) * cos(k*sum_x + dim*k*k*cos(a)*t_scaled);
          m[3] = sin(a) * sin(k*sum_x + dim*k*k*cos(a)*t_scaled);
          m[4] = cos(a);
        }
      else
        {
          double b = dim*k*k*damping*t_scaled;
          double d = sqrt(sin(a)*sin(a)  +  exp(2*b)*cos(a)*cos(a));
          double g = (1/damping) * log((d  +  exp(b)*cos(a))
                                       /(1  +  cos(a)));

          m[2] = (1/d) * sin(a) * cos(k*sum_x + g);
          m[3] = (1/d) * sin(a) * sin(k*sum_x + g);
          m[4] = (1/d) * cos(a) * exp(b);
        }

      return m;
    }

    void LLGWaveSolution::initialise_from_problem(const Problem* problem_pt)
    {
      const LLGProblem* llg_pt
        = checked_dynamic_cast<const LLGProblem*>(problem_pt);
      dim = llg_pt->dim();
      damping = llg_pt->mag_parameters_pt()->Gilbert_damping;
    }
  }

  namespace MManipulation
  {

    /// \short Compute the effective damping constant (alpha) for the
    /// previous time step (see Albuquerque2001).
    double alt_effective_damping_used(const LLGProblem& problem,
                                      std::deque<double>& previous_energies)
    {
      // Integral over all space of (dm/dt)^2 used in last step
      double dmdt_squared = integral_of_dmdt_squared(problem); //??ds

      // If no change then damping is undefined
      if(dmdt_squared  == 0) return nan("");

      // Forumla from Albuquerque2001 & dAquino2005
      double dEdt = alt_dEnergydt(problem, previous_energies);
      double effective_alpha = - dEdt / dmdt_squared;

      return effective_alpha;
    }


    /// \short Compute the effective damping constant (alpha) for the
    /// previous time step (see Albuquerque2001).
    double effective_damping_used(const LLGProblem& problem)
    {
      // Integral over all space of (dm/dt)^2 used in last step
      double dmdt_squared = integral_of_dmdt_squared(problem);

      // If no change then damping is undefined
      if(dmdt_squared  == 0) return nan("");

      // Forumla from Albuquerque2001 & dAquino2005
      double dEdt = dEnergydt(problem);
      double effective_alpha = - dEdt / dmdt_squared;

      return effective_alpha;
    }

    double effective_damping_used_3(const LLGProblem& problem)
    {
      // Integral over all space of (dm/dt)^2 used in last step, use the
      // default integration scheme (reduced if we are using it, otherwise
      // Gaussian).
      double dmdt_squared = integral_of_dmdt_squared(problem, 0);

      // If no change then damping is undefined
      if(dmdt_squared  == 0) return nan("");

      const double dt=  problem.time_pt()->dt();
      double dEdt = (problem.Previous_energies[0]
                     - problem.Previous_energies[1])/dt;
      double effective_alpha = - dEdt / dmdt_squared;

      return effective_alpha;
    }


    double exchange_energy(const LLGProblem& problem,
                           const Integral* quadrature_pt)
    {
      ExchangeEnergyFunction f;
      return problem.integrate_over_problem(&f, quadrature_pt);
    }


    double zeeman_energy(const LLGProblem& problem,
                         const Integral* quadrature_pt)
    {
      ZeemanEnergyFunction f;
      return problem.integrate_over_problem(&f, quadrature_pt);
    }

    double crystalline_anisotropy_energy(const LLGProblem& problem,
                                         const Integral* quadrature_pt)
    {
      CrystallineAnisotropyEnergyFunction f;
      return problem.integrate_over_problem(&f, quadrature_pt);
    }


    double magnetostatic_energy(const LLGProblem& problem,
                                const Integral* quadrature_pt)
    {
      MagnetostaticEnergyFunction f;
      return problem.integrate_over_problem(&f, quadrature_pt);
    }

    double integral_of_dmdt_squared(const LLGProblem& problem,
                                    const Integral* quadrature_pt)
    {
      DmdtSquaredFunction f;
      return problem.integrate_over_problem(&f, quadrature_pt);
    }

    double dEnergydt(const LLGProblem& problem)
    {
      dExchangeEnergydtFunction de_exdt;
      double I_de_exdt = problem.integrate_over_problem(&de_exdt);

      dZeemanEnergydtFunction de_zeedt;
      double I_de_zeedt = problem.integrate_over_problem(&de_zeedt);

      dCrystallineAnisotropydtEnergyFunction de_cadt;
      double I_de_cadt = problem.integrate_over_problem(&de_cadt);

      dMagnetostaticEnergydtFunction de_ms;
      double I_de_ms = problem.integrate_over_problem(&de_ms);

      return I_de_exdt + I_de_zeedt + I_de_cadt + I_de_ms;
    }

    double alt_dEnergydt(const LLGProblem& problem,
                         const std::deque<double>& previous_energies)
    {
      // // Make a BDF2 time stepper to look up weights from (because I'm
      // // lazy...)
      // BDF<2> bdf;
      // TimeStepper* node_ts_pt = problem.mesh_pt()->
      //   finite_element_pt(0)->node_pt(0)->time_stepper_pt();
      // bdf.time_pt() = node_ts_pt->time_pt();
      // bdf.set_weights();

      // // Calculate first derivative
      // double deriv = 0.0;
      // for(unsigned t=0;t<bdf.ntstorage();t++)
      //   {
      //     deriv += bdf.weight(1,t) * previous_energies[t];
      //   }
#warning "Not working now"
      return -1;
    }
  }

  namespace MeshCreationHelpers
  {


    unsigned min_node_number(Vector<Node*> nodes,
                             std::map<Node*, unsigned>& node_number_lookup)
    {
      Vector<unsigned> node_numbers;
      for(unsigned j=0; j<nodes.size(); j++)
        {
          unsigned node_num = node_number_lookup.find(nodes[j])->second;
          node_numbers.push_back(node_num);
        }

      return *std::min_element(node_numbers.begin(),
                               node_numbers.end());
    }


    void brick2tet(const Mesh& brick_mesh,
                   ElementFactoryFctPt element_factory_fpt,
                   TetMeshBase& out_mesh)
    {
      std::string err = "Disabled because can't work with private Boundary_node_pt.";
      throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                          OOMPH_CURRENT_FUNCTION);
      // Could use add function but would be very slow. Disable for now
      // since not in use...

// #ifdef PARANOID
//       if(brick_mesh.finite_element_pt(0)->dim() != 3)
//         {
//           std::string err = "Only for bricks! (3D)";
//           throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
//                               OOMPH_CURRENT_FUNCTION);
//         }
//       if(brick_mesh.finite_element_pt(0)->nnode() != 8
//          || brick_mesh.finite_element_pt(0)->nnode_1d() != 2)
//         {
//           std::string err = "Only for bricks w/ nnode1d = 2! (nnode = 8)";
//           throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
//                               OOMPH_CURRENT_FUNCTION);
//         }
//       if(brick_mesh.finite_element_pt(0)->nnode() != 8)
//         {
//           std::string err = "Only for bricks w/ nnode1d = 2! (nnode = 8)";
//           throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
//                               OOMPH_CURRENT_FUNCTION);
//         }
// #endif

//       // Copy number of boundaries over first
//       out_mesh.set_nboundary(brick_mesh.nboundary());


//       // From "How to Subdivide Pyramids, Prisms and Hexahedra into
//       // Tetrahedra" by Julien Dompierre Paul Labbé Marie-Gabrielle Vallet
//       // Ricardo Camarero: Can subdivide by creating tets with these nodes
//       // from the brick.

//       // Note to go from their ordering to ours

//       // 1) convert One-based -> Zero-based
//       // {{0, 1, 2, 5},
//       //  {0, 2, 7, 5},
//       //  {0, 2, 3, 7},
//       //  {0, 5, 7, 4},
//       //  {3, 7, 5, 6}};
//       // In emacs we can do this automatically with:
//       // M-x query-replace-regexp
//       // \([0-9]+\)
//       // \,(- \#1 1)

//       // 2) 2 <-> 3 and 6 <->7 to convert between brick node orderings:
//       // {{0, 1, 3, 5},
//       //  {0, 3, 6, 5},
//       //  {0, 3, 2, 6},
//       //  {0, 5, 6, 4},
//       //  {2, 6, 5, 7}};

//       // unsigned brick2tet_map[4][6][4] = {

//       //   // Zero diags through "VI7"
//       //   {{0, 1, 3, 5}, // Nodes in tet one
//       //    {0, 3, 6, 5}, // Nodes in tet two
//       //    {0, 3, 2, 6}, // etc...
//       //    {0, 5, 6, 4},
//       //    {2, 6, 5, 7},
//       //    {0, 0, 0, 0}},

//       //   // One
//       //   {{0, 5, 6, 4},
//       //    {0, 1, 6, 5},
//       //    {1, 7, 6, 5},
//       //    {0, 6, 3, 2},
//       //    {0, 6, 1, 3},
//       //    {1, 6, 7, 3}},

//       //   // two
//       //   {{0, 4, 5, 7},
//       //    {0, 2, 6, 7},
//       //    {0, 6, 4, 7},
//       //    {0, 1, 3, 5},
//       //    {0, 2, 7, 3},
//       //    {0, 7, 5, 3}},

//       //   // three
//       //   {{0, 3, 2, 7},
//       //    {0, 2, 6, 7},
//       //    {0, 6, 4, 7},
//       //    {0, 5, 7, 4},
//       //    {1, 5, 7, 0},
//       //    {1, 7, 3, 0}}
//       // };

//       unsigned brick2tet_map[4][6][4] = {

//         // Zero diags through "VI7"
//         {{0, 3, 1, 5}, // Nodes in tet one
//          {0, 6, 3, 5}, // Nodes in tet two
//          {0, 2, 3, 6}, // etc...
//          {0, 6, 5, 4},
//          {2, 5, 6, 7},
//          {0, 0, 0, 0}},

//         // One
//         {{0, 6, 5, 4},
//          {0, 6, 1, 5},
//          {1, 6, 7, 5},
//          {0, 3, 6, 2},
//          {0, 1, 6, 3},
//          {1, 7, 6, 3}},

//         // two
//         {{0, 5, 4, 7},
//          {0, 6, 2, 7},
//          {0, 4, 6, 7},
//          {0, 3, 1, 5},
//          {0, 7, 2, 3},
//          {0, 5, 7, 3}},

//         // three
//         {{0, 2, 3, 7},
//          {0, 6, 2, 7},
//          {0, 4, 6, 7},
//          {0, 7, 5, 4},
//          {1, 7, 5, 0},
//          {1, 3, 7, 0}}
//       };


//       // oomph-lib face indicies are weird
//       int faces[6] = {-3, -2, -1, 1, 2, 3};


//       // Copy nodes directly into mesh, add to boundaries if needed
//       for(unsigned nd=0, nnd=brick_mesh.nnode(); nd<nnd; nd++)
//         {
//           Node* nd_pt = brick_mesh.node_pt(nd);
//           out_mesh.add_node_pt(nd_pt);
//         }

//       // For each node: if on a boundary then loop over boundaries adding to
//       // the mesh's boundary node lists. ??ds move into mesh?
//       for(unsigned nd=0, nnd=out_mesh.nnode(); nd<nnd; nd++)
//         {
//           Node* nd_pt = out_mesh.node_pt(nd);

//           if(nd_pt->is_on_boundary())
//             {
//               std::set<unsigned>* boundaries_pt;
//               nd_pt->get_boundaries_pt(boundaries_pt);

//               std::set<unsigned>::const_iterator it;
//               for(it=boundaries_pt->begin(); it!=boundaries_pt->end(); it++)
//                 {
//                   out_mesh.Boundary_node_pt[*it].push_back(nd_pt);
//                 }
//             }
//         }


//       // Create a reverse lookup ready for element creation
//       std::map<Node*, unsigned> node_number_lookup;
//       for(unsigned nd=0, nnd=brick_mesh.nnode(); nd<nnd; nd++)
//         {
//           Node* nd_pt = brick_mesh.node_pt(nd);
//           node_number_lookup[nd_pt] = nd;
//         }


//       // Create new elements
//       for(unsigned ele=0, nele=brick_mesh.nelement(); ele<nele; ele++)
//         {
//           const FiniteElement* qele_pt = brick_mesh.finite_element_pt(ele);

// #ifdef PARANOID
//           if(dynamic_cast<const FaceElement*>(qele_pt) != 0)
//             {
//               throw OomphLibError("Function not implemented for face elements",
//                                   OOMPH_EXCEPTION_LOCATION,
//                                   OOMPH_CURRENT_FUNCTION);
//             }
// #endif

//           // vertex I 7 in the paper is equivalent to our node 7, nodes
//           // diagonally opposite to it are equivalent to our nodes 1,2,6.
//           Vector<Node*> VI7_nodes;
//           VI7_nodes.push_back(qele_pt->node_pt(7));
//           VI7_nodes.push_back(qele_pt->node_pt(1));
//           VI7_nodes.push_back(qele_pt->node_pt(2));
//           VI7_nodes.push_back(qele_pt->node_pt(6));


//           unsigned ndiag = 0;
//           for(unsigned face_i=0; face_i<6; face_i++)
//             {

//               unsigned face = faces[face_i];

//               // Get list of nodes in face
//               Vector<Node*> Face_nodes;
//               for(unsigned j=0; j<4; j++)
//                 {
//                   unsigned nd = qele_pt->get_bulk_node_number(face, j);
//                   Face_nodes.push_back(qele_pt->node_pt(nd));
//                 }

//               // Get the first node in the face (by the brick mesh ordering).
//               unsigned first_node_num = min_node_number(Face_nodes,
//                                                         node_number_lookup);
//               Node* first_node_pt = brick_mesh.node_pt(first_node_num);

//               bool through_VI7 = (std::find(VI7_nodes.begin(), VI7_nodes.end(),
//                                             first_node_pt)
//                                   != VI7_nodes.end());

//               if(through_VI7) {ndiag++;}
//             }


//           // Create new elements and assign nodes
// #ifdef PARANOID
//           double tvol = 0.0;
// #endif
//           for(unsigned j=0; j<6; j++)
//             {
//               // if ndiag == 0 we only need 5 elements so skip the last one
//               if(ndiag == 0 && j != 6) {continue;}

//               // Create the new T element and assign nodes
//               FiniteElement* tele_pt = element_factory_fpt();
//               tele_pt->node_pt(0) = qele_pt->node_pt(brick2tet_map[ndiag][j][0]);
//               tele_pt->node_pt(1) = qele_pt->node_pt(brick2tet_map[ndiag][j][1]);
//               tele_pt->node_pt(2) = qele_pt->node_pt(brick2tet_map[ndiag][j][2]);
//               tele_pt->node_pt(3) = qele_pt->node_pt(brick2tet_map[ndiag][j][3]);

// #ifdef PARANOID
//               tvol += tele_pt->size();
// #endif

//               // Put it into the element list
//               out_mesh.add_element_pt(tele_pt);
//             }

// #ifdef PARANOID
//           // Check that volumes are unchanged
//           double qvol = qele_pt->size();
//           if(std::abs(tvol - qvol) > 1e-8)
//             {
//               std::string err = "Element volume changed by more than 1e-8";
//               throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
//                                   OOMPH_CURRENT_FUNCTION);
//             }
// #endif
//         }


//       // Use TetMeshBase to sort out boundary element stuff
//       out_mesh.setup_boundary_element_info();



// #ifdef PARANOID
//       // Check jacobians are not inverted
//       Shape psi(4);
//       DShape dpsidx(4, 3);
//       Vector<double> s(3, 0.3);
//       for(unsigned ele=0, nele=out_mesh.nelement(); ele<nele; ele++)
//         {
//           double J = out_mesh.finite_element_pt(ele)->dshape_eulerian(s,psi,dpsidx);
//           J++; // suppress unused variable warning
//         }
// #endif

    }

    /// File-private helper function to check that nodal coords match in
    /// all dimensions except one.
    bool periodic_coords_match(const Node* nd1_pt, const Node* nd2_pt,
                               const unsigned& dir)
    {
#ifdef PARANOID
      if(nd1_pt->ndim() != nd2_pt->ndim())
        {
          std::string err = "Nodes have different dimensions!";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif

      const unsigned dim = nd1_pt->ndim();
      for(unsigned j=0; j<dim; j++)
        {
          if(j != dir)
            {
              if(std::abs(nd1_pt->x(j) - nd2_pt->x(j)) > 1e-10)
                {
                  return false;
                }
            }
        }

      return true;
    }

    /// File-private helper function to tell if two nodes are linked via
    /// copying.
    bool copied_equivalents(Node* nd1_pt, Node* nd2_pt)
    {
      return
        (   (nd1_pt == nd2_pt)
            || (nd1_pt->copied_node_pt() == nd2_pt)
            || (nd2_pt->copied_node_pt() == nd1_pt)
            || ((nd1_pt->copied_node_pt() != 0)
                && (nd1_pt->copied_node_pt() == nd2_pt->copied_node_pt()))
            );
    }


    /// Helper function to make two boundaries of a simple mesh periodic.
    void make_boundaries_periodic(Mesh* mesh_pt, const unsigned& b1,
                                  const unsigned& b2,
                                  const unsigned& direction)
    {
      const unsigned nbn = mesh_pt->nboundary_node(b1);

#ifdef PARANOID
      if(nbn != mesh_pt->nboundary_node(b2))
        {
          std::string err = "Boundaries must have same number of nodes";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif

      Vector<std::pair<Node*, Node*> > dodgy_node_pairs;

      // Handle case where nodes are in the opposite order on opposite
      // boundaries. Note that we can't handle more complex variations of
      // ordering here, but these cases will get caught by the paranoid
      // check below.
      bool opposite_direction_flag =
        periodic_coords_match(mesh_pt->boundary_node_pt(b1, 0),
                              mesh_pt->boundary_node_pt(b2, nbn-1),
                              direction);

      // Loop over nodes setting up copys
      for(unsigned nd=0; nd<nbn; nd++)
        {
          Node* nd1_pt = mesh_pt->boundary_node_pt(b1, nd);

          // Second node index depends on flag
          unsigned nd2_index = opposite_direction_flag ? (nbn-nd-1) : nd;
          Node* nd2_pt = mesh_pt->boundary_node_pt(b2, nd2_index);

          // Check the coordinates look ok for them to be periodic
#ifdef PARANOID
          if(!periodic_coords_match(nd1_pt, nd2_pt, direction))
            {
          std::string err = "Nodal positions don't match up.";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
            OOMPH_CURRENT_FUNCTION);
        }
#endif

          // Catch the case that they are already copies of each other or
          // both copies of the same node. This is ok, just leave
          // everything as it is.
          if(copied_equivalents(nd1_pt, nd2_pt))
            {
              // do nothing
            }

          // If both nodes are already copies then we need to do
          // something weird: link both the nodes to one of the copied
          // nodes and record pointers to both copied nodes so that we
          // can check this is consistent later.
          else if(nd1_pt->is_a_copy() && nd2_pt->is_a_copy())
            {
              Node* old_copied_pt = nd1_pt->copied_node_pt();

              nd1_pt->clear_copied_pointers();
              nd1_pt->make_periodic(nd2_pt);

              dodgy_node_pairs.push_back
                (std::make_pair(nd1_pt->copied_node_pt(), old_copied_pt));
            }

          // Otherwise it's easy: just make a node which isn't already a
          // copy into a copy.
          else if(nd1_pt->is_a_copy())
            {
              nd2_pt->make_periodic(nd1_pt);
            }
          else
            {
              nd1_pt->make_periodic(nd2_pt);
            }
        }


      // Finally check that the nodes we flagged as dodgy match up
      const unsigned ni = dodgy_node_pairs.size();
      for(unsigned i=0; i<ni; i++)
        {
          if(!copied_equivalents(dodgy_node_pairs[i].first,
                                 dodgy_node_pairs[i].second))
            {
              std::string err = "Inconsistent!";
              throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                  OOMPH_CURRENT_FUNCTION);
            }
        }
    }

    /// Helper function to make two boundaries of a mesh periodic. This
    /// version uses a brute force search to find the appropriate nodes to
    /// link together so it should be robust but is O(N^2). In practice it
    /// seems to be un-noticably fast for meshes that I've used so far.
    void slow_make_boundaries_periodic(Mesh* mesh_pt, const unsigned& b1,
                                       const unsigned& b2,
                                       const unsigned& direction)
    {
      const unsigned nbn = mesh_pt->nboundary_node(b1);

#ifdef PARANOID
      if(nbn != mesh_pt->nboundary_node(b2))
        {
          std::string err = "Boundaries must have same number of nodes";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif

      Vector<std::pair<Node*, Node*> > dodgy_node_pairs;

      // Loop over nodes setting up copys
      for(unsigned nd=0; nd<nbn; nd++)
        {
          Node* nd1_pt = mesh_pt->boundary_node_pt(b1, nd);

          // Search for matching node on opposite boundary
          Node* nd2_pt = 0;
          bool success = false;
          for(unsigned nd2=0; nd2<nbn; nd2++)
            {
              nd2_pt = mesh_pt->boundary_node_pt(b2, nd2);
              if(periodic_coords_match(nd1_pt, nd2_pt, direction))
                 {
                   success = true;
                   break;
                 }
            }
          if(!success)
            {
              std::string err = "Couldn't find matching node.";
              throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }


          // Catch the case that they are already copies of each other or
          // both copies of the same node. This is ok, just leave
          // everything as it is.
          if(copied_equivalents(nd1_pt, nd2_pt))
            {
              // do nothing
            }

          // If both nodes are already copies then we need to do
          // something weird: link both the nodes to one of the copied
          // nodes and record pointers to both copied nodes so that we
          // can check this is consistent later.
          else if(nd1_pt->is_a_copy() && nd2_pt->is_a_copy())
            {
              Node* old_copied_pt = nd1_pt->copied_node_pt();

              nd1_pt->clear_copied_pointers();
              nd1_pt->make_periodic(nd2_pt);

              dodgy_node_pairs.push_back
                (std::make_pair(nd1_pt->copied_node_pt(), old_copied_pt));
            }

          // Otherwise it's easy: just make a node which isn't already a
          // copy into a copy.
          else if(nd1_pt->is_a_copy())
            {
              nd2_pt->make_periodic(nd1_pt);
            }
          else
            {
              nd1_pt->make_periodic(nd2_pt);
            }
        }


      // Finally check that the nodes we flagged as dodgy match up
      const unsigned ni = dodgy_node_pairs.size();
      for(unsigned i=0; i<ni; i++)
        {
          if(!copied_equivalents(dodgy_node_pairs[i].first,
                                 dodgy_node_pairs[i].second))
            {
              std::string err = "Inconsistent!";
              throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                  OOMPH_CURRENT_FUNCTION);
            }
        }
    }
  }
}
