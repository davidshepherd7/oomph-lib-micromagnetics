/*
  description of file goes here
*/

#include "generic.h"

#include "../../my_general_header.h"
#include "../../semi_implicit_problem.h"
#include "../../midpoint_method.h"
#include "../../magnetics_helpers.h"

#include "../../micromag.h"

// Floating point error checks
#include <fenv.h>

using namespace oomph;
using namespace MathematicalConstants;
using namespace StringConversion;

int main(int argc, char *argv[])
{
  // Start MPI if necessary
  MPI_Helpers::init(argc,argv);

#ifdef PARANOID
  // Enable some floating point error checkers
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
#endif

  // Read and process command line arguments
  MyCliArgs args;
  args.parse(argc, argv);

  TimeStepper* ts_pt = args.time_stepper_pt;
  double lx = 1.0;
  unsigned nx = 5 * std::pow(2, args.refinement);

  // Create problem
  SemiImplicitHybridMicromagneticsProblem problem;

    // problem(new SimpleRectangularQuadMesh<QMagnetostaticFieldElement<2,2> >
    //         (nx, nx, lx, lx, ts_pt),
    //         new SimpleRectangularQuadMesh<QMagnetostaticFieldElement<2,2> >
    //         (nx, nx, lx, lx, ts_pt),
    //         new SimpleRectangularQuadMesh<QSemiImplicitMicromagElement<2,2> >
    //         (nx, nx, lx, lx, ts_pt));


  Mesh* phi_1_mesh_pt = new SimpleRectangularQuadMesh<QMagnetostaticFieldElement<2,2> >
    (nx, nx, lx, lx, ts_pt);
  Mesh* phi_mesh_pt = new SimpleRectangularQuadMesh<QMagnetostaticFieldElement<2,2> >
    (nx, nx, lx, lx, ts_pt);
  Mesh* llg_mesh_pt = new SimpleRectangularQuadMesh<QSemiImplicitMicromagElement<2,2> >
    (nx, nx, lx, lx, ts_pt);

  problem.add_time_stepper_pt(args.time_stepper_pt);
  problem.set_bulk_mesh_pt(llg_mesh_pt);

  GenericPoissonProblem* phi_1_problem_pt = new GenericPoissonProblem;
  phi_1_problem_pt->set_bulk_mesh(phi_1_mesh_pt);
  phi_1_problem_pt->Flux_mesh_factory = &Factories::surface_mesh_factory<QMagnetostaticFieldFluxElement<2,2> >;
  problem.set_phi_1_problem_pt(phi_1_problem_pt);

  GenericPoissonProblem* phi_problem_pt = new GenericPoissonProblem;
  phi_problem_pt->set_bulk_mesh(phi_mesh_pt);
  problem.set_phi_problem_pt(phi_problem_pt);


  // Set up phi_1 problem
  // ============================================================
  {
    bool pin_phi1 = true;

    Vector<unsigned> neu_bound;
    for(unsigned b=0, nb=phi_1_mesh_pt->nboundary(); b < nb; b++)
      {
        neu_bound.push_back(b);
      }

    // phi_1 b.c.s are all Neumann but with flux determined elsewhere (by m)
    phi_1_problem_pt->set_neumann_boundaries(neu_bound, 0);

    if(pin_phi1)
      {
        // Pin a node which isn't involved in the boundary element method (we
        // have to pin something to avoid a singular Jacobian, can't be a
        // boundary node or things will go wrong with BEM).
        Node* pinned_phi_1_node_pt = phi_1_mesh_pt->get_some_non_boundary_node();
        pinned_phi_1_node_pt->pin(0);
        pinned_phi_1_node_pt->set_value(0,0.0);
      }
    else
      {
        std::cout << "Warning: not pinning phi1 at any point, technically J is singular..."
                  << " you might be ok..."
                  << std::endl;
      }

  // Finish off the problem
  phi_1_problem_pt->build();

  // Things will go wrong if the nodes of all meshes are not in
  // the same place:
#ifdef PARANOID
  if((phi_1_mesh_pt->nnode() != phi_mesh_pt->nnode())
     || (phi_1_mesh_pt->nnode() != llg_mesh_pt->nnode()))
    {
      std::ostringstream error_msg;
      error_msg << "Mesh nodes must be the same.";
      throw OomphLibError(error_msg.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

  unsigned dim = phi_mesh_pt->node_pt(0)->ndim();
  for(unsigned j=0; j<dim; j++)
    {
      for(unsigned nd=0, nnode= phi_1_mesh_pt->nnode(); nd<nnode; nd++)
        {
          if((phi_1_mesh_pt->node_pt(nd)->x(j) !=
              phi_mesh_pt->node_pt(nd)->x(j))
             ||
             (phi_1_mesh_pt->node_pt(nd)->x(j) !=
              llg_mesh_pt->node_pt(nd)->x(j)))
            {
              std::ostringstream error_msg;
              error_msg << "Mesh nodes must be in the same places.";
              throw OomphLibError(error_msg.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
        }
    }
#endif

  // Assign micromagnetics element pointers
  // ??ds dodgy...
  for(unsigned e=0, ne=phi_1_mesh_pt->nelement(); e < ne; e++)
    {
      MagnetostaticFieldEquations* ele_pt = checked_dynamic_cast<MagnetostaticFieldEquations*>
        (phi_1_mesh_pt->element_pt(e));

      MicromagEquations* m_ele_pt = checked_dynamic_cast<MicromagEquations*>
        (llg_mesh_pt->element_pt(e));

      ele_pt->set_micromag_element_pt( m_ele_pt);
    }
  }


  // BEM handler:
  // ============================================================
  BoundaryElementHandler<MicromagFaceElement<QSemiImplicitMicromagElement<2,2> > >*
    bem_handler_pt = new   BoundaryElementHandler<MicromagFaceElement<QSemiImplicitMicromagElement<2,2> > >
;
  {
    // Construct the BEM (must be done before pinning phi values)
    bem_handler_pt->set_bem_all_boundaries(phi_1_mesh_pt);
    // both zero because they are in seperate problems
    bem_handler_pt->input_index() = 0;
    bem_handler_pt->output_index() = 0;

    // Create an integration scheme
    bem_handler_pt->integration_scheme_pt() =
      Factories::variable_order_integrator_factory(problem.phi_1_mesh_pt()->finite_element_pt(0));

    bem_handler_pt->input_corner_data_pt() = 0; //??Ds

    bem_handler_pt->build();
  }
  problem.bem_handler_pt() = bem_handler_pt;


  // Do the rest (mag parameters, phi etc.)
  // ============================================================

  // // Assign timestepper and mesh from input arguments.
  // problem.add_time_stepper_pt(args.time_stepper_pt);
  // problem.bulk_mesh_pt() = args.mesh_pt;
  // problem.bulk_mesh_pt()->setup_boundary_element_info();

  // Set magnetic parameters
  double gilbert_damping = 0.5, hk = 0.0; // (normalised hk)
  problem.mag_parameters_pt()->set_simple_llg_parameters();
  problem.mag_parameters_pt()->gilbert_damping() = gilbert_damping;
  // problem.mag_parameters_pt()->exchange_constant() = 0.0;
  problem.mag_parameters_pt()->magnetostatic_debug_coeff() = 0.0;
  problem.mag_parameters_pt()->k1() = hk * mag_parameters::mu0/2;

  // Set applied field
  problem.applied_field_fct_pt() = args.h_app_fct_pt;

  // Set exact solution if we have one
  if((args.h_app_name == "minus_z") && (args.initial_m_name == "z"))
    {
      problem.Compare_with_mallinson = true;
    }

  // Finished setup, now we can build the problem
  problem.build();

  // Initialise problem and output
  problem.initialise_dt(args.dt);
  problem.set_initial_condition(args.initial_m_fct_pt);

  problem.Doc_info.set_directory(args.outdir);
  problem.Doc_info.Args_pt = &args;
  problem.Doc_info.output_jacobian = args.output_jacobian;

  problem.initial_doc();

  // All ready: step until completion
  double dt = args.dt;
  while(problem.time() < args.tmax)
    {
      std::cout << "step number = " << problem.Doc_info.number()
                << ", time = " << problem.time()
                << ", dt = " << dt
                << ", |m| error = " << 1 - problem.mean_nodal_magnetisation_length()
                << std::endl;

      // The Newton step itself, adaptive if requested
      if(args.adaptive_flag())
        {
          dt = problem.adaptive_unsteady_newton_solve(dt, args.tol);
        }
      else
        {
          problem.semi_implicit_step(dt);
        }

      // Output
      problem.doc_solution();
    }


  problem.final_doc();

  // Shut down oomph-lib's MPI
  MPI_Helpers::finalize();

  return 0;
}
