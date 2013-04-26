// Floating point error checks
#include <fenv.h>

#include "generic.h"

// meshes
#include "meshes/simple_rectangular_tri_mesh.h"
#include "meshes/tetgen_mesh.h"
#include "meshes/triangle_mesh.h"

#include "../../implicit_llg_problem.h"



using namespace oomph;
using namespace MathematicalConstants;
using namespace StringConversion;

namespace Inputs
{

  void initial_m_2d(const double& t, const Vector<double> &x, Vector<double> &m)
  {
    m.assign(3,0.0);

    double l = 5.0;
    m[0] = sin(x[0]*2*Pi/l) + sin(x[1]*2*Pi/l);
    m[1] = cos(x[0]*2*Pi/l) + cos(x[1]*2*Pi/l);
    m[2] = 1.0 - m[0] - m[1];

    // m[0] = x[0]/2 + x[1]/2;
    // m[1] = (1 - m[0]);
    // m[2] = 0.0;

    // m[0] = x[0];
    // m[1] = 0;
    // m[2] = 1 - x[0];

    VectorOps::normalise(m);
  }

  void initial_m_3d(const double& t, const Vector<double> &x, Vector<double> &m)
  {
    m.assign(3,0.0);

    double l = 5;
    m[0] = sin(x[0]*2*Pi*l) + sin(x[1]*2*Pi*l) + sin(x[2]*2*Pi*l);
    m[1] = cos(x[0]*2*Pi*l) + cos(x[1]*2*Pi*l) + cos(x[2]*2*Pi*l);
    m[2] = 1.0 - m[0] - m[1];

    // m[0] = x[0]/2 + x[1]/2;
    // m[1] = (1 - m[0]);
    // m[2] = 0.0;

    VectorOps::normalise(m);
  }

  // Applied fields
  void applied_field(const double& t, const Vector<double> &x,
                     Vector<double> &h_app)
  {
    h_app.assign(3,0.0);
    h_app[2] = 1;
  }

  void x_field(const double& t, const Vector<double> &x,
                     Vector<double> &h_app)
    {
      h_app.assign(3, 0.0);
      h_app[0] = 1;
    }

  void y_field(const double& t, const Vector<double> &x,
                     Vector<double> &h_app)
    {
      h_app.assign(3, 0.0);
      h_app[1] = 1;
    }

  void z_field(const double& t, const Vector<double> &x,
                     Vector<double> &h_app)
    {
      h_app.assign(3, 0.0);
      h_app[2] = 1;
    }


}


int main(int argc, char** argv)
{

  // Start MPI
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::init(argc,argv);
#endif

  // Default values
  unsigned nx(20), ny(20), refines(1); //refinement for different grid types
  double dt(1e-4), tmax(100);
  std::string outdir("results"), prec(""), mesh_type("");
  std::string field("default");

  double eps = 0;

  // Default to adaptive BDF2
  bool use_midpoint = false;

  // amg params
  unsigned amg_smoother_type(9);
  unsigned amg_v_cycles(1);

  // Set up and read command line arguments
  CommandLineArgs::setup(argc,argv);
  CommandLineArgs::specify_command_line_flag("-nx", &nx);
  CommandLineArgs::specify_command_line_flag("-ny", &ny);
  CommandLineArgs::specify_command_line_flag("-r", &refines);
  CommandLineArgs::specify_command_line_flag("-dt", &dt);
  CommandLineArgs::specify_command_line_flag("-tmax", &tmax);
  CommandLineArgs::specify_command_line_flag("-outdir", &outdir);
  CommandLineArgs::specify_command_line_flag("-prec", &prec);
  CommandLineArgs::specify_command_line_flag("-amgsmth", &amg_smoother_type);
  CommandLineArgs::specify_command_line_flag("-amgvcyc", &amg_v_cycles);
  CommandLineArgs::specify_command_line_flag("-mesh", &mesh_type);
  CommandLineArgs::specify_command_line_flag("-field", &field);

  CommandLineArgs::parse_and_assign();
  CommandLineArgs::output();

  bool adaptive_flag = (eps !=0 );

  // Enable some floating point error checkers
  //  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  // Create problem
  const unsigned dim = 2;
  ImplicitLLGProblem<TMicromagElement<dim,2> > problem;

  // Create timestepper
  TimeStepper* ts_pt = 0;
  if(use_midpoint)
    {
      ts_pt = new BDF<1>(adaptive_flag);
    }
  else
    {
      ts_pt = new BDF<2>(adaptive_flag);
    }
  problem.add_time_stepper_pt(ts_pt);

  // Create mesh
  Mesh* mesh_pt(0);
  double lx = 1.0, ly = lx;
  if((mesh_type == "") || (mesh_type == "structured_rect"))
    {
      if(dim != 2) { std::cerr << "WRONG DIM." <<std::endl; return 10; }
      mesh_pt = new SimpleRectangularTriMesh<TMicromagElement<dim,2> >
        (nx,ny,lx,ly,ts_pt);
    }
  else if(mesh_type == "sphere")
    {
      if(dim != 3) { std::cerr << "WRONG DIM." <<std::endl; return 10; }
      mesh_pt = new TetgenMesh<TMicromagElement<dim,2> >
        ("sphere."+to_string(refines)+".node",
         "sphere."+to_string(refines)+".ele",
         "sphere."+to_string(refines)+".poly",
         ts_pt);
    }
  else if(mesh_type == "unstructured_rect")
    {
      if(dim != 2) { std::cerr << "WRONG DIM." <<std::endl; return 10; }
      mesh_pt = new TriangleMesh<TMicromagElement<dim,2> >
        ("square." + to_string(refines) + ".node",
         "square." + to_string(refines) + ".ele",
         "square." + to_string(refines) + ".poly",
         ts_pt);
    }

  mesh_pt->setup_boundary_element_info();
  problem.bulk_mesh_pt() = mesh_pt;

  // Set magnetic params, fiddle with exchange strength
  // problem.mag_parameters_pt()->set_mumag4();
  // problem.mag_parameters_pt()->exchange_constant() = mag_parameters::mu0/5;
  // // problem.mag_parameters_pt()->gilbert_damping() = 0.1;
  problem.mag_parameters_pt()->set_simple_llg_parameters();

  // Set applied field
  if(field == "x")
    problem.applied_field_fct_pt() = &Inputs::x_field;
  else if(field == "y")
    problem.applied_field_fct_pt() = &Inputs::y_field;
  else if(field == "z")
    problem.applied_field_fct_pt() = &Inputs::z_field;
  else if(field == "default")
    problem.applied_field_fct_pt() = &Inputs::applied_field;

  // Finished setup, now we can build the problem
  problem.build();

  // If we specify a preconditioner then use GMRES (id for
  // non-preconditioned).
  if(prec != "")
    {

      // Set up linear solver: GMRES
      problem.linear_solver_pt() = new GMRES<CRDoubleMatrix>;
      GMRES<CRDoubleMatrix>* it_solver_pt =
        dynamic_cast< GMRES<CRDoubleMatrix>* >(problem.linear_solver_pt());
      it_solver_pt->tolerance() = 1e-8;
      it_solver_pt->set_preconditioner_RHS();
      it_solver_pt->open_convergence_history_file_stream(outdir+"/gmres");

      // Set max iter low to speed up tests
      it_solver_pt->max_iter() = 150;

      // Set up a preconditioner
      if(prec == "amg")
        {
#ifdef OOMPH_HAS_HYPRE
          /// AMG preconditioner + options
          HyprePreconditioner* amg_pt = new HyprePreconditioner;
          it_solver_pt->preconditioner_pt() = amg_pt;
          amg_pt->hypre_method() = HyprePreconditioner::BoomerAMG;
          amg_pt->set_amg_iterations(amg_v_cycles);

          amg_pt->enable_hypre_error_messages();

          //??ds complex and simple smoother type numbers overlap at 6. In
          // fact this whole thing of having seperate enumerations and a
          // boolean is stupid... and why isn't it a real enumeration
          // instead of an int?

          // Set the smoother type
          if(amg_smoother_type < 6)
            {
              amg_pt->amg_using_simple_smoothing();
              amg_pt->amg_simple_smoother() = amg_smoother_type;
            }
          else if(amg_smoother_type >= 6)
            {
              amg_pt->amg_using_complex_smoothing();
              amg_pt->amg_complex_smoother() = amg_smoother_type;
            }

          // If we are using Euclid then we have even more settings to play with:
          if(amg_smoother_type == 9)
            {
              amg_pt->AMGEuclidSmoother_level = 2;
              //amg_pt->AMGEuclidSmoother_print_level = 1;

              // amg_pt->AMGEuclidSmoother_use_block_jacobi = false;
              // amg_pt->AMGEuclidSmoother_use_row_scaling = false;
              // amg_pt->AMGEuclidSmoother_use_ilut = false;
              // amg_pt->AMGEuclidSmoother_drop_tol = 0;
            }
#else // If no Hypre then give an error
          throw OomphLibError("Don't have Hypre.",
                              OOMPH_CURRENT_FUNCTION,OOMPH_EXCEPTION_LOCATION);
#endif
        }
      else if(prec == "euclid")
        {
#ifdef OOMPH_HAS_HYPRE
          HyprePreconditioner* euclid_pt = new HyprePreconditioner;
          it_solver_pt->preconditioner_pt() = euclid_pt;
          euclid_pt->hypre_method() = HyprePreconditioner::Euclid;
#else // If no Hypre then give an error
          throw OomphLibError("Don't have Hypre.",
                              OOMPH_CURRENT_FUNCTION,OOMPH_EXCEPTION_LOCATION);
#endif
        }
      else if(prec == "ilu0")
        {
          /// ilu0 preconditioner
          it_solver_pt->preconditioner_pt() = new ILUZeroPreconditioner<CRDoubleMatrix>;
        }
      else if(prec == "exact")
        {
          // Exact preconditioner (superlu)
          it_solver_pt->preconditioner_pt() = new SuperLUPreconditioner;
        }
      else if(prec == "id")
        {
          // Identity precondtioner
          it_solver_pt->preconditioner_pt() = new IdentityPreconditioner;
        }
      else
        {
          std::cout << "Don't recognise the preconditioner type " << prec
                    << std::endl;
          exit(1);
        }
    }

  // Initialise problem
  problem.initialise_dt(dt);

  if(problem.dim() == 2)
    {
      problem.set_initial_condition(Inputs::initial_m_2d);
    }
  else if(problem.dim() == 3)
    {
      problem.set_initial_condition(Inputs::initial_m_3d);
    }


  // Set up outputs
  DocInfo doc_info;
  doc_info.set_directory(outdir);
  problem.doc_solution(doc_info);

  ConvergenceData conv_data;
  problem.convergence_data_pt() = &conv_data;
  conv_data.write_headers(outdir+"/convergence_data");

  // // Dump blocked Jacobian
  // DummyBlockPreconditioner<CRDoubleMatrix> blocker;
  // blocker.add_mesh(problem.bulk_mesh_pt());

  // // put phi blocks at end so we can ignore them... ??ds hacky!
  // Vector<unsigned> block_ordering(5,3);
  // block_ordering[2] = 0;
  // block_ordering[3] = 1;
  // block_ordering[4] = 2;
  // blocker.set_dof_to_block_map(block_ordering);

  // blocker.Preconditioner::setup(&problem, &J);
  // blocker.output_blocks_to_files(outdir+"/j_t0");


  // Code for adjusting solver based on dt
  double iterative_solver_threshold = 0;// 0.01;
  // gmres with Euclid smoothed amg seems very efficient below dt~0.01
  SuperLUSolver spare_LUSolver;
  LinearSolver* spare_iterative_solver_pt = problem.linear_solver_pt();


  // Timestep
  if(adaptive_flag)
    {
      // Timestep to end
      while(problem.time() < tmax)
        {
          std::cout << "step number = " << doc_info.number()
                    << ", time = " << problem.time()
                    << ", dt = " << dt << std::endl;
          // Step
          if(dt < iterative_solver_threshold)
            {
              problem.linear_solver_pt() = spare_iterative_solver_pt;
            }
          else
            {
              problem.linear_solver_pt() = &spare_LUSolver;
            }

          dt = problem.adaptive_unsteady_newton_solve(dt,eps);

          // Output
          // conv_data.output_this_newton_step(outdir+"/convergence_data");
        }
    }
  else
    {
      unsigned nstep = int(tmax/dt);
      for(unsigned t=0; t < nstep; t++)
        {
          std::cout << "Timestep " << t << std::endl;

          // Step
          problem.unsteady_newton_solve(dt);

          // Output
          problem.doc_solution(doc_info);
        }
    }


  // Dump Jacobian (we do this after newton solves so that we only get here
  // if the newton method can converge).
  DoubleVector dummy;
  CRDoubleMatrix J;
  problem.get_jacobian(dummy, J);
  J.sparse_indexed_output(outdir+"/jacobian_tmax");


  // if(prec != "")
  //   {
  //     // Output iteration counts
  //     std::ofstream trace;
  //     trace.open("iteration_counts",std::ios::app);
  //     trace << nx
  //           << "," << dt
  //           << "," << amg_smoother_type
  //           << "," << amg_v_cycles
  //           << "," << conv_data.average_newton_iterations()
  //           << "," << conv_data.average_linear_solver_iterations() << std::endl;
  //     trace.close();
  //   }

#ifdef OOMPH_HAS_MPI
  // Shut down oomph-lib's MPI
  MPI_Helpers::finalize();
#endif

  return 0;
}
