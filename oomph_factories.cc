
#include "oomph_factories.h"

// for llg block prec, fix? ??ds
#include "llg_factories.h"
#include "llg_preconditioners.h"


// mesh factories
#include "../../src/generic/mesh.h"

// solver factories
#include "../../src/generic/linear_solver.h"
#include "../../src/generic/iterative_linear_solver.h"
#include "../../src/generic/sum_of_matrices.h"

#ifdef OOMPH_HAS_HYPRE
#include "../../src/generic/hypre_solver.h"
#endif

// timestepper factories
#include "../../src/generic/timesteppers.h"
#include "../../src/generic/midpoint_method.h"
#include "../../src/generic/explicit_timesteppers.h"
#include "tr.h"
#include "cayley_integrators.h"
#include "bdf2_with_predictor.h"


// Preconditioner factories
#include "../../src/generic/preconditioner.h"
#include "../../src/generic/general_purpose_preconditioners.h"
#include "../../src/generic/general_purpose_block_preconditioners.h"
#include "sum_of_matrices_preconditioner.h"


namespace oomph
{
  namespace Factories
  {



    /// Make an explicit time stepper
    ExplicitTimeStepper* explicit_time_stepper_factory(const std::string& ts_name)
    {
      if(ts_name == "rk4")
        {
          return new RungeKutta<4>;
        }
      else if(ts_name == "rk2")
        {
          return new RungeKutta<2>;
        }
      else if(ts_name == "lrk4")
        {
          return new LowStorageRungeKutta<4>;
        }
      else if (ts_name == "ebdf3")
        {
          return new EBDF3;
        }
      else if (ts_name == "euler")
        {
          return new Euler;
        }
      else if(ts_name == "cay-euler")
        {
          return new CayleyEuler;
        }
      else if(ts_name == "cay-rk2")
        {
          return new CayleyRK2;
        }
      else
        {
          // throw error
          unrecognised_name(ts_name, OOMPH_CURRENT_FUNCTION);
          return 0;
        }
    }


    /// Make a timestepper from an input argument. Assumption: this will be
    /// passed into a problem, which will delete the pointer when it's
    /// done.
    TimeStepper* time_stepper_factory(const std::string& ts_name,
                                      const std::string& mp_pred_name,
                                      const int& mp_update_pinned)
    {

      // Always make timestepper adaptive, we can control adaptivity by
      // calling adaptive or non adaptive newton solve.
      bool adaptive_flag = true;

      if(ts_name == "bdf1")
        {
          return new BDF<1>(adaptive_flag);
        }
      else if(ts_name == "bdf2")
        {
          return new BDF<2>(adaptive_flag);
        }
      else if(ts_name == "bdf2-pred")
        {
          BDF2Pred* bdf2_pt = new BDF2Pred(adaptive_flag);
          ExplicitTimeStepper* pred_pt = explicit_time_stepper_factory(mp_pred_name);
          bdf2_pt->set_predictor_pt(pred_pt);
          return bdf2_pt;
        }
      else if((ts_name == "midpoint") || (ts_name == "old-imr"))
        {
          MidpointMethod* mp_pt = new MidpointMethod(adaptive_flag);
          ExplicitTimeStepper* pred_pt = explicit_time_stepper_factory(mp_pred_name);
          mp_pt->set_predictor_pt(pred_pt);
          return mp_pt;
        }
      else if((ts_name == "midpoint-bdf") || (ts_name == "imr"))
        {
          MidpointMethodByBDF* mp_pt = new MidpointMethodByBDF(adaptive_flag);
          ExplicitTimeStepper* pred_pt = explicit_time_stepper_factory(mp_pred_name);
          mp_pt->set_predictor_pt(pred_pt);
          if(mp_update_pinned != -1)
            {
              mp_pt->Update_pinned = bool(mp_update_pinned);
            }
          return mp_pt;
        }
      else if(ts_name == "steady")
        {
          // 2 steps so that we have enough space to do reasonable time
          // derivative estimates in e.g. energy derivatives.
          return new Steady<3>;
        }
      else if(ts_name == "tr")
        {
          // 2 steps so that we have enough space to do reasonable time
          // derivative estimates in e.g. energy derivatives.
          return new TR(adaptive_flag);
        }
      else
        {
          return 0;
        }
    }


    /// Create a vector of meshes with the names given in mesh details and
    /// shifted as also specified in mesh_details.
    Vector<Mesh*> multimesh_factory(MeshFactoryFctPt underlying_factory,
                                    Vector<ShiftedMeshDetails>& mesh_details,
                                    int refinement_level,
                                    TimeStepper* time_stepper_pt,
                                    double scaling,
                                    double rotate_xy_angle,
                                    unsigned nnode1d)
    {
#ifdef PARANOID
      if(underlying_factory == 0)
        {
          std::string err = "Null underlying mesh factory pointer";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif

      Vector<Mesh*> mesh_pts;

      const unsigned nj = mesh_details.size();
      for(unsigned j=0; j<nj; j++)
        {
          // Build it without scaling
          Mesh* mesh_pt = underlying_factory(mesh_details[j].mesh_name,
                                             refinement_level,
                                             time_stepper_pt,
                                             1.0,
                                             rotate_xy_angle,
                                             nnode1d);

          // Shift it
          MeshCreationHelpers::shift_mesh(mesh_details[j].xshift,
                                          mesh_details[j].yshift,
                                          mesh_details[j].zshift,
                                          mesh_pt);

          // Scale it. We do this after shifting so that we are scaling the
          // entire multi-mesh, including gaps between meshes. This is much
          // more intuitive.
          mesh_pt->scale_mesh(scaling);

          // Add to list
          mesh_pts.push_back(mesh_pt);
        }

      return mesh_pts;
    }


    /// Create a pair of meshes near to each other (shifted along x).
    Vector<Mesh*> simple_multimesh_factory(MeshFactoryFctPt underlying_factory,
                                           const std::string& mesh_name,
                                           int refinement_level,
                                           TimeStepper* time_stepper_pt,
                                           double xshift,
                                           double scaling,
                                           unsigned nnode1d)
    {
      Vector<ShiftedMeshDetails> inputs(2);
      inputs[0].mesh_name = mesh_name;
      inputs[0].xshift = -xshift;


      inputs[1].mesh_name = mesh_name;
      inputs[1].xshift = +xshift;

      return multimesh_factory(underlying_factory,
                               inputs, refinement_level,
                               time_stepper_pt, scaling, 0.0, nnode1d);
    }


    /// Create lots of meshes near to each other (spaced out along x, y).
    Vector<Mesh*> simple_many_multimesh_factory
    (MeshFactoryFctPt underlying_factory, const std::string& mesh_name,
     int refinement_level, TimeStepper* time_stepper_pt,
     double xspacing, double yspacing, double scaling, unsigned nnode1d)
    {
      // Idea is to just make a square (or rect.) grid of meshes, spaced
      // according to xspacing and yspacing.

      // Decided number of meshes
      unsigned n_along_one_direction = 4;
      unsigned n_total = n_along_one_direction*n_along_one_direction;
      Vector<ShiftedMeshDetails> inputs(n_total);

      // Keep track of shift distances
      double basex = 0.0;
      double basey = 0.0;

      // Loop through setting the shift distances for each one
      for(unsigned j=0; j<n_along_one_direction; j++)
        {
          for(unsigned i=0; i<n_along_one_direction; i++)
            {
              // Store the mesh details
              inputs[j*n_along_one_direction+i].mesh_name = mesh_name;
              inputs[j*n_along_one_direction+i].xshift = basex;
              inputs[j*n_along_one_direction+i].yshift = basey;

              // Next column: move along y
              basey += yspacing;
            }

          // Next row: zero the y shift and move along x
          basex += xspacing;
          basey = 0.0;
        }

      // Construct the meshes using the general multimesh factory
      return multimesh_factory(underlying_factory,
                               inputs, refinement_level,
                               time_stepper_pt, scaling, 0.0, nnode1d);
    }


    /// Private (to .cc file) helper function, have to do this to deal with
    /// oomph-libs templated iterative solvers without a massive list of
    /// ifs.
    template<class T>
    LinearSolver* linear_solver_factory_helper(const std::string& _solver_name,
                                               const double& krylov_tol,
                                               const unsigned& max_iter,
                                               const bool& throw_on_convergence_fail)
    {
      const std::string solver_name = to_lower(_solver_name);

      LinearSolver* solver_pt;

      if(solver_name == "fdlu")
        {
          solver_pt = new FD_LU;
        }
      else if(solver_name == "superlu")
        {
          solver_pt = new SuperLUSolver;
        }
      else if(solver_name == "gmres-left-prec")
        {
          GMRES<T>* gmres_pt = new GMRES<T>;
          gmres_pt->set_preconditioner_LHS();
          solver_pt = gmres_pt;
        }
      else if(solver_name == "gmres")
        {
          GMRES<T>* gmres_pt = new GMRES<T>;
          gmres_pt->set_preconditioner_RHS();
          solver_pt = gmres_pt;
        }
      else if(solver_name == "cg")
        {
          solver_pt = new CG<T>;
        }
      else
        {
          std::string err("Unrecognised solver name ");
          err += solver_name;
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }


      IterativeLinearSolver* its_pt
        = dynamic_cast<IterativeLinearSolver*>(solver_pt);
      if(its_pt != 0)
        {
          its_pt->tolerance() = krylov_tol;

          its_pt->max_iter() = max_iter;

          // set error behaviour
          if(throw_on_convergence_fail)
            {
              its_pt->enable_error_after_max_iter();
            }
          else
            {
              its_pt->disable_error_after_max_iter();
            }
        }

      return solver_pt;
    }


    LinearSolver* linear_solver_factory(const std::string& solver_name,
                                        const std::string& _matrix_type,
                                        const double& krylov_tol,
                                        const unsigned& max_iter,
                                        const bool& throw_on_convergence_fail)
    {

      const std::string matrix_type = to_lower(_matrix_type);

      if(matrix_type == "cr")
        {
          return linear_solver_factory_helper<CRDoubleMatrix>(solver_name,
                                                              krylov_tol,
                                                              max_iter,
                                                              throw_on_convergence_fail);
        }
      else if(matrix_type == "som")
        {
          return linear_solver_factory_helper<SumOfMatrices>(solver_name,
                                                             krylov_tol,
                                                             max_iter,
                                                             throw_on_convergence_fail);
        }
      else
        {
          std::string err("Unrecognised matrix type ");
          err += matrix_type;
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }


    Preconditioner* preconditioner_factory(const std::string &_prec_name)
    {
      const std::string prec_name = to_lower(_prec_name);
      Preconditioner* prec_pt = 0;

      // First check if we want a sum of matrice preconditioner:
      if(has_prefix("som-main-", prec_name))
        {
          Preconditioner* ul_prec = preconditioner_factory
            (rest_of_name("som-main-", prec_name));

          MainMatrixOnlyPreconditioner* mm_prec_pt = new MainMatrixOnlyPreconditioner;
          mm_prec_pt->set_underlying_prec_pt(ul_prec);

          prec_pt = mm_prec_pt;
        }
      // Make a preconditioner which only acts on the main matrix and
      // diagonals of added matrices of a sum of matrices.
      else if(has_prefix("som-maindiag-", prec_name))
        {
          Preconditioner* ul_prec = preconditioner_factory
            (rest_of_name("som-maindiag-", prec_name));

          MainMatrixAndDiagsPreconditioner* mm_prec_pt
            = new MainMatrixAndDiagsPreconditioner;
          mm_prec_pt->set_underlying_prec_pt(ul_prec);

          prec_pt = mm_prec_pt;
        }

      // AMG with optimal poisson paramters
      // ============================================================
      else if(prec_name == "poisson-amg")
        {
#ifdef OOMPH_HAS_HYPRE
          HyprePreconditioner* amg_pt = new HyprePreconditioner;
          amg_pt->hypre_method() = HyprePreconditioner::BoomerAMG;

          // Use good Poisson settings for 3D (2D problems should be ok
          // with the same ones I hope...).
          Hypre_default_settings::set_defaults_for_3D_poisson_problem(amg_pt);

          prec_pt = amg_pt;

#else // If no Hypre then give a warning and use exact
          OomphLibWarning("Don't have Hypre, using exact preconditioner.",
                          OOMPH_CURRENT_FUNCTION,OOMPH_EXCEPTION_LOCATION);
          prec_pt = preconditioner_factory("exact");
#endif
        }

      // General purpose AMG (default parameters)
      // ============================================================
      else if(prec_name == "amg")
        {
#ifdef OOMPH_HAS_HYPRE
          HyprePreconditioner* amg_pt = new HyprePreconditioner;
          amg_pt->hypre_method() = HyprePreconditioner::BoomerAMG;
          prec_pt = amg_pt;

#else // If no Hypre then give a warning and use exact
          OomphLibWarning("Don't have Hypre, using exact preconditioner.",
                          OOMPH_CURRENT_FUNCTION,OOMPH_EXCEPTION_LOCATION);
          prec_pt = preconditioner_factory("exact");
#endif
        }
      else if(prec_name == "identity")
        { prec_pt = new IdentityPreconditioner; }

      else if(prec_name == "none")
        { prec_pt = 0; }

      else if(prec_name == "exact")
        { prec_pt = new SuperLUPreconditioner; }

      else if(prec_name == "blockexact")
        {
          prec_pt = new ExactBlockPreconditioner<CRDoubleMatrix>;
        }

      else if(prec_name == "blocklt")
        {
          BlockTriangularPreconditioner<CRDoubleMatrix>* bp_pt
            = new BlockTriangularPreconditioner<CRDoubleMatrix>;
          bp_pt->lower_triangular();
          prec_pt = bp_pt;
        }

      else if(prec_name == "blockut")
        {
          BlockTriangularPreconditioner<CRDoubleMatrix>* bp_pt
            = new BlockTriangularPreconditioner<CRDoubleMatrix>;
          bp_pt->upper_triangular();
          prec_pt = bp_pt;
       }

      else if(has_prefix("ilu-", prec_name))
        {
// ??ds make more robust and move to function!

#ifdef OOMPH_HAS_HYPRE
          HyprePreconditioner* hp_pt = new HyprePreconditioner;

          // Use Euclid (ILU)
          hp_pt->use_Euclid();

          // Use other algorithms
          // hp_pt->enable_euclid_rowscale();
          // hp_pt->enable_euclid_using_BJ() = 0.0;
          // hp_pt->euclid_using_ILUT();

          // Debugging info
          hp_pt->euclid_print_level() = 0;

          prec_pt = hp_pt;

          std::string postfix = split_string(prec_name, '-')[1];
          if(postfix == "high")
            {
              hp_pt->euclid_level() = 3;
              hp_pt->euclid_droptol() = 1e-5;
            }
          else
            {
              // Use level specified
              hp_pt->euclid_level() = atoi(postfix.c_str());

              // Simple parameters
              hp_pt->euclid_droptol() = 0.0;
            }


#else // If no Hypre can only do ilu 0
          if(prec_name == "ilu-0")
            {
              prec_pt = new ILUZeroPreconditioner<CRDoubleMatrix>;
            }
          else
            {
              throw OomphLibError("Don't have Hypre, can't do ilu non-zero fill in",
                                  OOMPH_CURRENT_FUNCTION,OOMPH_EXCEPTION_LOCATION);
            }
#endif
        }
      else
        {
          unrecognised_name(prec_name, OOMPH_CURRENT_FUNCTION);
        }

      return prec_pt;
    }


    /// \short Construct a list of times to output the full solution at
    /// based on command line input in label.
    Vector<double> doc_times_factory(const double& doc_interval,
                                     const double &t_max)
    {
      Vector<double> doc_times;

      if(doc_interval == 0)
        {
          // Do nothing: empty vector = output at every step.
        }

      // Otherwise we have a number giving the interval between doc
      // times.
      else
        {
          // Add an output time every "doc_interval" time units until we get
          // to t_max.
          double doc_t = 0.0;
          while(doc_t < t_max)
            {
              doc_times.push_back(doc_t);
              doc_t += doc_interval;
            }
        }

      return doc_times;
    }


    // void my_problem_factory_helper(TimeStepper* ts_pt,
    //                                ExplicitTimeStepper* expl_ts_pt,
    //                                LinearSolver* linear_solver_pt,
    //                                const double& newton_tol,
    //                                const double& error_norm_limit,
    //                                bool disable_explicit_solver_optimisations,
    //                                const std::string& outdir,
    //                                const std::string& args_info,
    //                                const std::string& output_jacobian,
    //                                const Vector<double>& doc_times,
    //                                MyProblem& new_problem)
    //   {
    //     // Assign time steppers to problem. At most one of these two will be a
    //     // real timestepper, the other will be null or a dummy.
    //     new_problem.add_time_stepper_pt(time_stepper_pt);
    //     new_problem.set_explicit_time_stepper_pt(explicit_time_stepper_pt);

    //     // Assign general purpose parameters
    //     new_problem.linear_solver_pt() = solver_pt;
    //     new_problem.newton_solver_tolerance() = newton_tol;
    //     // new_problem.Use_fd_jacobian = use_fd_jacobian; //??ds
    //     new_problem.Error_norm_limit = error_norm_limit;
    //     new_problem.Disable_explicit_solver_optimisations =
    //       disable_explicit_solver_optimisations;

    //     // Assign doc parameters
    //     new_problem.Doc_info.copy_args_string(args_info);
    //     new_problem.Doc_info.set_directory(outdir);
    //     new_problem.Doc_info.output_jacobian = output_jacobian;
    //     new_problem.set_doc_times(doc_times);
    //   }


    // LLGProblem* llg_problem_factory(LLGArgs* args_pt)
    // {
    //   LLGProblem* problem_pt = new LLGProblem;

    //   my_problem_factory_helper( *problem_pt)

    // }


    Integral* gauss_integration_factory(const unsigned& dim,
                                        const unsigned& nnode_1d,
                                        const ElementGeometry::ElementGeometry& geom)
    {
      using namespace ElementGeometry;

      if(nnode_1d == 2)
        {
          if(dim ==1 && geom == T) // T line
            {
              return new TGauss<1, 2>;
            }
          else if(dim ==1 && geom == Q) // Q line
            {
              return new Gauss<1, 2>;
            }
          else if((dim == 2) && geom == T) //triangle
            {
              return new TGauss<2, 2>;
            }
          else if((dim == 2) && geom == Q) // square
            {
              return new Gauss<2, 2>;
            }
          else if((dim == 3) && geom == T) // tet
            {
              return new TGauss<3, 2>;
            }
          else if((dim == 3) && geom == Q) // brick
            {
              return new Gauss<3, 2>;
            }
        }
      if(nnode_1d == 3)
        {
          if(dim ==1 && geom == T) // T line
            {
              return new TGauss<1, 3>;
            }
          else if(dim ==1 && geom == Q) // Q line
            {
              return new Gauss<1, 3>;
            }
          else if((dim == 2) && geom == T) //triangle
            {
              return new TGauss<2, 3>;
            }
          else if((dim == 2) && geom == Q) // square
            {
              return new Gauss<2, 3>;
            }
          else if((dim == 3) && geom == T) // tet
            {
              return new TGauss<3, 3>;
            }
          else if((dim == 3) && geom == Q) // brick
            {
              return new Gauss<3, 3>;
            }
        }
      else if(nnode_1d == 4)
        {
          if(dim ==1 && geom == T) // T line
            {
              return new TGauss<1, 4>;
            }
          else if(dim ==1 && geom == Q) // Q line
            {
              return new Gauss<1, 4>;
            }
          else if((dim == 2) && geom == T) //triangle
            {
              return new TGauss<2, 4>;
            }
          else if((dim == 2) && geom == Q) // square
            {
              return new Gauss<2, 4>;
            }
          else if((dim == 3) && geom == T) // tet
            {
              // Note: there's no scheme for <3,4> so we use <3,5> instead
              return new TGauss<3, 5>;
            }
          else if((dim == 3) && geom == Q) // brick
            {
              return new Gauss<3, 4>;
            }
        }
      std::string err("Cannot determine element type.\n");
      err += "Only implemented for 1 < nnode1d < 4, is that the problem?";
      err += "I got dim = " + to_string(dim) + ", nnode1d = " + to_string(nnode_1d);
      err += " and element geometry = " + to_string(geom);
      throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }


  }
}
