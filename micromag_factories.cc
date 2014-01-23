
#include "micromag_factories.h"


#include "llg_problem.h"


#include "boundary_element_handler.h"

// mesh factories
#include "../../src/generic/mesh.h"
#include "multi_mesh.h"

// solver factories
#include "../../src/generic/linear_solver.h"
#include "../../src/generic/iterative_linear_solver.h"
#include "../../src/generic/hypre_solver.h"
#include "../../src/generic/sum_of_matrices.h"

// timestepper factories
#include "../../src/generic/timesteppers.h"
#include "../../src/generic/midpoint_method.h"
#include "../../src/generic/explicit_timesteppers.h"

// Preconditioner factories
#include "../../src/generic/preconditioner.h"
#include "../../src/generic/general_purpose_preconditioners.h"
#include "../../src/generic/general_purpose_block_preconditioners.h"
#include "sum_of_matrices_preconditioner.h"


namespace oomph
{
  namespace Factories
  {

    void bem_handler_factory(BoundaryElementHandler& new_bem_handler,
                             const BemBoundaryData& bem_boundaries,
                             const unsigned& phi_index,
                             const unsigned& phi_1_index,
                             const CornerDataInput& input_corner_data,
                             int use_hierarchical_bem,
                             bool disable_corner_angles,
                             int use_numerical_integration)
    {
      // Figure out what defaults to use if any bool options are -1
      // ============================================================

      FiniteElement* bulk_fe_pt = bem_boundaries[0].second->finite_element_pt(0);

      // Use H-lib if possible (have it and surface mesh is triangular)
      if(use_hierarchical_bem == -1)
        {
#ifdef OOMPH_HAS_HLIB
          if((bulk_fe_pt->nodal_dimension() == 3)
             && (bulk_fe_pt->nnode_1d() == 2)
             && (bulk_fe_pt->nnode() == 4))
            {
              use_hierarchical_bem = true;
            }
          else
            {
              use_hierarchical_bem = false;
            }
#else
          use_hierarchical_bem = false;
#endif
        }

      // Use analytical integation if possible, numerical otherwise
      if(use_numerical_integration == -1)
        {
          if((bulk_fe_pt->nodal_dimension() == 3)
             && (bulk_fe_pt->nnode_1d() == 2)
             && (bulk_fe_pt->nnode() == 4))
            {
              use_numerical_integration = false;
            }
          else
            {
              use_numerical_integration = true;
            }
        }


      // Next assign the parameters
      // ============================================================

      // Check that we can do hierarchical bem, if so set the parameter.
      if(use_hierarchical_bem)
        {
#ifndef OOMPH_HAS_HLIB
          std::string err = "Hlib library required for hierarchical bem matrix";
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
#endif
        }
      new_bem_handler.Use_hierarchical_bem = use_hierarchical_bem;

      // Get the first finite element we can find in the a bulk mesh on
      // which we are going to construct our bem mesh. Assuimg that all
      // elements in all the meshes are the same type...
      FiniteElement* sample_fele_pt =
        bem_boundaries[0].second->finite_element_pt(0);

      // Figure out which element type we should use in the bem mesh
      // (based on the element type used in the bulk mesh) and store the
      // function needed to create them.
      new_bem_handler.Bem_element_factory_fpt = LLGFactories::
        bem_element_factory_factory(sample_fele_pt);

      // Create an integration scheme ??ds move this outside somewhere...
      new_bem_handler.integration_scheme_pt() = LLGFactories::
        variable_order_integrator_factory(sample_fele_pt);

      // Set indexes to look in for phi/phi1 variables
      new_bem_handler.set_input_index(phi_1_index);
      new_bem_handler.set_output_index(phi_index);

      // Copy in the list of boundaries to operate on
      new_bem_handler.Bem_boundaries = bem_boundaries;

      // Set debug parameters
      new_bem_handler.Debug_disable_corner_contributions = disable_corner_angles;
      new_bem_handler.Use_numerical_integration = use_numerical_integration;

      // Now build it
      new_bem_handler.build(input_corner_data);
    }


    BoundaryElementHandler* bem_handler_factory
    (const BemBoundaryData& bem_boundaries,
     const unsigned& phi_index,
     const unsigned& phi_1_index,
     const CornerDataInput& input_corner_data,
     int use_hierarchical_bem,
     bool disable_corner_angles,
     int use_numerical_integration)
    {
      // Create with new, fill in with factory
      BoundaryElementHandler* bem_handler_pt = new BoundaryElementHandler;
      bem_handler_factory(*bem_handler_pt,
                          bem_boundaries, phi_index, phi_1_index,
                          input_corner_data, use_hierarchical_bem, disable_corner_angles,
                          use_numerical_integration);
      return bem_handler_pt;
    }


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
    TimeStepper* time_stepper_factory
    (const std::string& ts_name, const std::string& mp_pred_name)
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
      else if(ts_name == "midpoint")
        {
          MidpointMethod* mp_pt = new MidpointMethod(adaptive_flag);
          ExplicitTimeStepper* pred_pt = explicit_time_stepper_factory(mp_pred_name);
          mp_pt->set_predictor_pt(pred_pt);
          return mp_pt;
        }
      else if(ts_name == "midpoint-bdf")
        {
          MidpointMethodByBDF* mp_pt = new MidpointMethodByBDF(adaptive_flag);
          ExplicitTimeStepper* pred_pt = explicit_time_stepper_factory(mp_pred_name);
          mp_pt->set_predictor_pt(pred_pt);
          return mp_pt;
        }
      else if(ts_name == "steady")
        {
          // 2 steps so that we have enough space to do reasonable time
          // derivative estimates in e.g. energy derivatives.
          return new Steady<3>;
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
                                             nnode1d);

          // Shift it
          shift_mesh(mesh_details[j].xshift,
                     mesh_details[j].yshift,
                     mesh_details[j].zshift,
                     mesh_pt);

          // Scale it. We do this after shifting so that we are scaling the
          // entire multi-mesh, including gaps between meshes. This is much
          // more intuitive.
          scale_mesh(scaling, mesh_pt);

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
                               time_stepper_pt, scaling, nnode1d);
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
                               time_stepper_pt, scaling, nnode1d);
    }



    LinearSolver* linear_solver_factory(const std::string& _solver_name)
    {
      const std::string solver_name = to_lower(_solver_name);

      LinearSolver* solver_pt;

      if(solver_name == "superlu")
        { solver_pt = new SuperLUSolver; }
      else if(solver_name == "gmres")
        {
          IterativeLinearSolver* its_pt = new GMRES<CRDoubleMatrix>;
          its_pt->max_iter() = 200;
          solver_pt = its_pt;
        }
      else if(solver_name == "cg")
        {
          IterativeLinearSolver* its_pt = new CG<CRDoubleMatrix>;
          its_pt->max_iter() = 200;
          solver_pt = its_pt;
        }
      else if(solver_name == "fdlu")
        { solver_pt = new FD_LU; }
      else if(solver_name == "som-gmres")
        {
          IterativeLinearSolver* its_pt = new GMRES<SumOfMatrices>;
          its_pt->max_iter() = 200;
          solver_pt = its_pt;
        }
      else
        {
          std::string err("Unrecognised solver name ");
          err += solver_name;
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      return solver_pt;
    }


    Preconditioner* preconditioner_factory(const std::string &_prec_name)
    {
      const std::string prec_name = to_lower(_prec_name);
      Preconditioner* prec_pt = 0;

      // AMG with optimal poisson paramters
      // ============================================================
      if(prec_name == "poisson-amg")
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

      else if(prec_name == "ilu0")
        { prec_pt = new ILUZeroPreconditioner<CRDoubleMatrix>; }

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


      // if it starts with blockllg then call that factory
      else if(split_string(prec_name, '-')[0] == "blockllg")
        {
          prec_pt = block_llg_factory(prec_name);
        }

      else
        {
          std::string err("Unrecognised preconditioner name ");
          err += prec_name;
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      return prec_pt;
    }



    Preconditioner* block_llg_factory(const std::string &_prec_name)
    {
      // Parse the parameter string
      Vector<std::string> parameters = split_string(to_lower(_prec_name), '-');

#ifdef PARANOID
      if(parameters[0] != "blockllg")
        {
          std::string error_msg = _prec_name + " is not a block llg preconditioner";
          error_msg += "(should begin with blockllg-).";
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      if(parameters.size() != 5)
        {
          std::string error_msg = "Not enough parameters in llg block string.";
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      std::cout << "block llg parameters = " << parameters << std::endl;

      // Pick the basic block structure
      GeneralPurposeBlockPreconditioner<CRDoubleMatrix>* bp_pt = 0;
      if(parameters[1] == "blockexact")
        {
          bp_pt = new ExactBlockPreconditioner<CRDoubleMatrix>;
        }
      else if(parameters[1] == "uppertriangular")
        {
          BlockTriangularPreconditioner<CRDoubleMatrix>* tri_prec_pt
            = new BlockTriangularPreconditioner<CRDoubleMatrix>;
          tri_prec_pt->upper_triangular();
          bp_pt = tri_prec_pt;
        }
      else if(parameters[1] == "lowertriangular")
        {
          BlockTriangularPreconditioner<CRDoubleMatrix>* tri_prec_pt
            = new BlockTriangularPreconditioner<CRDoubleMatrix>;
          tri_prec_pt->lower_triangular();
          bp_pt = tri_prec_pt;
        }
      else if(parameters[1] == "blockdiagonal")
        {
          bp_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;
        }
      else
        {
          throw OomphLibError("Unrecognised block structure setting "
                              + parameters[1], OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }


      // Create + set reordering of blocks: the two dofs given (x = m_x
      // etc.) are put into the first block, the other magnetisation dof is
      // put into the second block. Poisson dofs are in their own blocks
      // (empty for now).
      Vector<unsigned> a(5);
      if(parameters[3] == "xy")
        {
          a[0] = 2; a[1] = 3; // poisson blocks
          a[2] = 0; a[3] = 0; // first m block
          a[4] = 1; // second m block
        }
      else if(parameters[3] == "xz")
        {
          a[0] = 2; a[1] = 3; // poisson blocks
          a[2] = 0; a[4] = 0; // first m block
          a[3] = 1; // second m block
        }
      else if(parameters[3] == "yz")
        {
          a[0] = 2; a[1] = 3; // poisson blocks
          a[4] = 0; a[4] = 0; // first m block
          a[2] = 1; // second m block
        }
      else
        {
          throw OomphLibError("Unrecognised block swapping setting "
                              + parameters[3], OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      bp_pt->set_dof_to_block_map(a);


      // Pick the solver to use on the individual blocks (or on the entire
      // thing if we are using "blockexact").
      if(parameters[4] == "exact")
        {
          // Do nothing--default GeneralPurposeBlockPreconditioner solver
          // is SuperLU.
        }
      // else if(parameters[4] == "amg")
      //   {
      //     bp_pt->set_subsidiary_preconditioner_function();
      //   }
      else
        {
          throw OomphLibError("Unrecognised block structure setting "
                              + parameters[4], OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }


      // Now create the upper-left block sub preconditioner.
      // ============================================================
      GeneralPurposeBlockPreconditioner<CRDoubleMatrix>* sub_bp_pt = 0;
      if(parameters[2] == "blockexact")
        {
          sub_bp_pt = new ExactBlockPreconditioner<CRDoubleMatrix>;
        }
      else if(parameters[2] == "uppertriangular")
        {
          BlockTriangularPreconditioner<CRDoubleMatrix>* tri_prec_pt
            = new BlockTriangularPreconditioner<CRDoubleMatrix>;
          tri_prec_pt->upper_triangular();
          sub_bp_pt = tri_prec_pt;
        }
      else if(parameters[2] == "lowertriangular")
        {
          BlockTriangularPreconditioner<CRDoubleMatrix>* tri_prec_pt
            = new BlockTriangularPreconditioner<CRDoubleMatrix>;
          tri_prec_pt->lower_triangular();
          sub_bp_pt = tri_prec_pt;
        }
      else if(parameters[2] == "blockdiagonal")
        {
          sub_bp_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;
        }
      else if(parameters[2] == "blockantidiagonal")
        {
          sub_bp_pt = new BlockAntiDiagonalPreconditioner<CRDoubleMatrix>;
        }
      else
        {
          std::string err = "Unrecognised sub-preconditioner block structure";
          err += " setting " + parameters[2];
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      // The blocks we want to use it on are just the ones which are in the
      // zeroth block of the master. Find out which ones these are.
      Vector<unsigned> sub_bp_mapping;
      for(unsigned j=0; j<a.size(); j++)
        {if(a[j] == 0) sub_bp_mapping.push_back(j);}

      // And set the master pointer along with this mapping
      sub_bp_pt->turn_into_subsidiary_block_preconditioner(bp_pt,
                                                           sub_bp_mapping);

      // Provide our subsidiary preconditioner to the master for use on
      // block 0, which is the first m block.
      bp_pt->set_subsidiary_preconditioner_pt(sub_bp_pt, 0);

      return bp_pt;
    }



    Vector<unsigned> dof_to_block_factory(const std::string& _name)
    {
      const std::string name = to_lower(_name);

      // Make an element to look up indicies from
      TMicromagElement<2,2> dummy_ele;

      const unsigned ndof = dummy_ele.ndof_types(); //??ds unsafe?
      Vector<unsigned> dof_to_block(ndof);

      if(name == "none")
        {
          // identity mapping
          for(unsigned j=0; j<ndof; j++)
            {
              dof_to_block[j] = j;
            }
        }

      // All m values in one block, others left alone.
      // [0, 1, 2, 2, 2, 3, 4]
      else if(name == "group-m")
        {
          unsigned k = 0;

          // Normal phi/phi1
          dof_to_block[dummy_ele.phi_index_micromag()] = k++;
          dof_to_block[dummy_ele.phi_1_index_micromag()] = k++;

          // m all in one block
          for(unsigned j=0; j<3; j++)
            {
              int index = dummy_ele.m_index_micromag(j);
              dof_to_block[index] = k;
            }
          k++;

          // boundary phi/phi1 ??ds assume they are at the end...
          dof_to_block[5] = k++;
          dof_to_block[6] = k++;
        }
      else if(name == "group-m-phi-phi-boundary")
        {
          unsigned k = 0;

          // All phi into one block
          dof_to_block[dummy_ele.phi_index_micromag()] = k;
          dof_to_block[5] = k;
          k++;

          // Simiarly for phi1
          dof_to_block[dummy_ele.phi_1_index_micromag()] = k;
          dof_to_block[6] = k;
          k++;

          // m all in one block
          for(unsigned j=0; j<3; j++)
            {
              int index = dummy_ele.m_index_micromag(j);
              dof_to_block[index] = k;
            }
          k++;
        }
      else
        {
          std::string err = "Unrecognised blocking name";
          err += name;
          throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }

      return dof_to_block;
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


    // LLGProblem* llg_problem_factory(MMArgs* args_pt)
    // {
    //   LLGProblem* problem_pt = new LLGProblem;

    //   my_problem_factory_helper( *problem_pt)

    // }

  }
}
