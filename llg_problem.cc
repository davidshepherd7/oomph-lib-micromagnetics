

#include "llg_problem.h"
#include "micromag_factories.h"
#include "micromag_types.h"

// Meshes for mesh factory
#include "../../src/meshes/simple_rectangular_quadmesh.h"
#include "../../src/meshes/simple_rectangular_tri_mesh.h"
#include "../../src/meshes/simple_cubic_tet_mesh.h"
#include "../../src/meshes/simple_cubic_mesh.h"
#include "../../src/meshes/tetgen_mesh.h"
#include "../../src/meshes/triangle_mesh.h"

#include "./multi_mesh.h"
#include "./single_element_mesh.h"

#include "micromagnetics_element.h"
#include "magnetostatic_field_flux_element.h"

namespace oomph
{

  /// Function that does the real work of the constructors.
  void LLGProblem::build(Vector<Mesh*>& bulk_mesh_pts)
  {
#ifdef PARANOID
    if((applied_field_fct_pt() == 0)
       || (Residual_calculator_pt == 0))
      {
        std::ostringstream error_msg;
        error_msg
          << "Must the following pointers to non-null values before calling build():\n"
          << "applied_field_fct_pt() (= " << applied_field_fct_pt() << ")\n"
          << "Residual_calculator_pt (= " << Residual_calculator_pt << ")"
          << std::endl;

        throw OomphLibError(error_msg.str(), OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

    if(Analytic_ms_fct_pt != 0 && !Disable_ms)
      {
        std::string err = "Other ms must be disabled to use analytical ms.";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif

    // Call the underlying build to deal with adding meshes and time stepper
    MyProblem::build(bulk_mesh_pts);

    // Finish off element build, at this point we should have only micromag
    // elements in the meshes (so we can loop over all meshes) but we don't
    // have a global mesh yet.
    for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
      {
        for(unsigned i=0; i<mesh_pt(msh)->nelement(); i++)
          {
            MicromagEquations* elem_pt =checked_dynamic_cast<MicromagEquations*>
              (mesh_pt(msh)->element_pt(i));

            // Set whether the Jacobian should be finite differenced
            elem_pt->Use_fd_jacobian = Use_fd_jacobian;

            // Set values for magnetic parameters
            elem_pt->magnetic_parameters_pt() = mag_parameters_pt();

            // Set pointer for an applied field
            elem_pt->applied_field_pt() = applied_field_fct_pt();

            // Set the residual calculation function
            elem_pt->Residual_calculator_pt = Residual_calculator_pt;

            if(Analytic_ms_fct_pt != 0)
              {
                AnalyticalMagnetostatics* ams_pt = new AnalyticalMagnetostatics;
                ams_pt->Magnetostatic_field_fct_pt = Analytic_ms_fct_pt;
                elem_pt->Ms_calc_pt = ams_pt;
              }
          }
      }

    // For debugging we might want to pin m values on the boundary
    if(Pin_boundary_m)
      {
        // Loop over all meshes in problem
        for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
          {
            Mesh* mesh_pt = this->mesh_pt(msh);
            for(unsigned b=0, nb=mesh_pt->nboundary(); b<nb; b++)
              {
                for(unsigned nd=0, nnd=mesh_pt->nboundary_node(b); nd<nnd; nd++)
                  {
                    Node* nd_pt = mesh_pt->boundary_node_pt(b, nd);
                    for(unsigned j=0; j<3; j++)
                      {
                        nd_pt->pin(m_index(j));
                      }
                  }
              }

          }
      }

    // Set up bem stuff if we are doing it
    if(implicit_ms_flag())
      {
        // Figure out how to build the flux meshes that we're going to need
        // for neumann boundaries.
        Flux_mesh_factory_pt = LLGFactories::mm_flux_mesh_factory_factory
          (bulk_mesh_pts[0]->finite_element_pt(0));

        // Loop over all meshes in problem
        for(unsigned msh=0, nmsh=bulk_mesh_pts.size(); msh<nmsh; msh++)
          {
            Vector<unsigned> boundaries;
            for(unsigned b=0; b<bulk_mesh_pts[msh]->nboundary(); b++)
              {boundaries.push_back(b);}

            // Set up neumann condition on phi_1 boundary values (using flux mesh)
            Flux_mesh_pt = flux_mesh_factory(bulk_mesh_pts[msh], boundaries);

            // Add to global mesh
            this->add_sub_mesh(Flux_mesh_pt);


            // Pin a phi_1 value which isn't involved in the boundary element
            // method (we have to pin something to avoid a singular Jacobian,
            // can't be a boundary node or things will go wrong with BEM).
#ifdef OOMPH_HAS_MPI
            // In parallel we need to make sure that only one node is
            // pinned in total
            std::string err = "Not implemented!";
            throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                OOMPH_CURRENT_FUNCTION);

            // Check that processor id is 0, if so then pin as for serial,
            // otherwise do nothing.

            //??ds Could be problems when nodes duplicated? Not sure how
            //all that works yet?
#else
            Node* pinned_phi_1_node_pt = bulk_mesh_pts[msh]->get_some_non_boundary_node();
            pinned_phi_1_node_pt->pin(phi_1_index());
            pinned_phi_1_node_pt->set_value(phi_1_index(), 0.0);
#endif

          }

      }
    // Otherwise pin all phi and phi_1 dofs to zero
    else
      {
        oomph_info << "Pinning phi values in main problem's meshes." << std::endl;

        // Loop over all meshes in problem
        for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
          {
            for(unsigned nd=0, nnode=mesh_pt(msh)->nnode(); nd<nnode; nd++)
              {
                Node* nd_pt = mesh_pt(msh)->node_pt(nd);

                nd_pt->pin(phi_index());
                nd_pt->pin(phi_1_index());
                nd_pt->set_value(phi_index(),0.0);
                nd_pt->set_value(phi_1_index(),0.0);
              }
          }
      }


    // If we are using llg residual with midpoint method then we need to
    // swap residuals over in the explicit predictor time steps. Put in
    // the class to do this here.
    if(Residual_calculator_pt->use_gilbert_form()
       && dynamic_cast<MidpointMethodBase*>(time_stepper_pt()) != 0)
      {
        MidpointMethodBase* mp_pt = checked_dynamic_cast<MidpointMethodBase*>
          (time_stepper_pt());

        // Create and set up our residual swapping timestepper.
        ResidualSwappingExplicitTimestepper* rsts_pt
          = new ResidualSwappingExplicitTimestepper;
        rsts_pt->underlying_time_stepper_pt = mp_pt->predictor_pt();
        rsts_pt->residual_pt = checked_dynamic_cast<LLGResidualCalculator*>
          (Residual_calculator_pt);

        mp_pt->set_predictor_pt(rsts_pt);
      }



    // Build the global mesh
    this->build_global_mesh();

    // Number the equations
    this->assign_eqn_numbers();

    // Write out some stuff
    mag_parameters_pt()->output(*oomph_info.stream_pt());
    oomph_info << "LLG Number of equations: " << ndof() << std::endl;
    oomph_info << "Number of sub meshes: " << this->nsub_mesh() << std::endl;

  }

  void LLGProblem::build_decoupled_ms(Vector<Mesh*>& llg_mesh_pts,
                                      Vector<Mesh*>& phi_mesh_pts,
                                      Vector<Mesh*>& phi_1_mesh_pts)
  {

    // Throughout this function we need to be careful not to use each
    // problem's global mesh pointer because they are built only during
    // this function.

    // Set up phi_1 problem
    // ============================================================

    //??ds move outside?
    Phi_1_problem_pt = new GenericPoissonProblem;
    phi_1_problem_pt()->set_flux_mesh_factory(Phi_1_flux_mesh_factory_fct_pt);
    phi_1_problem_pt()->newton_solver_tolerance() = newton_solver_tolerance();

    // phi_1 b.c.s are all Neumann but with flux determined elsewhere (by m).
    for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
      {
        Mesh* mesh_pt = phi_1_mesh_pts[msh];
        for(unsigned b=0, nb=mesh_pt->nboundary(); b < nb; b++)
          {
            phi_1_problem_pt()->add_neumann_boundary(mesh_pt, b, 0);
          }
      }

    // Finish off the problem
    phi_1_problem_pt()->build(phi_1_mesh_pts);


    // Check that all three types of mesh match up
#ifdef PARANOID
    if((llg_mesh_pts.size() != phi_mesh_pts.size())
       || (llg_mesh_pts.size() != phi_1_mesh_pts.size()))
      {
        std::string err = "All mesh lists must be the same size!";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif


    // Phi problem
    // ============================================================

    Phi_problem_pt = new GenericPoissonProblem;
    phi_problem_pt()->newton_solver_tolerance() = newton_solver_tolerance();

    // phi b.c.s are all Dirichlet, determined by a vector of BEM
    // values. Create these vectors and link them into the phi problem.
    for(unsigned msh=0, nmsh=phi_mesh_pts.size(); msh<nmsh; msh++)
      {
        Mesh* mesh_pt = phi_mesh_pts[msh];
        for(unsigned b=0, nb=mesh_pt->nboundary(); b < nb; b++)
          {
            // Create the vector
            LinearAlgebraDistribution* dist_pt =
              new LinearAlgebraDistribution(MPI_Helpers::communicator_pt(),
                                            mesh_pt->nboundary_node(b), false);
            DoubleVector* db_vec_pt = new DoubleVector(dist_pt);

            // Plug it into the phi problem and our list of bem controlled
            // values.
            Phi_boundary_values_pts.push_back(db_vec_pt);
            phi_problem_pt()->add_dirichlet_boundary_by_vector(mesh_pt, b,
                                                               db_vec_pt);
          }

      }

    phi_problem_pt()->build(phi_mesh_pts);


    // Coupling between problems
    // ============================================================

    // Assign bulk elements pointers to each other (and check that it is
    // safe).
    for(unsigned msh=0; msh<llg_mesh_pts.size(); msh++)
      {
        Mesh* phi_mesh_pt = phi_mesh_pts[msh];
        Mesh* llg_mesh_pt = llg_mesh_pts[msh];
        Mesh* phi_1_mesh_pt = phi_1_mesh_pts[msh];

#ifdef PARANOID
        // Things will go wrong if the nodes of all meshes are not in
        // the same place:
        if((phi_1_mesh_pt->nnode() != phi_mesh_pt->nnode())
           || (phi_1_mesh_pt->nnode() != llg_mesh_pt->nnode()))
          {
            std::ostringstream error_msg;
            error_msg << "Mesh nnodes must be the same.";
            throw OomphLibError(error_msg.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }

        if((phi_1_mesh_pt->nelement() != phi_mesh_pt->nelement())
           || (phi_1_mesh_pt->nelement() != llg_mesh_pt->nelement()))
          {
            std::ostringstream error_msg;
            error_msg << "Mesh nelements must be the same.";
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
                    Vector<double> phi1x(dim), phix(dim), llgx(dim);
                    phi_1_mesh_pt->node_pt(nd)->position(phi1x);
                    phi_mesh_pt->node_pt(nd)->position(phix);
                    llg_mesh_pt->node_pt(nd)->position(llgx);

                    std::ostringstream error_msg;
                    error_msg << "Mesh nodes must be in the same places."
                              << " Problem in node " << nd
                              << "\nphi1 position " << phi1x
                              << "\nphi position " << phix
                              << "\nllg position " << llgx;
                    throw OomphLibError(error_msg.str(),
                                        OOMPH_CURRENT_FUNCTION,
                                        OOMPH_EXCEPTION_LOCATION);
                  }
              }
          }
#endif

        // Assign the various elements pointers to each other
        for(unsigned e=0, ne=phi_1_mesh_pt->nelement(); e < ne; e++)
          {

            // Get the element pointers
            MagnetostaticFieldEquations* phi_1_ele_pt =
              checked_dynamic_cast<MagnetostaticFieldEquations*>
              (phi_1_mesh_pt->element_pt(e));
            MagnetostaticFieldEquations* phi_ele_pt =
              checked_dynamic_cast<MagnetostaticFieldEquations*>
              (phi_mesh_pt->element_pt(e));
            MicromagEquations* m_ele_pt =
              checked_dynamic_cast<MicromagEquations*>
              (llg_mesh_pt->element_pt(e));

            phi_1_ele_pt->set_micromag_element_pt(m_ele_pt);
            phi_ele_pt->set_micromag_element_pt(m_ele_pt);
            checked_dynamic_cast<SemiImplicitMagnetostaticsCalculator*>
              (m_ele_pt->Ms_calc_pt)
              ->magnetostatic_field_element_pt() = phi_ele_pt;
          }
      }

  }


  /// \short Error for adaptive timestepper (rms of nodal error determined by
  /// comparison with explicit timestepper result).
  double LLGProblem::global_temporal_error_norm()
  {
    double global_error = 0.0;

    // Loop over the nodes
    for(unsigned i=0, ni=mesh_pt()->nnode(); i<ni; i++)
      {
        Node* nd_pt = mesh_pt()->node_pt(i);
        for(unsigned j=0; j<3; j++)
          {
            // Get timestepper's error estimate for this direction of m
            // at this point.
            double error = nd_pt->time_stepper_pt()->
              temporal_error_in_value(nd_pt, m_index(j));

            //Add the square of the individual error to the global error
            global_error += error*error;
          }
      }

    // Divide by the number of data points
    global_error /= 3*double(mesh_pt()->nnode());

    return std::sqrt(global_error);
  }


  /// \short Loop over all nodes in bulk mesh and get magnetisations
  Vector<Vector<double> > LLGProblem::get_nodal_magnetisations(unsigned i_time) const
  {
    // Get the space needed
    unsigned nnode = mesh_pt()->nnode();
    Vector< Vector<double> > m_list(nnode, Vector<double>(3, 0.0));

    for(unsigned nd=0, nnd=mesh_pt()->nnode(); nd<nnd; nd++)
      {
        for(unsigned j=0; j<3; j++)
          {
            m_list[nd][j] =
              mesh_pt()->node_pt(nd)->value(i_time, m_index(j));
          }
      }

    return m_list;
  }

  /// \short Solve for the magnetostatic field.
  void LLGProblem::magnetostatics_solve()
  {
#ifdef PARANOID
    if(!Decoupled_ms)
      {
        std::string err = "Requested a decoupled magnetostatics solve";
        err += "  but problem is not decoupled!";
        throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif

    oomph_info << std::endl
               << "Decoupled BEM solve" << std::endl
               << "--------------------------" <<std::endl;

    // solve for phi1
    oomph_info << "solving phi1" << std::endl;
    phi_1_problem_pt()->newton_solve();

    // update boundary values of phi
    oomph_info << "solving BEM" << std::endl;
    double t_start = TimingHelpers::timer();
    Bem_handler_pt->get_bem_values(Phi_boundary_values_pts);
    double t_end = TimingHelpers::timer();
    oomph_info << "BEM time taken: " << t_end - t_start << std::endl;

    // solve for phi
    oomph_info << "solving phi" << std::endl;
    phi_problem_pt()->newton_solve();

    oomph_info << "mean field is " << average_magnetostatic_field() << std::endl;

    Decoupled_ms_has_been_calculated = true;
  }


  /// Linearly extrapolate phi
  void LLGProblem::extrapolate_phi(const double& new_dt, const double& prev_dt)
  {
    // Don't change phi_1 values because we don't need them except for
    // calculating phi.

    //??ds
    double dtnp1 = time_stepper_pt()->time_pt()->dt();
    double dtn = time_stepper_pt()->time_pt()->dt(1);

    const unsigned phi_index = phi_problem_pt()->poisson_dof_number();

    // Loop over all meshes in problem
    for(unsigned msh=0, nmsh=phi_problem_pt()->nsub_mesh(); msh<nmsh; msh++)
      {
        Mesh* mesh_pt = phi_problem_pt()->mesh_pt(msh);
        for(unsigned nd=0, nnd=mesh_pt->nnode(); nd<nnd; nd++)
          {
            Node* nd_pt = mesh_pt->node_pt(nd);
            double phi_nm1 = nd_pt->value(2, phi_index);
            double phi_n = nd_pt->value(1, phi_index);

            double phi_np1 = ((dtnp1 + dtn)/dtn)*phi_n - (dtnp1/dtn)*phi_nm1;

            nd_pt->set_value(0, phi_index, phi_np1);
          }
      }

  }


  /// \short Abs of mean difference of actual m and m given by a function
  /// at the middle of each element.
  double LLGProblem::compare_m_with_function(const InitialM::InitialMFctPt fct_pt) const
  {
    double diff = 0.0;

    // Compare at middle of element
    Vector<double> s(3,0.0);
    for(unsigned j=0; j<dim(); j++) s[j] = 0.5;

    // Sum the difference over all bulk elements in problem
    unsigned bulk_ele_count = 0;

    // Loop over all meshes in problem
    for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
      {
        // Skip non-bulk meshes
        if((mesh_pt(msh)->nnode() == 0)
           || (mesh_pt(msh)->node_pt(0)->ndim() != Dim)) continue;

        for(unsigned e=0, ne=mesh_pt(msh)->nelement(); e < ne; e++)
          {
            MicromagEquations* ele_pt = checked_dynamic_cast<MicromagEquations*>
              (mesh_pt(msh)->element_pt(e));

            Vector<double> numerical_m(3,0.0);
            ele_pt->interpolated_m_micromag(s, numerical_m);

            Vector<double> x(dim(),0.0);
            ele_pt->interpolated_x(s,x);
            Vector<double> exact_m = fct_pt(time(), x);

            for(unsigned j=0; j<3; j++)
              {
                diff += std::abs(numerical_m[j] - exact_m[j]);
              }

            bulk_ele_count++;
          }
      }

    // Divide to get the mean
    diff /= (3.0 * double(bulk_ele_count));

    return diff;
  }


  double LLGProblem::mean_norm_m_error() const
  {
    double temp_error = 0;

    // Loop over all nodes in problem
    for(unsigned nd=0, nnd=mesh_pt()->nnode(); nd<nnd; nd++)
      {
        Node* nd_pt = mesh_pt()->node_pt(nd);
        Vector<double> m_values(3, 0.0);
        for(unsigned j=0; j<3; j++) m_values[j] = nd_pt->value(m_index(j));

        double m_2norm = VectorOps::two_norm(m_values);
        temp_error += std::abs(m_2norm - 1);
      }

    return temp_error/double(mesh_pt()->nnode());
  }

  /// \short Get a vector of the mean of the nodal magnetisations
  Vector<double> LLGProblem::mean_magnetisation() const
  {
    Vector<Vector<double> > ms = get_nodal_magnetisations();
    Vector<double> mean_m(3, 0.0);

    unsigned n_ms = ms.size();
    for(unsigned i=0; i<n_ms; i++)
      {
        mean_m[0] += ms[i][0];
        mean_m[1] += ms[i][1];
        mean_m[2] += ms[i][2];
      }

    mean_m[0] /= n_ms;
    mean_m[1] /= n_ms;
    mean_m[2] /= n_ms;

    return mean_m;
  }


  void LLGProblem::norm_m_error(double &m_error_avg, double &m_error_stddev) const
  {
    // Get mean from other function
    m_error_avg = mean_norm_m_error();

    // Calculate std deviation
    // ============================================================
    double temp_stddev = 0.0;
    for(unsigned nd=0, nnd=mesh_pt()->nnode(); nd<nnd; nd++)
      {
        Node* nd_pt = mesh_pt()->node_pt(nd);
        Vector<double> m_values(3,0.0);
        for(unsigned j=0; j<3; j++) m_values[j] = nd_pt->value(m_index(j));

        double m_2norm = VectorOps::two_norm(m_values);
        temp_stddev += std::pow( m_2norm - 1 - m_error_avg, 2);
      }

    temp_stddev /= double(mesh_pt()->nnode());
    m_error_stddev = std::sqrt(temp_stddev);
  }


  /// \short Functions for constructing things needed for LLG problems
  namespace LLGFactories
  {
    /// \short Make a mesh as specified by an input argument. Refined
    /// according to the given refinement level (in some way appropriate
    /// for that mesh type). Assumption: this will be passed into a
    /// problem, which will delete the pointer when it's done.
    Mesh* mesh_factory(const std::string& _mesh_name,
                       int refinement_level,
                       TimeStepper* time_stepper_pt,
                       double scaling_factor,
                       unsigned nnode1d)
    {
      // Ignore case in mesh names
      const std::string mesh_name = to_lower(_mesh_name);

      // Refinement always roughly the same for structured meshes
      unsigned nx = 5 * std::pow(2, refinement_level-1);

      // Make the mesh and store a pointer to it
      Mesh* mesh_pt = 0;
      if(mesh_name == "sq_square" && nnode1d == 2)
        {
          double lx = 1.0;
          mesh_pt = new SimpleRectangularQuadMesh<QMicromagElement<2,2> >
            (nx, nx, lx, lx, time_stepper_pt);
        }
      else if(mesh_name == "st_square" && nnode1d == 2)
        {
          double lx = 1.0;
          mesh_pt = new SimpleRectangularTriMesh<TMicromagElement<2,2> >
            (nx, nx, lx, lx, time_stepper_pt);

          mesh_pt->setup_boundary_element_info();

          // Turn off triangle refinement dump stuff (breaks Micromag
          // elements).
          checked_dynamic_cast<TriangleMeshBase*>(mesh_pt)->
            disable_triangulateio_restart();
         }
      else if(mesh_name == "single-element" && nnode1d == 2)
        {
          mesh_pt = new SingleElementMesh<QMicromagElement<2,2> >(time_stepper_pt);
        }
      else if(mesh_name == "ut_square" && nnode1d == 2)
        {
          mesh_pt = new TriangleMesh<TMicromagElement<2, 2> >
            ("./meshes/square." + to_string(refinement_level) + ".node",
             "./meshes/square." + to_string(refinement_level) + ".ele",
             "./meshes/square." + to_string(refinement_level) + ".poly",
             time_stepper_pt);

          // Turn off triangle refinement dump stuff (breaks Micromag
          // elements).
          checked_dynamic_cast<TriangleMeshBase*>(mesh_pt)->
            disable_triangulateio_restart();
        }
      else if(mesh_name == "st_cubeoid" && nnode1d == 2)
        {
          // nmag cubeoid
          double lx = 1, ly = lx, lz = 3*lx;
          unsigned ny = nx, nz = std::ceil(lz/lx) * nx;
          mesh_pt = new SimpleCubicTetMesh<TMicromagElement<3, 2> >
            (nx, ny, nz, lx, ly, lz, time_stepper_pt);

          mesh_pt->setup_boundary_element_info();
        }
      else if(mesh_name == "sqt_cubeoid" && nnode1d == 2)
        {
          Mesh* qmesh_pt = mesh_factory("sq_cubeoid", refinement_level,
                                        time_stepper_pt,
                                        1, nnode1d);
          //??ds memory leak, fix? Can't delete this mesh or nodes will
          //go...

          TetMeshBase* tmesh_pt = new TetMeshBase;
          ElementFactoryFctPt factory_fpt =
            MeshCreationHelpers::new_element<TMicromagElement<3, 2> >;

          MeshCreationHelpers::brick2tet(*qmesh_pt, factory_fpt, *tmesh_pt);

          mesh_pt = tmesh_pt;
        }
      else if(mesh_name == "ut_cubeoid" && nnode1d == 2)
        {
          mesh_pt = new TetgenMesh<TMicromagElement<3, 2> >
            ("./meshes/cubeoid." + to_string(refinement_level) + ".node",
             "./meshes/cubeoid." + to_string(refinement_level) + ".ele",
             "./meshes/cubeoid." + to_string(refinement_level) + ".face",
             time_stepper_pt);
        }
      else if(mesh_name == "ut_mumag4" && nnode1d == 2)
        {
          mesh_pt = new TetgenMesh<TMicromagElement<3, 2> >
            ("./meshes/mumag4." + to_string(refinement_level) + ".node",
             "./meshes/mumag4." + to_string(refinement_level) + ".ele",
             "./meshes/mumag4." + to_string(refinement_level) + ".face",
             time_stepper_pt);
        }
      else if(mesh_name == "st_mumag4" && nnode1d == 2)
        {
          mesh_pt = new SimpleCubicTetMesh<TMicromagElement<3, 2> >
            (5*nx, std::ceil(1.25*nx), 1, 500, 125, 3, time_stepper_pt);

          mesh_pt->setup_boundary_element_info();
        }
      else if(mesh_name == "sq_mumag4" && nnode1d == 2)
        {
          unsigned this_nx = refinement_level;

          mesh_pt = new SimpleCubicMesh<QMicromagElement<3, 2> >
            (5*this_nx, std::ceil(1.25*this_nx), 2, 500, 125, 3, time_stepper_pt);

          mesh_pt->setup_boundary_element_info();
        }
      else if(mesh_name == "sqt_mumag4" && nnode1d == 2)
        {
          Mesh* qmesh_pt = mesh_factory("sq_mumag4", refinement_level,
                                        time_stepper_pt,
                                        1, nnode1d);
          //??ds memory leak, fix? Can't delete this mesh or nodes will
          //go...

          // Convert to tet mesh
          TetMeshBase* tmesh_pt = new TetMeshBase;
          ElementFactoryFctPt factory_fpt =
            MeshCreationHelpers::new_element<TMicromagElement<3, 2> >;
          MeshCreationHelpers::brick2tet(*qmesh_pt, factory_fpt, *tmesh_pt);

          mesh_pt = tmesh_pt;
        }
      else if(mesh_name == "sq_cubeoid" && nnode1d == 2)
        {
          double lx = 1, ly = lx, lz = 3*lx;
          mesh_pt = new SimpleCubicMesh<QMicromagElement<3, 2> >
            (nx, nx, int(lz/lx)*nx, lx, ly, lz, time_stepper_pt);
        }
      else if(mesh_name == "ut_sphere" && nnode1d == 2)
        {
          mesh_pt = new TetgenMesh<TMicromagElement<3, 2> >
            ("./meshes/sphere." + to_string(refinement_level) + ".node",
             "./meshes/sphere." + to_string(refinement_level) + ".ele",
             "./meshes/sphere." + to_string(refinement_level) + ".face",
             time_stepper_pt);
        }
      else
        {
          throw OomphLibError("Unrecognised mesh name " + mesh_name,
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      // Scale the mesh as requested
      scale_mesh(scaling_factor, mesh_pt);

      // This should go inside an element factory but our meshes don't
      // allow that :(
      for(unsigned ele=0, nele=mesh_pt->nelement(); ele<nele; ele++)
        {
          MicromagEquations* ele_pt = checked_dynamic_cast<MicromagEquations*>
            (mesh_pt->element_pt(ele));
          ele_pt->Ms_calc_pt = new ImplicitMagnetostaticsCalculator;
        }

      // Done: pass out the mesh pointer
      return mesh_pt;
    }

    /// Pick and create a residual calculator to use
    LLGResidualCalculator* residual_calculator_factory(const std::string& residual)
    {
      if(residual == "llg")
        return new LLGResidualCalculator(true);
      else if(residual == "ll")
        return new LLGResidualCalculator(false);
      else
        throw OomphLibError("Unrecognised residual "+residual,
                            OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
    }


    /// \short Create a variable order quadrature object based on the
    /// dimension and shape of the element. Only works for
    Integral* variable_order_integrator_factory(const FiniteElement* const el_pt)
    {
      if((el_pt->nodal_dimension() == 2) && (el_pt->nvertex_node() == 3))
        {
          return new TVariableOrderGaussLegendre<1>;
        }
      else if((el_pt->nodal_dimension() == 2) && (el_pt->nvertex_node() == 4))
        {
          return new QVariableOrderGaussLegendre<1>;
        }
      else if((el_pt->nodal_dimension() == 3) && (el_pt->nvertex_node() == 4))
        {
          return new TVariableOrderGaussLegendre<2>;
        }
      else if((el_pt->nodal_dimension() == 3) && (el_pt->nvertex_node() == 8))
        {
          return new QVariableOrderGaussLegendre<2>;
        }
      else
        {
          std::string err("Cannot determine element type.\n");
          err += "Maybe it is a higher order element (NNODE_1D > 2)?\n";
          err += "Variable order quadratures are not supported for this case.";
          throw OomphLibError(err, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }



    /// \short Return a function which will create the appropriate BEM face
    /// element for the bulk element pointer given (should work for a
    /// pointer to any bulk element type i.e., field or llg).
    BEMElementFactoryFctPt bem_element_factory_factory
    (const FiniteElement* bulk_ele_pt)
    {
      if(dynamic_cast<const TElement<2, 2>*>(bulk_ele_pt) != 0)
        {
          return &bem_element_factory<TMicromagBEMElement<2,2> >;
        }
      else if(dynamic_cast<const TElement<3, 2>*>(bulk_ele_pt) != 0)
        {
          return &bem_element_factory<TMicromagBEMElement<3,2> >;
        }

      else if(dynamic_cast<const QElement<2,2>*>(bulk_ele_pt) != 0)
        {
          return &bem_element_factory<QMicromagBEMElement<2,2> >;
        }
      else if(dynamic_cast<const QElement<3,2>*>(bulk_ele_pt) != 0)
        {
          return &bem_element_factory<QMicromagBEMElement<3,2> >;
        }

      else
        {
          throw OomphLibError("Unrecognised element type",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }

    /// \short Return a factory function which will create the appropriate
    /// "flux mesh" for the bulk element pointer given.
    FluxMeshFactoryFctPt
    mm_flux_mesh_factory_factory(const FiniteElement* bulk_ele_pt)
    {
      if(dynamic_cast<const TMicromagElement<2, 2>*>(bulk_ele_pt) != 0)
        {
          return Factories::surface_mesh_factory
            <MicromagFluxElement<TMicromagElement<2, 2> > >;
        }
      else if(dynamic_cast<const TMicromagElement<3, 2>*>(bulk_ele_pt) != 0)
        {
          return Factories::surface_mesh_factory
            <MicromagFluxElement<TMicromagElement<3, 2> > >;
        }

      else if(dynamic_cast<const QMicromagElement<2,2>*>(bulk_ele_pt) != 0)
        {
          return Factories::surface_mesh_factory
            <MicromagFluxElement<QMicromagElement<2,2> > >;
        }
      else if(dynamic_cast<const QMicromagElement<3,2>*>(bulk_ele_pt) != 0)
        {
          return Factories::surface_mesh_factory
            <MicromagFluxElement<QMicromagElement<3,2> > >;
        }

      else
        {
          throw OomphLibError("Unrecognised element type",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

    }

  }

  /// \short A namespace full of functions that take some "dynamic"
  /// (i.e. can be calculated at runtime) input and create an instance
  /// of the appropriate object, using "new" (Factory Method
  /// design pattern).
  ///
  /// Typically these objects are passed straight into other classes and
  /// will be deleted by the destructor of that class. If not it is your
  /// responsibility to make sure the objects are deleted.
  namespace SemiImplicitFactories
  {

    /// \short Make a mesh of Micromag elements as specified by an
    /// input argument. Refined according to the given refinement level (in
    /// some way appropriate for that mesh type).
    Mesh* llg_mesh_factory(const std::string& _mesh_name,
                           int refinement_level,
                           TimeStepper* time_stepper_pt,
                           double scaling_factor,
                           unsigned nnode1d)
    {
      // Construct mesh full of normal llg elements
      Mesh* mesh_pt = LLGFactories::mesh_factory(_mesh_name, refinement_level,
                                                 time_stepper_pt,
                                                 scaling_factor,
                                                 nnode1d);

      // Change to semi implicit. This should go inside an element factory
      // but our meshes don't allow that :(
      for(unsigned ele=0, nele=mesh_pt->nelement(); ele<nele; ele++)
        {
          MicromagEquations* ele_pt = checked_dynamic_cast<MicromagEquations*>
            (mesh_pt->element_pt(ele));
          delete ele_pt->Ms_calc_pt; ele_pt->Ms_calc_pt = 0;
          ele_pt->Ms_calc_pt = new SemiImplicitMagnetostaticsCalculator;
        }

      // Done: pass out the mesh pointer
      return mesh_pt;
    }


    /// \short Make a mesh of MagnetostaticField elements as specified by an
    /// input argument. Refined according to the given refinement level (in
    /// some way appropriate for that mesh type).
    Mesh* phi_mesh_factory(const std::string& _mesh_name,
                           int refinement_level,
                           TimeStepper* time_stepper_pt,
                           double scaling_factor,
                           unsigned nnode1d)
    {
      //??ds add scaling factor

      // Ignore case in mesh names
      const std::string mesh_name = to_lower(_mesh_name);

      unsigned nx = 5 * std::pow(2, refinement_level-1);


      // Make the mesh and store a pointer to it
      Mesh* mesh_pt;
      if(mesh_name == "sq_square" && nnode1d == 2)
        {
          double lx = 1.0;
          mesh_pt = new SimpleRectangularQuadMesh<QMagnetostaticFieldElement<2,2> >
            (nx, nx, lx, lx, time_stepper_pt);
        }
      else if(mesh_name == "st_square" && nnode1d == 2)
        {
          double lx = 1.0;
          mesh_pt = new SimpleRectangularTriMesh<TMagnetostaticFieldElement<2,2> >
            (nx, nx, lx, lx, time_stepper_pt);

          mesh_pt->setup_boundary_element_info();

          // Turn off triangle refinement dump stuff (breaks Micromag
          // elements).
          checked_dynamic_cast<TriangleMeshBase*>(mesh_pt)->
            disable_triangulateio_restart();
        }
      else if(mesh_name == "ut_square" && nnode1d == 2)
        {
          mesh_pt = new TriangleMesh<TMagnetostaticFieldElement<2, 2> >
            ("./meshes/square." + to_string(refinement_level) + ".node",
             "./meshes/square." + to_string(refinement_level) + ".ele",
             "./meshes/square." + to_string(refinement_level) + ".poly",
             time_stepper_pt);

          // Turn off triangle refinement dump stuff (breaks Micromag
          // elements).
          checked_dynamic_cast<TriangleMeshBase*>(mesh_pt)->
            disable_triangulateio_restart();
        }
      else if(mesh_name == "st_cubeoid" && nnode1d == 2)
        {
          double lx = 1, ly = lx, lz = 3*lx;
          unsigned ny = nx, nz = std::ceil(lz/lx) * nx;
          mesh_pt = new SimpleCubicTetMesh<TMagnetostaticFieldElement<3, 2> >
            (nx, ny, nz, lx, ly, lz, time_stepper_pt);

          mesh_pt->setup_boundary_element_info();
        }
      else if(mesh_name == "ut_cubeoid" && nnode1d == 2)
        {
          mesh_pt = new TetgenMesh<TMagnetostaticFieldElement<3, 2> >
            ("./meshes/cubeoid." + to_string(refinement_level) + ".node",
             "./meshes/cubeoid." + to_string(refinement_level) + ".ele",
             "./meshes/cubeoid." + to_string(refinement_level) + ".face",
             time_stepper_pt);
        }
      else if(mesh_name == "ut_mumag4" && nnode1d == 2)
        {
          mesh_pt = new TetgenMesh<TMagnetostaticFieldElement<3, 2> >
            ("./meshes/mumag4." + to_string(refinement_level) + ".node",
             "./meshes/mumag4." + to_string(refinement_level) + ".ele",
             "./meshes/mumag4." + to_string(refinement_level) + ".face",
             time_stepper_pt);
        }
      else if(mesh_name == "sq_mumag4" && nnode1d == 2)
        {
          unsigned this_nx = refinement_level;

          mesh_pt = new SimpleCubicMesh<QMagnetostaticFieldElement<3, 2> >
            (5*this_nx, std::ceil(1.25*this_nx), 2, 500, 125, 3, time_stepper_pt);

          mesh_pt->setup_boundary_element_info();
        }
      else if(mesh_name == "sq_cubeoid" && nnode1d == 2)
        {
          double lx = 1, ly = lx, lz = 3*lx;
          mesh_pt = new SimpleCubicMesh<QMagnetostaticFieldElement<3, 2> >
            (nx, nx, int(lz/lx)*nx, lx, ly, lz, time_stepper_pt);

          mesh_pt->setup_boundary_element_info();
        }
      else if(mesh_name == "ut_sphere" && nnode1d == 2)
        {
          mesh_pt = new TetgenMesh<TMagnetostaticFieldElement<3, 2> >
            ("./meshes/sphere." + to_string(refinement_level) + ".node",
             "./meshes/sphere." + to_string(refinement_level) + ".ele",
             "./meshes/sphere." + to_string(refinement_level) + ".face",
             time_stepper_pt);
        }
      else
        {
          throw OomphLibError("Unrecognised mesh name " + mesh_name,
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      // Scale the mesh as requested
      scale_mesh(scaling_factor, mesh_pt);

      // Done: pass out the mesh pointer
      return mesh_pt;
    }



    /// \short Return a factory function which will create the appropriate
    /// "flux mesh" for the bulk element pointer given.
    GenericPoissonProblem::FluxMeshFactoryFctPt
    phi_1_flux_mesh_factory_factory(const FiniteElement* bulk_phi_1_ele_pt)
    {
      if(dynamic_cast<const TMagnetostaticFieldElement<2, 2>*>(bulk_phi_1_ele_pt) != 0)
        {
          return Factories::surface_mesh_factory<TMagnetostaticFieldFluxElement<2, 2> >;
        }
      else if(dynamic_cast<const TMagnetostaticFieldElement<3, 2>*>(bulk_phi_1_ele_pt) != 0)
        {
          return Factories::surface_mesh_factory<TMagnetostaticFieldFluxElement<3, 2> >;
        }

      else if(dynamic_cast<const QMagnetostaticFieldElement<2,2>*>(bulk_phi_1_ele_pt) != 0)
        {
          return Factories::surface_mesh_factory<QMagnetostaticFieldFluxElement<2,2> >;
        }
      else if(dynamic_cast<const QMagnetostaticFieldElement<3,2>*>(bulk_phi_1_ele_pt) != 0)
        {
          return Factories::surface_mesh_factory<QMagnetostaticFieldFluxElement<3,2> >;
        }

      else
        {
          throw OomphLibError("Unrecognised element type",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

    }

  }
}
