

#include "llg_problem.h"

// Meshes for mesh factory
#include "../../src/meshes/simple_rectangular_quadmesh.h"
#include "../../src/meshes/simple_rectangular_tri_mesh.h"
#include "../../src/meshes/simple_cubic_tet_mesh.h"
#include "../../src/meshes/simple_cubic_mesh.h"
#include "../../src/meshes/tetgen_mesh.h"
#include "../../src/meshes/triangle_mesh.h"


namespace oomph
{

  /// Function that does the real work of the constructors.
  void LLGProblem::build(Vector<Mesh*>& bulk_mesh_pts)
  {
#ifdef PARANOID
    if((applied_field_fct_pt() == 0)
       || (this->time_stepper_pt() == 0)
       || (Residual_calculator_pt == 0))
      {
        std::ostringstream error_msg;
        error_msg
          << "Must the following pointers to non-null values before calling build():\n"
          << "applied_field_fct_pt() (= " << applied_field_fct_pt() << ")\n"
          << "this->time_stepper_pt() (= " << this->time_stepper_pt() << ")\n"
          << "Residual_calculator_pt (= " << Residual_calculator_pt << ")"
          << std::endl;

        throw OomphLibError(error_msg.str(), OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

    // Call the underlying build to deal with adding meshes
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
          }
      }

    // Set up bem stuff if we are doing it
    if(Use_implicit_ms)
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
            Node* pinned_phi_1_node_pt = bulk_mesh_pts[msh]->get_some_non_boundary_node();
            pinned_phi_1_node_pt->pin(phi_1_index());
            pinned_phi_1_node_pt->set_value(phi_1_index(), 0.0);
          }

      }
    // Otherwise pin all phi and phi_1 dofs to zero
    else
      {
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



    // Build the global mesh
    this->build_global_mesh();

    // Number the equations
    this->assign_eqn_numbers();

    // // ??ds For if we want to swap solver for large dt. For now swap if
    // // we are using any iterative solver.
    // if(dynamic_cast<IterativeLinearSolver*>(linear_solver_pt()) != 0)
    //   {
    //     Swap_solver_large_dt = true;
    //   }
    // My_linear_solver_pt = linear_solver_pt();

    // Write out some stuff
    mag_parameters_pt()->output(std::cout);
    oomph_info << "LLG Number of equations: " << ndof() << std::endl;
    oomph_info << "Number of sub meshes: " << this->nsub_mesh() << std::endl;



    // Create bem_handler if we are doing fully implicit magnetostatics
    // ============================================================
    if(Use_implicit_ms)
      {
        Bem_handler_pt = new BoundaryElementHandler;

        // Figure out which element type we should use in the bem mesh
        // (based on the element type used in the bulk mesh) and store the
        // function needed to create them.
        Bem_handler_pt->Bem_element_factory = LLGFactories::
          bem_element_factory_factory(bulk_mesh_pts[0]->finite_element_pt(0));

        // Create an integration scheme ??ds move this outside somewhere...
        Bem_handler_pt->integration_scheme_pt() = LLGFactories::
          variable_order_integrator_factory(bulk_mesh_pts[0]->finite_element_pt(0));

        // Set indexes to look phi/phi1 variables
        Bem_handler_pt->set_input_index(phi_1_index());
        Bem_handler_pt->set_output_index(phi_index());

        // Loop over all meshes in problem adding to bem list
        for(unsigned msh=0, nmsh=bulk_mesh_pts.size(); msh<nmsh; msh++)
          {
            Bem_handler_pt->set_bem_all_boundaries(bulk_mesh_pts[msh]);
          }

        // ??ds no corners apart from boundary joins for now
        Bem_handler_pt->input_corner_data_pt() = 0;

        // Now build it
        Bem_handler_pt->build();

        // // Calculate the (initial) bem boundary conditions on phi
        // // ??ds might not need this? ok to wait until after a step?
        // maybe_update_bem_boundary_conditions();
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

  /// \short ??ds
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



  //======================================================================
  /// Set up the initial conditions
  //======================================================================
  void LLGProblem::
  set_initial_condition(const InitialM::InitialMFctPt initial_m_pt)
  {
    // Backup time in global Time object
    double backed_up_time=this->time_pt()->time();

    // Past history needs to be established for t=time0-deltat, ...
    // Then provide current values (at t=time0) which will also form
    // the initial guess for the first solve at t=time0+deltat

    // Get M indicies
    Vector<unsigned> m_index_micromag(3,0);
    MicromagEquations* elem_pt = dynamic_cast<MicromagEquations*>
      (this->mesh_pt()->element_pt(0));
    for(unsigned i=0; i<3; i++)
      {
        m_index_micromag[i] = elem_pt->m_index_micromag(i);
      }

    // Set continuous times at previous timesteps:
    int nprev_steps=this->time_stepper_pt()->nprev_values();
    Vector<double> prev_time(nprev_steps+1);
    for (int t=nprev_steps;t>=0;t--)
      {
        prev_time[t]=this->time_pt()->time(t);
      }

    // Loop over current & previous timesteps
    for (int t=nprev_steps;t>=0;t--)
      {
        // Continuous time
        double time = prev_time[t];
        std::cout << "setting IC at time =" << time << std::endl;

        // Loop over the nodes to set initial values everywhere
        for(unsigned n=0, nnd=mesh_pt()->nnode(); n<nnd; n++)
          {
            // Get initial value of m from inputs
            Vector<double> x(Dim, 0.0);
            mesh_pt()->node_pt(n)->position(t,x);
            Vector<double> m = initial_m_pt(time,x);

            // Set initial condition on m
            for(unsigned i=0; i<3; i++)
              {
                mesh_pt()->node_pt(n)->set_value
                  (t,m_index_micromag[i],m[i]);
              }
          }
      }

    // Reset backed up time for global timestepper
    this->time_pt()->time()=backed_up_time;

    // Do the energy calculations, don't try to calculate an effective damping
    calculate_energies(false);
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
                       unsigned nnode1d)
    {
      // Ignore case in mesh names
      const std::string mesh_name = to_lower(_mesh_name);

      // Make the mesh and store a pointer to it
      Mesh* mesh_pt;
      if(mesh_name == "sq_square" && nnode1d == 2)
        {
          double lx = 1.0;
          unsigned nx = 5 * std::pow(2, refinement_level);
          mesh_pt = new SimpleRectangularQuadMesh<QMicromagElement<2,2> >
            (nx, nx, lx, lx, time_stepper_pt);
        }
      else if(mesh_name == "ut_square" && nnode1d == 2)
        {
          mesh_pt = new TriangleMesh<TMicromagElement<2, 2> >
            ("./meshes/square." + to_string(refinement_level) + ".node",
             "./meshes/square." + to_string(refinement_level) + ".ele",
             "./meshes/square." + to_string(refinement_level) + ".poly",
             time_stepper_pt);
        }
      else if(mesh_name == "st_cubeoid" && nnode1d == 2)
        {
          // nmag cubeoid
          double lx = 30, ly = lx, lz = 100;
          unsigned nx = 2 * std::pow(2, refinement_level);
          unsigned ny = nx, nz = std::ceil(lz/lx) * nx;
          mesh_pt = new SimpleCubicTetMesh<TMicromagElement<3, 2> >
            (nx, ny, nz, lx, ly, lz, time_stepper_pt);
        }
      else if(mesh_name == "ut_cubeoid" && nnode1d == 2)
        {
          mesh_pt = new TetgenMesh<TMicromagElement<3, 2> >
            ("./meshes/cubeoid." + to_string(refinement_level) + ".node",
             "./meshes/cubeoid." + to_string(refinement_level) + ".ele",
             "./meshes/cubeoid." + to_string(refinement_level) + ".face",
             time_stepper_pt);
        }
      else if(mesh_name == "st_cubeoid" && nnode1d == 2)
        {
          double lx = 30, ly = lx, lz = 100;
          unsigned nx = std::pow(2, refinement_level);
          mesh_pt = new SimpleCubicTetMesh<TMicromagElement<3, 2> >
            (nx, nx, int(lz/lx)*nx, lx, ly, lz, time_stepper_pt);
        }
      else if(mesh_name == "sq_cubeoid" && nnode1d == 2)
        {
          double lx = 30, ly = lx, lz = 100;
          unsigned nx = std::pow(2, refinement_level);
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
}
