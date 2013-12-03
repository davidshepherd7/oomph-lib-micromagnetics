
#include "semi_implicit_problem.h"

// Quadratures for the quadrature factory
#include "variable_order_quadrature.h"


// Meshes for mesh factory
#include "../../src/meshes/simple_rectangular_quadmesh.h"
#include "../../src/meshes/simple_rectangular_tri_mesh.h"
#include "../../src/meshes/simple_cubic_tet_mesh.h"
#include "../../src/meshes/simple_cubic_mesh.h"
#include "../../src/meshes/tetgen_mesh.h"
#include "../../src/meshes/triangle_mesh.h"


namespace oomph
{
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
          mesh_pt = new SimpleRectangularQuadMesh<QSemiImplicitMicromagElement<2,2> >
            (nx, nx, lx, lx, time_stepper_pt);
        }
      else if(mesh_name == "ut_square" && nnode1d == 2)
        {
          mesh_pt = new TriangleMesh<TSemiImplicitMicromagElement<2, 2> >
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
          mesh_pt = new SimpleCubicTetMesh<TSemiImplicitMicromagElement<3, 2> >
            (nx, ny, nz, lx, ly, lz, time_stepper_pt);
        }
      else if(mesh_name == "ut_cubeoid" && nnode1d == 2)
        {
          mesh_pt = new TetgenMesh<TSemiImplicitMicromagElement<3, 2> >
            ("./meshes/cubeoid." + to_string(refinement_level) + ".node",
             "./meshes/cubeoid." + to_string(refinement_level) + ".ele",
             "./meshes/cubeoid." + to_string(refinement_level) + ".face",
             time_stepper_pt);
        }
      else if(mesh_name == "st_cubeoid" && nnode1d == 2)
        {
          double lx = 30, ly = lx, lz = 100;
          unsigned nx = std::pow(2, refinement_level);
          mesh_pt = new SimpleCubicTetMesh<TSemiImplicitMicromagElement<3, 2> >
            (nx, nx, int(lz/lx)*nx, lx, ly, lz, time_stepper_pt);
        }
      else if(mesh_name == "sq_cubeoid" && nnode1d == 2)
        {
          double lx = 30, ly = lx, lz = 100;
          unsigned nx = std::pow(2, refinement_level);
          mesh_pt = new SimpleCubicMesh<QSemiImplicitMicromagElement<3, 2> >
            (nx, nx, int(lz/lx)*nx, lx, ly, lz, time_stepper_pt);
        }
      else if(mesh_name == "ut_sphere" && nnode1d == 2)
        {
          mesh_pt = new TetgenMesh<TSemiImplicitMicromagElement<3, 2> >
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

      // For some reason we need to call this manually
      mesh_pt->setup_boundary_element_info();

      // Done: pass out the mesh pointer
      return mesh_pt;
    }


    /// \short Make a mesh of MagnetostaticField elements as specified by an
    /// input argument. Refined according to the given refinement level (in
    /// some way appropriate for that mesh type).
    Mesh* phi_mesh_factory(const std::string& _mesh_name,
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
          mesh_pt = new SimpleRectangularQuadMesh<QMagnetostaticFieldElement<2,2> >
                             (nx, nx, lx, lx, time_stepper_pt);
        }
      else if(mesh_name == "ut_square" && nnode1d == 2)
        {
          mesh_pt = new TriangleMesh<TMagnetostaticFieldElement<2, 2> >
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
          mesh_pt = new SimpleCubicTetMesh<TMagnetostaticFieldElement<3, 2> >
                             (nx, ny, nz, lx, ly, lz, time_stepper_pt);
        }
      else if(mesh_name == "ut_cubeoid" && nnode1d == 2)
        {
          mesh_pt = new TetgenMesh<TMagnetostaticFieldElement<3, 2> >
                             ("./meshes/cubeoid." + to_string(refinement_level) + ".node",
                              "./meshes/cubeoid." + to_string(refinement_level) + ".ele",
                              "./meshes/cubeoid." + to_string(refinement_level) + ".face",
                              time_stepper_pt);
        }
      else if(mesh_name == "st_cubeoid" && nnode1d == 2)
        {
          double lx = 30, ly = lx, lz = 100;
          unsigned nx = std::pow(2, refinement_level);
          mesh_pt = new SimpleCubicTetMesh<TMagnetostaticFieldElement<3, 2> >
                             (nx, nx, int(lz/lx)*nx, lx, ly, lz, time_stepper_pt);
        }
      else if(mesh_name == "sq_cubeoid" && nnode1d == 2)
        {
          double lx = 30, ly = lx, lz = 100;
          unsigned nx = std::pow(2, refinement_level);
          mesh_pt = new SimpleCubicMesh<QMagnetostaticFieldElement<3, 2> >
                             (nx, nx, int(lz/lx)*nx, lx, ly, lz, time_stepper_pt);
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

      // For some reason we need to call this manually
      mesh_pt->setup_boundary_element_info();

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

  /// \short Function to do the real work of the constructor.
  void SemiImplicitHybridMicromagneticsProblem::
  build(Vector<Mesh*>& llg_mesh_pts, Vector<Mesh*>& phi_mesh_pts,
        Vector<Mesh*>& phi_1_mesh_pts, bool pin_phi1)
  {
    // Throughout this function we need to be careful not to use each
    // problem's global mesh pointer until it as been built.


    // build llg (no extra work needed before the build)
    // ============================================================
    llg_sub_problem_pt()->build(llg_mesh_pts);


    // Set up phi_1 problem
    // ============================================================

    // phi_1 b.c.s are all Neumann but with flux determined elsewhere (by m).
    for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
      {
        Mesh* mesh_pt = phi_1_mesh_pts[msh];
        for(unsigned b=0, nb=mesh_pt->nboundary(); b < nb; b++)
          {
            phi_1_problem_pt()->add_neumann_boundary(mesh_pt, b, 0);
          }
      }

    // Doing this in base class instead
    // // Pin a node which isn't involved in the boundary element method (we
    // // have to pin something to avoid a singular Jacobian, can't be a
    // // boundary node or things will go wrong with BEM).
    // if(pin_phi1)
    //   {
    //     Node* pinned_phi_1_node_pt = phi_1_mesh_pts[0]->get_some_non_boundary_node();
    //     pinned_phi_1_node_pt->pin(0);
    //     pinned_phi_1_node_pt->set_value(0,0.0);
    //   }
    // else
    //   {
    //     std::cout << "Warning: not pinning phi1 at any point, technically J is singular..."
    //               << " you might be ok..."
    //               << std::endl;
    //   }

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
            error_msg << "Mesh nodes must be the same.";
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
                    std::ostringstream error_msg;
                    error_msg << "Mesh nodes must be in the same places.";
                    throw OomphLibError(error_msg.str(),
                                        OOMPH_CURRENT_FUNCTION,
                                        OOMPH_EXCEPTION_LOCATION);
                  }
              }
          }

        //??ds check that element e in each mesh has the nodes in the same
        //places?
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
            SemiImplicitMicromagEquations* m_ele_pt =
              checked_dynamic_cast<SemiImplicitMicromagEquations*>
              (llg_mesh_pt->element_pt(e));

            phi_1_ele_pt->set_micromag_element_pt(m_ele_pt);
            phi_ele_pt->set_micromag_element_pt(m_ele_pt);
            m_ele_pt->magnetostatic_field_element_pt() = phi_ele_pt;
          }

      }

    // BEM handler:
    // ============================================================
    // Construct the BEM (must be done before pinning phi values)

    // Set list of boundaries to use bem on
    for(unsigned msh=0, nmsh=phi_1_mesh_pts.size(); msh<nmsh; msh++)
      {
        for(unsigned b = 0; b < phi_1_mesh_pts[msh]->nboundary(); b++)
          {
            bem_handler_pt()->set_bem_boundary(b, phi_1_mesh_pts[msh]);
          }
      }

    // both indicies zero (because they are in separate problems)
    bem_handler_pt()->set_input_index(0);
    bem_handler_pt()->set_output_index(0);

    // Create an integration scheme ??ds move this outside somewhere...
    bem_handler_pt()->integration_scheme_pt() = LLGFactories::
      variable_order_integrator_factory(phi_1_mesh_pt()->finite_element_pt(0));

    // Only treat cases where the corners are rectangle-like for now
    bem_handler_pt()->input_corner_data_pt() = 0; //??Ds

    oomph_info << "Creating BEM handler" << std::endl;
    bem_handler_pt()->build();

    // Phi problem
    // ============================================================

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
  }



  void SemiImplicitHybridMicromagneticsProblem::
  set_initial_condition(const InitialM::InitialMFctPt initial_m_pt)
  {
    // Backup time in global Time object
    double backed_up_time=llg_sub_problem_pt()->time_pt()->time();

    // Past history needs to be established for t=time0-deltat, ...
    // Then provide current values (at t=time0) which will also form
    // the initial guess for the first solve at t=time0+deltat

    // Get M indicies
    Vector<unsigned> m_index_micromag(3,0);
    MicromagEquations* elem_pt = checked_dynamic_cast<MicromagEquations* >(llg_mesh_pt()->element_pt(0));
    for(unsigned i=0; i<3; i++)
      {
        m_index_micromag[i] = elem_pt->m_index_micromag(i);
      }

    // Find number of nodes in mesh
    unsigned num_nod = llg_mesh_pt()->nnode();

    // Set continuous times at previous timesteps:
    int nprev_steps=time_stepper_pt()->nprev_values();
    Vector<double> prev_time(nprev_steps+1);
    for (int t=nprev_steps;t>=0;t--)
      {
        prev_time[t]=llg_sub_problem_pt()->time_pt()->time(t);
      }

    // Loop over current & previous timesteps
    for (int t=nprev_steps;t>=0;t--)
      {
        // Continuous time
        double time = prev_time[t];
        std::cout << "setting IC at time =" << time << std::endl;

        // Loop over the nodes to set initial values everywhere
        for (unsigned n=0;n<num_nod;n++)
          {
            unsigned dim = llg_sub_problem_pt()->mesh_pt()->node_pt(n)->ndim();

            // Get initial value of m from inputs
            Vector<double> x(dim,0.0);
            llg_sub_problem_pt()->mesh_pt()->node_pt(n)->position(t,x);
            Vector<double> m = initial_m_pt(time,x);

            // Set initial condition on m
            for(unsigned i=0; i<3; i++)
              {
                llg_sub_problem_pt()->mesh_pt()->node_pt(n)->
                  set_value(t,m_index_micromag[i],m[i]);
              }
          }
      }

    // Reset backed up time for global timestepper
    llg_sub_problem_pt()->time_pt()->time()=backed_up_time;

    // Do the energy calculations, don't try to calculate an effective damping
    calculate_energies(false);
  }

  void SemiImplicitHybridMicromagneticsProblem::
  doc_solution_additional(std::ofstream &some_file) const
  {
    unsigned npts = 2;

    // Output llg solution into the main output file
    llg_sub_problem_pt()->mesh_pt()->output(some_file, npts);

    // Output the magnetostatic field data
    std::ofstream field_file((Doc_info.directory() + "/field"
                              + Doc_info.number_as_string() + ".dat").c_str());
    phi_problem_pt()->mesh_pt()->output(field_file, npts);
    field_file.close();

    std::ofstream phi1_file((Doc_info.directory() + "/phione"
                             + Doc_info.number_as_string() + ".dat").c_str());
    phi_1_problem_pt()->mesh_pt()->output(phi1_file, npts);
    phi1_file.close();
  }

  //============================================================
  //
  //============================================================
  Vector<double> SemiImplicitHybridMicromagneticsProblem::
  average_magnetostatic_field() const
  {
    const unsigned nodal_dim = checked_dynamic_cast<MagnetostaticFieldEquations*>
      (phi_problem_pt()->mesh_pt()->element_pt(0))->node_pt(0)->ndim();

    // Pick a point in the middle of the element
    const Vector<double> s(nodal_dim, 0.3);
    Vector<double> total_dphidx(nodal_dim,0.0);

    // Loop over all elements calculating the value in the middle of the element
    for(unsigned e=0, ne=phi_problem_pt()->mesh_pt()->nelement(); e < ne; e++)
      {
        MagnetostaticFieldEquations* ele_pt = checked_dynamic_cast<MagnetostaticFieldEquations*>
          (phi_problem_pt()->mesh_pt()->element_pt(e));

        // Get the shape function and eulerian coordinate derivative at
        // position s.
        unsigned n_node = ele_pt->nnode();
        Shape psi(n_node); DShape dpsidx(n_node,nodal_dim);
        ele_pt->dshape_eulerian(s,psi,dpsidx);

        // Interpolate grad phi
        Vector<double> interpolated_dphidx(nodal_dim,0.0);
        for(unsigned l=0;l<n_node;l++)
          {
            double phi_value = ele_pt->raw_nodal_value(l,0);
            for(unsigned i=0; i<nodal_dim; i++)
              {interpolated_dphidx[i] += phi_value*dpsidx(l,i);}
          }

        // Add this grad phi to the sum
        for(unsigned j=0; j<nodal_dim; j++)
          {
            total_dphidx[j] += interpolated_dphidx[j];
          }
      }

    // Divide sum by number of elements to get the average. Take the
    // negative to get the field.
    double nele = double(phi_problem_pt()->mesh_pt()->nelement());
    Vector<double> average_magnetostatic_field(3,0.0);
    for(unsigned j=0; j<nodal_dim; j++)
      {
        average_magnetostatic_field[j] = - total_dphidx[j] / nele;
      }

    return average_magnetostatic_field;
  }

}
