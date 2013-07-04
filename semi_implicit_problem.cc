
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
  /// (i.e. can be calculated at runtime) input and create a new instance
  /// of the appropriate object, using the new command (Factory Method
  /// design pattern).
  ///
  /// Typically these objects are passed straight into other classes and
  /// will be deleted by the destructor of that class. If not it is your
  /// responsibility to make sure the objects are deleted.
  namespace SemiImplicitFactories
  {

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
      Mesh* mesh_pt = 0;
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

      // For some reason we have to call this manually...
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
      Mesh* mesh_pt = 0;
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

      // For some reason we have to call this manually...
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

  }

  /// \short Function to do the real work of the constructor.
  void SemiImplicitHybridMicromagneticsProblem::build(bool pin_phi1)
  {
    // Set up phi_1 problem
    // ============================================================
    {

      Vector<unsigned> neu_bound;
      for(unsigned b=0, nb=phi_1_mesh_pt()->nboundary(); b < nb; b++)
        {
          neu_bound.push_back(b);
        }

      // phi_1 b.c.s are all Neumann but with flux determined elsewhere (by m)
      phi_1_problem_pt()->set_neumann_boundaries(neu_bound, 0);

      if(pin_phi1)
        {
          // Pin a node which isn't involved in the boundary element method (we
          // have to pin something to avoid a singular Jacobian, can't be a
          // boundary node or things will go wrong with BEM).
          Node* pinned_phi_1_node_pt = phi_1_mesh_pt()->get_some_non_boundary_node();
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
      phi_1_problem_pt()->build();

      // Things will go wrong if the nodes of all meshes are not in
      // the same place:
#ifdef PARANOID
      if((phi_1_mesh_pt()->nnode() != phi_mesh_pt()->nnode())
         || (phi_1_mesh_pt()->nnode() != llg_mesh_pt()->nnode()))
        {
          std::ostringstream error_msg;
          error_msg << "Mesh nodes must be the same.";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      unsigned dim = phi_mesh_pt()->node_pt(0)->ndim();
      for(unsigned j=0; j<dim; j++)
        {
          for(unsigned nd=0, nnode= phi_1_mesh_pt()->nnode(); nd<nnode; nd++)
            {
              if((phi_1_mesh_pt()->node_pt(nd)->x(j) !=
                  phi_mesh_pt()->node_pt(nd)->x(j))
                 ||
                 (phi_1_mesh_pt()->node_pt(nd)->x(j) !=
                  llg_mesh_pt()->node_pt(nd)->x(j)))
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
      for(unsigned e=0, ne=phi_1_mesh_pt()->nelement(); e < ne; e++)
        {
          MagnetostaticFieldEquations* ele_pt =
            checked_dynamic_cast<MagnetostaticFieldEquations*>
            (phi_1_mesh_pt()->element_pt(e));

          MicromagEquations* m_ele_pt = checked_dynamic_cast<MicromagEquations*>
            (llg_mesh_pt()->element_pt(e));

          ele_pt->set_micromag_element_pt(m_ele_pt);
        }

    }

    // BEM handler:
    // ============================================================
    {
      // Construct the BEM (must be done before pinning phi values)
      bem_handler_pt()->set_bem_all_boundaries(phi_1_mesh_pt());

      // both zero because they are in seperate problems
      bem_handler_pt()->set_input_index(0);
      bem_handler_pt()->set_output_index(0);

      // Create an integration scheme
      bem_handler_pt()->integration_scheme_pt() =
        SemiImplicitFactories::
        variable_order_integrator_factory(phi_1_mesh_pt()->finite_element_pt(0));

      bem_handler_pt()->input_corner_data_pt() = 0; //??Ds


      oomph_info << "Creating BEM handler" << std::endl;
      bem_handler_pt()->build();
    }

    // Now we can set up phi problem
    // ============================================================

    unsigned nboundary = phi_mesh_pt()->nboundary();
    Phi_boundary_values_pts.assign(nboundary, 0);
    for(unsigned b=0; b < nboundary; b++)
      {
        // Phi is determined by BEM
        LinearAlgebraDistribution* dist_pt =
          new LinearAlgebraDistribution(MPI_Helpers::communicator_pt(),
                                        phi_mesh_pt()->nboundary_node(b), false);

        Phi_boundary_values_pts[b] = new DoubleVector(dist_pt);
        phi_problem_pt()->set_dirichlet_boundary_by_vector(b, Phi_boundary_values_pts[b]);
      }
    phi_problem_pt()->build();

    // Assign micromagnetics element pointers
    // ??ds dodgy...
    for(unsigned e=0, ne=phi_mesh_pt()->nelement(); e < ne; e++)
      {
        MagnetostaticFieldEquations* ele_pt = checked_dynamic_cast<MagnetostaticFieldEquations*>
          (phi_mesh_pt()->element_pt(e));

        MicromagEquations* m_ele_pt = checked_dynamic_cast<MicromagEquations*>
          (llg_mesh_pt()->element_pt(e));

        ele_pt->set_micromag_element_pt(m_ele_pt);
      }


    // LLG problem:
    // ============================================================

    //??ds while we still have phi in MM elements pin them all
    unsigned phi_index = llg_element_pt()->phi_index_micromag();
    unsigned phi_1_index = llg_element_pt()->phi_1_index_micromag();
    for(unsigned nd=0, nnode=llg_mesh_pt()->nnode(); nd<nnode; nd++)
      {
        Node* nd_pt = llg_mesh_pt()->node_pt(nd);
        nd_pt->pin(phi_index);
        nd_pt->pin(phi_1_index);
      }

    // Assign phi element pointers
    // ??ds dodgy...
    for(unsigned e=0, ne=llg_mesh_pt()->nelement(); e < ne; e++)
      {
        SemiImplicitMicromagEquations* ele_pt
          = checked_dynamic_cast<SemiImplicitMicromagEquations*>
          (llg_mesh_pt()->element_pt(e));
        ele_pt->magnetostatic_field_element_pt() =
          checked_dynamic_cast<MagnetostaticFieldEquations*>(phi_mesh_pt()->element_pt(e));
      }


    // Build the LLG part of the problem. ??ds
    llg_sub_problem_pt()->build();
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
