

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
  void LLGProblem::build()
  {
#ifdef PARANOID
    if((bulk_mesh_pt() == 0)
       || (applied_field_fct_pt() == 0)
       || (this->time_stepper_pt() == 0)
       || (Residual_calculator_pt == 0))
      {
        std::ostringstream error_msg;
        error_msg
          << "Must the following pointers to non-null values before calling build():\n"
          << "bulk_mesh_pt() (= " << bulk_mesh_pt() << ")\n"
          << "applied_field_fct_pt() (= " << applied_field_fct_pt() << ")\n"
          << "this->time_stepper_pt() (= " << this->time_stepper_pt() << ")\n"
          << "Residual_calculator_pt (= " << Residual_calculator_pt << ")"
          << std::endl;

        throw OomphLibError(error_msg.str(), OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

    // Write out parameters data.
    mag_parameters_pt()->output(std::cout);

    // Cache the problem dimension
    this->Dim = this->ele_pt()->nodal_dimension();

    // Boundary conditions - our finite element discretisation requires
    // residual addition of m x dm/dn along the boundary (due to need to
    // reduce the order of the laplacian on m in exchange field).

    // // Create mesh of MicromagFluxElement<MicromagEquations> to add boundary
    // // contribution (if any).
    // surface_exchange_mesh_pt() = new Mesh;
    // for(unsigned b=0, nb=bulk_mesh_pt()->nboundary(); b < nb; b++)
    //   {
    //     create_surface_exchange_elements(b);
    //   }

    // // Figure out which residual function to use
    // if(use_llg_residual)
    //   {
    //     Residual_calculator_pt = new LLResidualCalculator;
    //   }
    // else
    //   {
    //     Residual_calculator_pt = new LLGResidualCalculator;
    //   }

    // Finish off elements
    for(unsigned i=0; i< bulk_mesh_pt()->nelement(); i++)
      {
        MicromagEquations* elem_pt = checked_dynamic_cast<MicromagEquations*>
          (bulk_mesh_pt()->element_pt(i));

        // Set whether the Jacobian should be finite differenced
        elem_pt->Use_fd_jacobian = Use_fd_jacobian;

        // Set values for magnetic parameters
        elem_pt->magnetic_parameters_pt() = mag_parameters_pt();

        // Set pointer for an applied field
        elem_pt->applied_field_pt() = applied_field_fct_pt();

        // Set the residual calculation function
        elem_pt->Residual_calculator_pt = Residual_calculator_pt;
      }

    // Pin all phi dofs...??ds remove them from element?
    for(unsigned nd=0, nnode=bulk_mesh_pt()->nnode(); nd<nnode; nd++)
      {
        Node* nd_pt = bulk_mesh_pt()->node_pt(nd);

        nd_pt->pin(phi_index());
        nd_pt->pin(phi_1_index());
        nd_pt->set_value(phi_index(),0.0);
        nd_pt->set_value(phi_1_index(),0.0);
      }

    // Build the global mesh
    this->add_sub_mesh(bulk_mesh_pt());
    // add_sub_mesh(surface_exchange_mesh_pt());
    this->build_global_mesh();

    // Set up base class (at the moment just does meshes in block
    // preconditioner).
    MyProblem::build();

    // // ??ds For if we want to swap solver for large dt. For now swap if
    // // we are using any iterative solver.
    // if(dynamic_cast<IterativeLinearSolver*>(linear_solver_pt()) != 0)
    //   {
    //     Swap_solver_large_dt = true;
    //   }
    // My_linear_solver_pt = linear_solver_pt();

    // Do equation numbering
    oomph_info << "LLG Number of equations: " << this->assign_eqn_numbers() << std::endl;
    oomph_info << "Number of sub meshes: " << this->nsub_mesh() << std::endl;
  }


  /// \short Error for adaptive timestepper (rms of nodal error determined by
  /// comparison with explicit timestepper result).
  double LLGProblem::global_temporal_error_norm()
  {
    double global_error = 0.0;

    //Find out how many nodes there are in the problem
    unsigned n_node = bulk_mesh_pt()->nnode();

    //Loop over the nodes and calculate the estimated error in the values
    for(unsigned i=0;i<n_node;i++)
      {
        for(unsigned j=0; j<3; j++)
          {
            // Get timestepper's error estimate for this direction of m
            // at this point.
            double error = bulk_mesh_pt()->node_pt(i)->time_stepper_pt()->
              temporal_error_in_value(bulk_mesh_pt()->node_pt(i), m_index(j));

            //Add the square of the individual error to the global error
            global_error += error*error;
          }
      }

    // Divide by the number of data points
    global_error /= 3*double(n_node);

    return std::sqrt(global_error);
  }


  /// \short Loop over all nodes in bulk mesh and get magnetisations
  Vector<Vector<double> > LLGProblem::get_nodal_magnetisations(unsigned i_time) const
  {
    unsigned nnode = bulk_mesh_pt()->nnode();
    Vector< Vector<double> > m_list(nnode, Vector<double>(3, 0.0));

    for(unsigned nd=0; nd<nnode; nd++)
      {
        for(unsigned j=0; j<3; j++)
          {
            m_list[nd][j] = bulk_mesh_pt()->node_pt(nd)->value(i_time, m_index(j));
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

    // Sum the difference over all elements
    for(unsigned e=0, ne=bulk_mesh_pt()->nelement(); e < ne; e++)
      {
        MicromagEquations* ele_pt = dynamic_cast<MicromagEquations* >
          (bulk_mesh_pt()->element_pt(e));

        Vector<double> numerical_m(3,0.0);
        ele_pt->interpolated_m_micromag(s,numerical_m);

        Vector<double> x(dim(),0.0);
        ele_pt->interpolated_x(s,x);
        double t = this->time();
        Vector<double> exact_m = fct_pt(t, x);

        for(unsigned j=0; j<3; j++)
          {
            diff += std::abs(numerical_m[j] - exact_m[j]);
          }
      }

    // Divide to get the mean
    diff /= (3 * bulk_mesh_pt()->nelement());

    return diff;
  }


  double LLGProblem::mean_norm_m_error() const
  {
    double temp_error = 0;

    unsigned nnode = bulk_mesh_pt()->nnode();
    for(unsigned nd=0; nd<nnode; nd++)
      {
        Node* nd_pt = bulk_mesh_pt()->node_pt(nd);
        Vector<double> m_values(3,0.0);
        for(unsigned j=0; j<3; j++) m_values[j] = nd_pt->value(m_index(j));

        double m_2norm = VectorOps::two_norm(m_values);
        temp_error += std::abs(m_2norm - 1);
        // std::cout << std::abs(m_2norm - 1) << std::endl;
      }

    return temp_error/double(nnode);
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
    double temp_stddev = 0.0;
    unsigned nnode=bulk_mesh_pt()->nnode();
    for(unsigned nd=0; nd<nnode; nd++)
      {
        Node* nd_pt = bulk_mesh_pt()->node_pt(nd);
        Vector<double> m_values(3,0.0);
        for(unsigned j=0; j<3; j++) m_values[j] = nd_pt->value(m_index(j));

        double m_2norm = VectorOps::two_norm(m_values);
        temp_stddev += std::pow( m_2norm - 1 - m_error_avg, 2);
      }

    temp_stddev /= nnode;
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

    // Find number of nodes in mesh
    unsigned num_nod = this->mesh_pt()->nnode();

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
        for (unsigned n=0;n<num_nod;n++)
          {
            unsigned dim = this->mesh_pt()->node_pt(n)->ndim();

            // Get initial value of m from inputs
            Vector<double> x(dim,0.0);
            this->mesh_pt()->node_pt(n)->position(t,x);
            Vector<double> m = initial_m_pt(time,x);

            // Set initial condition on m
            for(unsigned i=0; i<3; i++)
              this->mesh_pt()->node_pt(n)->set_value(t,m_index_micromag[i],m[i]);
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
      Mesh* mesh_pt = 0;
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

      // Nnode1d = 3
      else if(mesh_name == "sq_square" && nnode1d == 3)
        {
          double lx = 1.0;
          unsigned nx = 5 * std::pow(2, refinement_level);
          mesh_pt = new SimpleRectangularQuadMesh<QMicromagElement<2,3> >
            (nx, nx, lx, lx, time_stepper_pt);
        }
      else if(mesh_name == "ut_square" && nnode1d == 3)
        {
          mesh_pt = new TriangleMesh<TMicromagElement<2, 3> >
            ("./meshes/square." + to_string(refinement_level) + ".node",
             "./meshes/square." + to_string(refinement_level) + ".ele",
             "./meshes/square." + to_string(refinement_level) + ".poly",
             time_stepper_pt);
        }
      else if(mesh_name == "st_cubeoid" && nnode1d == 3)
        {
          // nmag cubeoid
          double lx = 30, ly = lx, lz = 100;
          unsigned nx = 2 * std::pow(2, refinement_level);
          unsigned ny = nx, nz = std::ceil(lz/lx) * nx;
          mesh_pt = new SimpleCubicTetMesh<TMicromagElement<3, 3> >
            (nx, ny, nz, lx, ly, lz, time_stepper_pt);
        }
      else if(mesh_name == "ut_cubeoid" && nnode1d == 3)
        {
          mesh_pt = new TetgenMesh<TMicromagElement<3, 3> >
            ("./meshes/cubeoid." + to_string(refinement_level) + ".node",
             "./meshes/cubeoid." + to_string(refinement_level) + ".ele",
             "./meshes/cubeoid." + to_string(refinement_level) + ".face",
             time_stepper_pt);
        }
      else if(mesh_name == "st_cubeoid" && nnode1d == 3)
        {
          double lx = 30, ly = lx, lz = 100;
          unsigned nx = std::pow(2, refinement_level);
          mesh_pt = new SimpleCubicTetMesh<TMicromagElement<3, 3> >
            (nx, nx, int(lz/lx)*nx, lx, ly, lz, time_stepper_pt);
        }
      else if(mesh_name == "sq_cubeoid" && nnode1d == 3)
        {
          double lx = 30, ly = lx, lz = 100;
          unsigned nx = std::pow(2, refinement_level);
          mesh_pt = new SimpleCubicMesh<QMicromagElement<3, 3> >
            (nx, nx, int(lz/lx)*nx, lx, ly, lz, time_stepper_pt);
        }
      else if(mesh_name == "ut_sphere" && nnode1d == 3)
        {
          mesh_pt = new TetgenMesh<TMicromagElement<3, 3> >
            ("./meshes/sphere." + to_string(refinement_level) + ".node",
             "./meshes/sphere." + to_string(refinement_level) + ".ele",
             "./meshes/sphere." + to_string(refinement_level) + ".face",
             time_stepper_pt);
        }

      // Otherwise fail
      else
        {
          throw OomphLibError("Unrecognised mesh name = " + mesh_name +
                              " or nnode1d = " + to_string(nnode1d),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      // For some reason we have to call this manually...
      mesh_pt->setup_boundary_element_info();

      // Done: pass out the mesh pointer
      return mesh_pt;
    }

  /// Pick and create a residual calculator to use
  ResidualCalculator* residual_calculator_factory(const std::string& residual)
    {
      if(residual == "llg")
        return new LLGResidualCalculator;
      else if(residual == "ll")
        return new LLResidualCalculator;
      else
        throw OomphLibError("Unrecognised residual "+residual,
                            OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
    }

  }
}
