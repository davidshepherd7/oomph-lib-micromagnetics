#ifndef OOMPH_SEMI_IMPLICIT_PROBLEM_H
#define OOMPH_SEMI_IMPLICIT_PROBLEM_H

/*
  description of file goes here
*/

#include "generic.h"
#include "./boundary_element_handler.h"
#include "./generic_poisson_problem.h"
#include "./implicit_llg_problem.h"

#include "./micromagnetics_boundary_element.h"
#include "./magnetostatic_field_flux_element.h"
#include "./micromag.h"

using namespace oomph;

namespace oomph
{

  class SemiImplicitHybridMicromagneticsProblem :
    public ImplicitLLGProblem
  {
  public:

    // Default constructor, who knows what will happen here... ??ds
    SemiImplicitHybridMicromagneticsProblem() :
      Bem_handler_pt(0), Phi_1_problem_pt(), Phi_problem_pt(),
      Phi_boundary_values_pts()
    {}


    void build(bool pin_phi1 = true)
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
            MagnetostaticFieldEquations* ele_pt = checked_dynamic_cast<MagnetostaticFieldEquations*>
              (phi_1_mesh_pt()->element_pt(e));

            MicromagEquations* m_ele_pt = checked_dynamic_cast<MicromagEquations*>
              (llg_mesh_pt()->element_pt(e));

            ele_pt->set_micromag_element_pt( m_ele_pt);
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
          Factories::variable_order_integrator_factory(phi_1_mesh_pt()->finite_element_pt(0));

        bem_handler_pt()->input_corner_data_pt() = 0; //??Ds

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
      for(unsigned nd=0, nnode=llg_mesh_pt()->nnode(); nd<nnode; nd++)
        {
          Node* nd_pt = llg_mesh_pt()->node_pt(nd);

          unsigned phi_index = llg_element_pt()->phi_index_micromag();
          unsigned phi_1_index = llg_element_pt()->phi_1_index_micromag();

          nd_pt->pin(phi_index);
          nd_pt->pin(phi_1_index);
        }

      // // Get timestepper from mesh
      // TimeStepper* ts_pt = llg_mesh_pt->node_pt(0)->time_stepper_pt();
      // llg_sub_problem_pt()->add_time_stepper_pt(ts_pt);

      // // Magnetic parameters
      // llg_sub_problem_pt()->mag_parameters_pt()->set_nmag_rectangle();
      // llg_sub_problem_pt()->applied_field_fct_pt() = applied_field_pt;

      // llg_sub_problem_pt()->build();

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
      if(llg_sub_problem_pt() == this)
        {
          ImplicitLLGProblem::build();
        }
      else
        {
          llg_sub_problem_pt()->build();
        }

    };

    /// Destructor
    ~SemiImplicitHybridMicromagneticsProblem()
    {
      // Kill boundary value storage vectors
      for(unsigned j=0; j<Phi_boundary_values_pts.size(); j++)
        {
          delete Phi_boundary_values_pts[j];
        }
    }

    /// Replacement for "newton_solve()" that does a few different solves.
    double semi_implicit_step(const double &dt, const double eps=0.0)
    {

      // solve for phi1
      std::cout << "solving phi1" << std::endl;
      phi_1_problem_pt()->newton_solve();

      // update boundary values of phi
      std::cout << "solving BEM" << std::endl;
      Bem_handler_pt->get_bem_values(Phi_boundary_values_pts);

      // solve for phi
      std::cout << "solving phi" << std::endl;
      phi_problem_pt()->newton_solve();

      // solve for m
      if(eps != 0.0)
        {
          std::cout << "solving LLG adaptively" << std::endl;
          return llg_sub_problem_pt()->adaptive_unsteady_newton_solve(dt, eps);
        }
      else
        {
          std::cout << "solving LLG" << std::endl;
          llg_sub_problem_pt()->unsteady_newton_solve(dt);
        }
      return dt;
    }

    /// Set up an initial M
    void set_initial_condition(const InitialM::InitialMFctPt initial_m_pt);

    /// Initialise timestep: only llg problem has a timestep.
    void initialise_dt(const double &dt)
    {llg_sub_problem_pt()->initialise_dt(dt);}

    /// Output
    void doc_solution();

    void average_magnetostatic_field(Vector<double> &average_magnetostatic_field) const;

    /// Access functions
    // =================================================================

    /// Get pointer to an LLG element for looking up info. All elements in
    /// mesh should be the same type otherwise preconditioning framework
    /// will fail so this should be safe.
    MicromagEquations* llg_element_pt() const
    {return checked_dynamic_cast<MicromagEquations*>
        (llg_sub_problem_pt()->bulk_mesh_pt()->element_pt(0));}

    /// Get access to magnetic parameters - only relevant in LLG problem so
    /// return the pointer from that problem.
    MagneticParameters* mag_parameters_pt()
    {return llg_sub_problem_pt()->mag_parameters_pt();}

    /// \short Const access function for LLG_problem.
    const ImplicitLLGProblem* llg_sub_problem_pt() const
    {return this;}

    /// \short Non-const acess function for LLG_problem.
    ImplicitLLGProblem* llg_sub_problem_pt() {return this;}

    const Mesh* bem_mesh_pt() const {return Bem_handler_pt->bem_mesh_pt();}

    const Mesh* llg_mesh_pt() const
    {
#ifdef PARANOID
      if(llg_sub_problem_pt()->bulk_mesh_pt() == 0)
        {
          std::string error_msg = "LLG mesh pointer is null!";
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return llg_sub_problem_pt()->bulk_mesh_pt();
    }

    const Mesh* phi_1_mesh_pt() const {return phi_1_problem_pt()->bulk_mesh_pt();}
    const Mesh* phi_mesh_pt() const {return phi_problem_pt()->bulk_mesh_pt();}

    const DenseMatrix<double>* bem_matrix_pt() const {return Bem_handler_pt->bem_matrix_pt();}

    // /// Set the list of sharp corners in the mesh to be a rectangle.
    // void set_rectangular_corners()
    // {Bem_handler_pt->Mesh_angles_type = "rectangular";}

    BoundaryElementHandler* &bem_handler_pt() {return Bem_handler_pt;}


    void set_phi_1_problem_pt(GenericPoissonProblem* p)
    { Phi_1_problem_pt = p;}

    GenericPoissonProblem* phi_1_problem_pt() const
    {
#ifdef PARANOID
      if(Phi_1_problem_pt == 0)
        {
          std::string error_msg = "Phi 1 problem pointer is null!";
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return Phi_1_problem_pt;
    }


    void set_phi_problem_pt(GenericPoissonProblem* p)
    { Phi_problem_pt = p;}

    GenericPoissonProblem* phi_problem_pt() const
    {
#ifdef PARANOID
      if(Phi_problem_pt == 0)
        {
          std::string error_msg = "Phi problem pointer is null!";
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return Phi_problem_pt;
    }


  private:

    /// Object to provide all BEM related capabilities.
    BoundaryElementHandler* Bem_handler_pt;

    /// Problem to solve for phi_1 (the BEM pre-calculation).
    GenericPoissonProblem* Phi_1_problem_pt;

    /// Problem to solve for phi (the magnetostatic potential).
    GenericPoissonProblem* Phi_problem_pt;

    // /// Problem to solve for the magnetisation change.
    // ImplicitLLGProblem LLG_problem;

    /// Intermediate storage for results of bem (ideally we would have it
    /// call a function to get the boundary values filled in but c++ member
    /// functions pointers are useless...)
    Vector<DoubleVector*> Phi_boundary_values_pts;
  };


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
  doc_solution()
  {
    using namespace StringConversion;

    unsigned npts = 2;

    // Output llg solution
    std::ofstream soln_file((Doc_info.directory() + "/soln"
                             + Doc_info.number_as_string() + ".dat").c_str());
    llg_sub_problem_pt()->mesh_pt()->output(soln_file,npts);
    soln_file.close();

    // Output the magnetostatic field data
    std::ofstream field_file((Doc_info.directory() + "/field"
                              + Doc_info.number_as_string() + ".dat").c_str());
    for(unsigned e=0, ne=phi_problem_pt()->mesh_pt()->nelement(); e < ne; e++)
      {
        MagnetostaticFieldEquations* ele_pt = checked_dynamic_cast<MagnetostaticFieldEquations*>
          (phi_problem_pt()->mesh_pt()->element_pt(e));
        ele_pt->output(field_file,npts);
      }
    field_file.close();

    std::ofstream phi1_file((Doc_info.directory() + "/phione"
                             + Doc_info.number_as_string() + ".dat").c_str());
    phi_1_problem_pt()->mesh_pt()->output(phi1_file,npts);
    phi1_file.close();

    // Write average magnetisations to a file
    Vector<double> m = llg_sub_problem_pt()->mean_magnetisation();

    std::ofstream avgs((Doc_info.directory() +"/averages").c_str(),
                       std::ios::app);
    avgs << llg_sub_problem_pt()->time();
    for(unsigned j=0; j<3; j++) avgs << " " << m[j];
    avgs << std::endl;
    avgs.close();


    // Write average field to a file
    Vector<double> hms;
    average_magnetostatic_field(hms);
    std::ofstream field_avgs((Doc_info.directory() +"/field_averages").c_str(),
                             std::ios::app);
    field_avgs << llg_sub_problem_pt()->time();
    for(unsigned j=0; j<3; j++) field_avgs << " " << hms[j];
    field_avgs << std::endl;
    field_avgs.close();


    // Get average (and standard deviation) of |m| - 1 and |m|.dm/dn
    double m_error_avg(0), m_error_stddev(0), orthogonality_error_avg(0),
      orthogonality_error_stddev(0);
    llg_sub_problem_pt()->norm_m_error(m_error_avg, m_error_stddev);
    // llg_sub_problem_pt()->orthogonality_m_error
      // (orthogonality_error_avg, orthogonality_error_stddev);

    // Write them to file
    std::ofstream errors((Doc_info.directory()+"/errors").c_str(),std::ios::app);
    errors << llg_sub_problem_pt()->time()
           << " " << m_error_avg
           << " " << m_error_stddev
           << " " << orthogonality_error_avg
           << " " << orthogonality_error_stddev
           << std::endl;
    errors.close();


    // Output convergence data if we have it
    if(llg_sub_problem_pt()->record_convergence_data())
      {
        llg_sub_problem_pt()->convergence_data_pt()
          ->output_this_newton_step("convergence_data");
      }

    // Finally increment the label ready for next time
    Doc_info.number()++;
  }

  //============================================================
  //
  //============================================================
  void SemiImplicitHybridMicromagneticsProblem::
  average_magnetostatic_field(Vector<double> &average_magnetostatic_field) const
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
    average_magnetostatic_field.assign(3,0.0);
    for(unsigned j=0; j<nodal_dim; j++)
      {
        average_magnetostatic_field[j] = - total_dphidx[j] / nele;
      }
  }

  // void SemiImplicitHybridMicromagneticsProblem::
  // doc_magnetostatic_field(DocInfo &doc_info) const
  // {
  //   // Number of plot points
  //   unsigned npts;
  //   npts=2;

  //   std::ofstream soln_file((doc_info.directory() + "/soln"
  //                            + doc_info.number_as_string() + ".dat").c_str());

  //   phi_problem_pt()->doc_solution(soln

  //   soln_file.close();
  // }


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
    /// \short Make a mesh of Micromag elements as specified by an
    /// input argument. Refined according to the given refinement level (in
    /// some way appropriate for that mesh type).
    Mesh* llg_mesh_factory(const std::string& _mesh_name,
                           int refinement_level,
                           TimeStepper* time_stepper_pt,
                           unsigned nnode1d = 2)
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
            ("../meshes/square." + to_string(refinement_level) + ".node",
             "../meshes/square." + to_string(refinement_level) + ".ele",
             "../meshes/square." + to_string(refinement_level) + ".poly",
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
            ("../meshes/cubeoid." + to_string(refinement_level) + ".node",
             "../meshes/cubeoid." + to_string(refinement_level) + ".ele",
             "../meshes/cubeoid." + to_string(refinement_level) + ".face",
             time_stepper_pt);
        }
      else if(mesh_name == "st_cubeoid" && nnode1d == 2)
        {
          double lx = 30, ly = lx, lz = 100;
          unsigned nx = std::pow(2, refinement_level);
          mesh_pt = new SimpleCubicTetMesh<TSemiImplicitMicromagElement<3, 2> >
            (nx, nx, int(lz/lx)*nx, lx, ly, lz, time_stepper_pt);
        }
      else if(mesh_name == "ut_sphere" && nnode1d == 2)
        {
          mesh_pt = new TetgenMesh<TSemiImplicitMicromagElement<3, 2> >
            ("../meshes/sphere." + to_string(refinement_level) + ".node",
             "../meshes/sphere." + to_string(refinement_level) + ".ele",
             "../meshes/sphere." + to_string(refinement_level) + ".face",
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
                           unsigned nnode1d = 2)
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
            ("../meshes/square." + to_string(refinement_level) + ".node",
             "../meshes/square." + to_string(refinement_level) + ".ele",
             "../meshes/square." + to_string(refinement_level) + ".poly",
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
            ("../meshes/cubeoid." + to_string(refinement_level) + ".node",
             "../meshes/cubeoid." + to_string(refinement_level) + ".ele",
             "../meshes/cubeoid." + to_string(refinement_level) + ".face",
             time_stepper_pt);
        }
      else if(mesh_name == "st_cubeoid" && nnode1d == 2)
        {
          double lx = 30, ly = lx, lz = 100;
          unsigned nx = std::pow(2, refinement_level);
          mesh_pt = new SimpleCubicTetMesh<TMagnetostaticFieldElement<3, 2> >
            (nx, nx, int(lz/lx)*nx, lx, ly, lz, time_stepper_pt);
        }
      else if(mesh_name == "ut_sphere" && nnode1d == 2)
        {
          mesh_pt = new TetgenMesh<TMagnetostaticFieldElement<3, 2> >
            ("../meshes/sphere." + to_string(refinement_level) + ".node",
             "../meshes/sphere." + to_string(refinement_level) + ".ele",
             "../meshes/sphere." + to_string(refinement_level) + ".face",
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


    typedef MicromagBEMElementEquations*
    (*BEMElementFactoryFctPt)(FiniteElement* const, const int&);


    /// \short very simple function: create a new face element of type
    /// ELEMENT.
    template<class ELEMENT>
    MicromagBEMElementEquations* bem_element_factory(FiniteElement* ele,
                                                     const int& face)
      {
        return new ELEMENT(ele, face);
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

  class SemiImplicitMMArgs : public MyCliArgs
  {
  public:

    SemiImplicitMMArgs() : llg_mesh_pt(0),
                           phi_1_mesh_pt(0), phi_mesh_pt(0),
                           initial_m_fct_pt(0), h_app_fct_pt(0),
                           phi_1_flux_mesh_factory_fct_pt(0),
                           bem_element_factory_fct_pt(0)
    {}


    virtual void set_flags()
    {
      MyCliArgs::set_flags();

      CommandLineArgs::specify_command_line_flag("-mesh", &mesh_name);
      mesh_name = "sq_square";

      CommandLineArgs::specify_command_line_flag("-initm", &initial_m_name);
      initial_m_name = "z";

      CommandLineArgs::specify_command_line_flag("-happ", &h_app_name);
      h_app_name = "minus_z";
    }


    void run_factories()
    {
      MyCliArgs::run_factories();

      to_lower(mesh_name);
      to_lower(initial_m_name);
      to_lower(h_app_name);

      llg_mesh_pt = SemiImplicitFactories::llg_mesh_factory
        (mesh_name, refinement, time_stepper_pt);

      // Make the two phi meshes
      phi_mesh_pt = SemiImplicitFactories::phi_mesh_factory
        (mesh_name, refinement, time_stepper_pt);
      phi_1_mesh_pt = SemiImplicitFactories::phi_mesh_factory
        (mesh_name, refinement, time_stepper_pt);

      // Pick the m and applied field function pointers
      initial_m_fct_pt = InitialM::initial_m_factory(initial_m_name);
      h_app_fct_pt = HApp::h_app_factory(h_app_name);

      // Pick the factory function for creating the phi 1 surface mesh
      phi_1_flux_mesh_factory_fct_pt =
        SemiImplicitFactories::phi_1_flux_mesh_factory_factory
        (phi_1_mesh_pt->finite_element_pt(0));

      // Pick the factory function for creating the BEM elements
      bem_element_factory_fct_pt =
        SemiImplicitFactories::bem_element_factory_factory
        (llg_mesh_pt->finite_element_pt(0));
    }

    /// Write out all args (in a parseable format) to a stream.
    virtual void dump_args(std::ostream& out_stream) const
    {
      MyCliArgs::dump_args(out_stream);

      out_stream
        << "mesh " << mesh_name << std::endl
        << "initial_m " << initial_m_name << std::endl
        << "h_app " << h_app_name << std::endl;
    }


    Mesh* llg_mesh_pt;
    Mesh* phi_1_mesh_pt;
    Mesh* phi_mesh_pt;
    InitialM::InitialMFctPt initial_m_fct_pt;
    HApp::HAppFctPt h_app_fct_pt;


    GenericPoissonProblem::FluxMeshFactoryFctPt phi_1_flux_mesh_factory_fct_pt;
    SemiImplicitFactories::BEMElementFactoryFctPt bem_element_factory_fct_pt;

    // Strings for input to factory functions
    std::string mesh_name;
    std::string initial_m_name;
    std::string h_app_name;
  };


} // End of oomph namespace

#endif
