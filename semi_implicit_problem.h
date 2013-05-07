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


    void build()
    {

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


      // Build the LLG part of the problem
      ImplicitLLGProblem::build();

//       // Linear solvers
//       // ============================================================

//       // CG + AMG preconditioner is excellent for poisson solves:
//       IterativeLinearSolver* CG_pt = new CG<CRDoubleMatrix>;
//       phi_problem_pt()->linear_solver_pt() = CG_pt;
// #ifdef OOMPH_HAS_HYPRE
//       HyprePreconditioner* AMG_pt = new HyprePreconditioner;
//       AMG_pt->hypre_method() = HyprePreconditioner::BoomerAMG;
//       CG_pt->preconditioner_pt() = AMG_pt;
// #endif

//       // Same for phi_1
//       IterativeLinearSolver* CG_1_pt = new CG<CRDoubleMatrix>;
//       Phi_1_problem_pt->linear_solver_pt() = CG_1_pt;
// #ifdef OOMPH_HAS_HYPRE
//       HyprePreconditioner* AMG_1_pt = new HyprePreconditioner;
//       AMG_1_pt->hypre_method() = HyprePreconditioner::BoomerAMG;
//       CG_1_pt->preconditioner_pt() = AMG_1_pt;
// #endif

//       // ??ds No idea yet for LLG solver.. something fancy?

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

    BoundaryElementHandlerBase* &bem_handler_pt() {return Bem_handler_pt;}


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
    BoundaryElementHandlerBase* Bem_handler_pt;

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

} // End of oomph namespace

#endif
