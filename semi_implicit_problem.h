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

  template<class FIELD_ELEMENT, class MM_ELEMENT>
  class SemiImplicitHybridMicromagneticsProblem
  {
  public:

    // Function pointer for applied field.
    typedef typename MM_ELEMENT::TimeSpaceToDoubleVectFctPt AppliedFieldFctPt;

    // Function pointer for initial magnetisation.
    typedef void (*InitialMFctPt)(const double& t, const Vector<double> &x,
                                  Vector<double> &m);

    /// Constructor
    SemiImplicitHybridMicromagneticsProblem
    (Mesh* phi_1_mesh_pt,
     Mesh* phi_mesh_pt,
     Mesh* llg_mesh_pt,
     AppliedFieldFctPt applied_field_pt,
     Vector<std::pair<Vector<double>,double> >* corner_data_pt=0,
     bool pin_phi1=true)
    {
      // Set up phi_1 problem
      // ============================================================

      Phi_1_problem.set_bulk_mesh(phi_1_mesh_pt);
      for(unsigned b=0, nb=phi_1_mesh_pt->nboundary(); b < nb; b++)
        {
          // Phi_1 b.c.s are all Neumann but with flux determined elsewhere (by m)
          Phi_1_problem.set_neumann_boundary(b, 0);
        }

      if(pin_phi1)
        {
          // Pin a node which isn't involved in the boundary element method (we
          // have to pin something to avoid a singular Jacobian, can't be a
          // boundary node or things will go wrong with BEM).
          Node* pinned_phi_1_node_pt = phi_1_mesh_pt->get_non_boundary_node();
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
      Phi_1_problem.build();

      // Things will go wrong if the nodse of all meshes are not in
      // the same place:
#ifdef PARANOID
      if((phi_1_mesh_pt->nnode() != phi_mesh_pt->nnode())
         || (phi_1_mesh_pt->nnode() != llg_mesh_pt->nnode()))
        {
          std::ostringstream error_msg;
          error_msg << "Mesh nodes must be the same.";
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
#endif

      // Assign micromagnetics element pointers
      // ??ds dodgy...
      for(unsigned e=0, ne=phi_1_mesh_pt->nelement(); e < ne; e++)
        {
          FIELD_ELEMENT* ele_pt = dynamic_cast<FIELD_ELEMENT*>
            (phi_1_mesh_pt->element_pt(e));

          MM_ELEMENT* m_ele_pt = dynamic_cast<MM_ELEMENT*>
            (llg_mesh_pt->element_pt(e));

          ele_pt->set_micromag_element_pt( m_ele_pt);
        }


      // BEM handler:
      // ============================================================

      // Construct the BEM (must be done before pinning phi values)
      Bem_handler.set_bem_all_boundaries(phi_1_mesh_pt);
      // both zero because they are in seperate problems
      Bem_handler.input_index() = 0;
      Bem_handler.output_index() = 0;

      // Create an integration scheme
      // ??ds this is still a pretty bad way to do it I think...
      FIELD_ELEMENT* el_pt = dynamic_cast<FIELD_ELEMENT*>
        (phi_1_mesh_pt->element_pt(0));
      std::cout << el_pt->nodal_dimension() <<" " << el_pt->nvertex_node()  << std::endl;
      if(el_pt->nodal_dimension() == 2)
        {
          if(el_pt->nvertex_node() == 3)
            {
              Bem_handler.integration_scheme_pt() = new TVariableOrderGaussLegendre<1>;
            }
          else if(el_pt->nvertex_node() == 4)
            {
              Bem_handler.integration_scheme_pt() = new QVariableOrderGaussLegendre<1>;
            }
        }
      else if(el_pt->nodal_dimension() == 3)
        {
          if(el_pt->nvertex_node() == 4)
            {
              Bem_handler.integration_scheme_pt() = new TVariableOrderGaussLegendre<2>;
            }
          else if(el_pt->nvertex_node() == 8)
            {
              Bem_handler.integration_scheme_pt() = new QVariableOrderGaussLegendre<2>;
            }
          else
            {
              throw OomphLibError("Cannot determine element type.",
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
        }
      else
        {
          throw OomphLibError("Cannot determine element type.",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

      Bem_handler.input_corner_data_pt() = corner_data_pt;

      Bem_handler.build();


      // Now we can set up phi problem
      // ============================================================

      Phi_problem.set_bulk_mesh(phi_mesh_pt);
      unsigned nboundary = phi_mesh_pt->nboundary();
      Phi_boundary_values_pts.assign(nboundary, 0);
      for(unsigned b=0; b < nboundary; b++)
        {
          // Phi is determined by BEM
          LinearAlgebraDistribution* dist_pt =
            new LinearAlgebraDistribution(0, phi_mesh_pt->nboundary_node(b), false);

          Phi_boundary_values_pts[b] = new DoubleVector(dist_pt);
          Phi_problem.set_dirichlet_boundary_by_vector(b, Phi_boundary_values_pts[b]);
        }
      Phi_problem.build();

      // Assign micromagnetics element pointers
      // ??ds dodgy...
      for(unsigned e=0, ne=phi_mesh_pt->nelement(); e < ne; e++)
        {
          FIELD_ELEMENT* ele_pt = dynamic_cast<FIELD_ELEMENT*>
            (phi_mesh_pt->element_pt(e));

          MM_ELEMENT* m_ele_pt = dynamic_cast<MM_ELEMENT*>
            (llg_mesh_pt->element_pt(e));

          ele_pt->set_micromag_element_pt(m_ele_pt);
        }

      // LLG problem:
      // ============================================================

      llg_sub_problem_pt()->bulk_mesh_pt() = llg_mesh_pt;

      //??ds while we still have phi in MM elements pin them all
      for(unsigned nd=0, nnode=llg_mesh_pt->nnode(); nd<nnode; nd++)
        {
          Node* nd_pt = llg_mesh_pt->node_pt(nd);

          unsigned phi_index = llg_element_pt()->phi_index_micromag();
          unsigned phi_1_index = llg_element_pt()->phi_1_index_micromag();

          nd_pt->pin(phi_index);
          nd_pt->pin(phi_1_index);
        }

      // Get timestepper from mesh
      TimeStepper* ts_pt = llg_mesh_pt->node_pt(0)->time_stepper_pt();
      llg_sub_problem_pt()->add_time_stepper_pt(ts_pt);

      // Magnetic parameters
      llg_sub_problem_pt()->mag_parameters_pt()->set_nmag_rectangle();
      llg_sub_problem_pt()->applied_field_fct_pt() = applied_field_pt;

      llg_sub_problem_pt()->build();

      // Assign phi element pointers
      // ??ds dodgy...
      for(unsigned e=0, ne=llg_mesh_pt->nelement(); e < ne; e++)
        {
          MM_ELEMENT* ele_pt = dynamic_cast<MM_ELEMENT*>
            (llg_mesh_pt->element_pt(e));
          ele_pt->magnetostatic_field_element_pt() =
            dynamic_cast<FIELD_ELEMENT*>(phi_mesh_pt->element_pt(e));
        }


      // Linear solvers
      // ============================================================

      // CG + AMG preconditioner is excellent for poisson solves:
      IterativeLinearSolver* CG_pt = new CG<CRDoubleMatrix>;
      Phi_problem.linear_solver_pt() = CG_pt;
#ifdef OOMPH_HAS_HYPRE
      HyprePreconditioner* AMG_pt = new HyprePreconditioner;
      AMG_pt->hypre_method() = HyprePreconditioner::BoomerAMG;
      CG_pt->preconditioner_pt() = AMG_pt;
#endif

      // Same for phi_1
      IterativeLinearSolver* CG_1_pt = new CG<CRDoubleMatrix>;
      Phi_1_problem.linear_solver_pt() = CG_1_pt;
#ifdef OOMPH_HAS_HYPRE
      HyprePreconditioner* AMG_1_pt = new HyprePreconditioner;
      AMG_1_pt->hypre_method() = HyprePreconditioner::BoomerAMG;
      CG_1_pt->preconditioner_pt() = AMG_1_pt;
#endif

      // ??ds No idea yet for LLG solver.. something fancy?

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
      Phi_1_problem.newton_solve();

      // update boundary values of phi
      std::cout << "solving BEM" << std::endl;
      Bem_handler.get_bem_values(Phi_boundary_values_pts);

      // solve for phi
      std::cout << "solving phi" << std::endl;
      Phi_problem.newton_solve();

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
    void set_initial_condition(const InitialMFctPt initial_m_pt);

    /// Initialise timestep: only llg problem has a timestep.
    void initialise_dt(const double &dt)
    {llg_sub_problem_pt()->initialise_dt(dt);}

    /// Output
    void doc_solution(DocInfo &doc_info, const unsigned &npts=2) const;

    void average_magnetostatic_field(Vector<double> &average_magnetostatic_field) const;

    /// Access functions
    // =================================================================

    /// Get pointer to an LLG element for looking up info. All elements in
    /// mesh should be the same type otherwise preconditioning framework
    /// will fail so this should be safe.
    MM_ELEMENT* llg_element_pt() const
    {return dynamic_cast<MM_ELEMENT*>
        (llg_sub_problem_pt()->bulk_mesh_pt()->element_pt(0));}

    /// Get access to magnetic parameters - only relevant in LLG problem so
    /// return the pointer from that problem.
    MagneticParameters* mag_parameters_pt()
    {return llg_sub_problem_pt()->mag_parameters_pt();}

    /// \short Const access function for LLG_problem.
    const ImplicitLLGProblem<MM_ELEMENT>* llg_sub_problem_pt() const
    {return &LLG_problem;}

    /// \short Non-const acess function for LLG_problem.
    ImplicitLLGProblem<MM_ELEMENT>* llg_sub_problem_pt() {return &LLG_problem;}

    const Mesh* bem_mesh_pt() const {return Bem_handler.bem_mesh_pt();}

    const DenseMatrix<double>* bem_matrix_pt() const {return Bem_handler.bem_matrix_pt();}


    GenericPoissonProblem<FIELD_ELEMENT>* phi_problem_pt()
    {return &Phi_problem;}

    // /// Set the list of sharp corners in the mesh to be a rectangle.
    // void set_rectangular_corners()
    // {Bem_handler.Mesh_angles_type = "rectangular";}

  private:

    /// Object to provide all BEM related capabilities.
    BoundaryElementHandler<MicromagFaceElement<FIELD_ELEMENT> > Bem_handler;

    /// Problem to solve for phi_1 (the BEM pre-calculation).
    GenericPoissonProblem<FIELD_ELEMENT,
                          MagnetostaticFieldFluxElement<FIELD_ELEMENT> >
    Phi_1_problem;

    /// Problem to solve for phi (the magnetostatic potential).
    GenericPoissonProblem<FIELD_ELEMENT> Phi_problem;

    /// Problem to solve for the magnetisation change.
    ImplicitLLGProblem<MM_ELEMENT> LLG_problem;

    /// Intermediate storage for results of bem (ideally we would have it
    /// call a function to get the boundary values filled in but c++ member
    /// functions pointers are useless...)
    Vector<DoubleVector*> Phi_boundary_values_pts;
  };


  template<class FIELD_ELEMENT, class MM_ELEMENT>
  void SemiImplicitHybridMicromagneticsProblem<FIELD_ELEMENT, MM_ELEMENT>::
  set_initial_condition(const InitialMFctPt initial_m_pt)
  {
    // Backup time in global Time object
    double backed_up_time=llg_sub_problem_pt()->time_pt()->time();

    // Past history needs to be established for t=time0-deltat, ...
    // Then provide current values (at t=time0) which will also form
    // the initial guess for the first solve at t=time0+deltat

    // Get M indicies
    Vector<unsigned> m_index_micromag(3,0);
    MM_ELEMENT* elem_pt = dynamic_cast<MM_ELEMENT* >(llg_sub_problem_pt()->mesh_pt()->element_pt(0));
    for(unsigned i=0; i<3; i++)
      {
        m_index_micromag[i] = elem_pt->m_index_micromag(i);
      }

    // Find number of nodes in mesh
    unsigned num_nod = llg_sub_problem_pt()->mesh_pt()->nnode();

    // Set continuous times at previous timesteps:
    int nprev_steps=llg_sub_problem_pt()->time_stepper_pt()->nprev_values();
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
            Vector<double> m(3,0.0), x(dim,0.0);
            llg_sub_problem_pt()->mesh_pt()->node_pt(n)->position(t,x);
            initial_m_pt(time,x,m);

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

  template<class FIELD_ELEMENT, class MM_ELEMENT>
  void SemiImplicitHybridMicromagneticsProblem<FIELD_ELEMENT, MM_ELEMENT>::
  doc_solution(DocInfo &doc_info, const unsigned &npts) const
  {
    using namespace StringConversion;

    // Output llg solution
    std::ofstream soln_file((doc_info.directory() + "/soln"
                             + doc_info.number_as_string() + ".dat").c_str());
    llg_sub_problem_pt()->mesh_pt()->output(soln_file,npts);
    soln_file.close();

    // Output the magnetostatic field data
    std::ofstream field_file((doc_info.directory() + "/field"
                              + doc_info.number_as_string() + ".dat").c_str());
    for(unsigned e=0, ne=Phi_problem.mesh_pt()->nelement(); e < ne; e++)
      {
        FIELD_ELEMENT* ele_pt = dynamic_cast<FIELD_ELEMENT*>
          (Phi_problem.mesh_pt()->element_pt(e));
        ele_pt->output(field_file,npts);
      }
    field_file.close();

    std::ofstream phi1_file((doc_info.directory() + "/phione"
                             + doc_info.number_as_string() + ".dat").c_str());
    Phi_1_problem.mesh_pt()->output(phi1_file,npts);
    phi1_file.close();

    // Write average magnetisations to a file
    Vector<double> m;
    llg_sub_problem_pt()->mean_magnetisation(m);

    std::ofstream avgs((doc_info.directory() +"/averages").c_str(),
                       std::ios::app);
    avgs << llg_sub_problem_pt()->time();
    for(unsigned j=0; j<3; j++) avgs << " " << m[j];
    avgs << std::endl;
    avgs.close();


    // Write average field to a file
    Vector<double> hms;
    average_magnetostatic_field(hms);
    std::ofstream field_avgs((doc_info.directory() +"/field_averages").c_str(),
                             std::ios::app);
    field_avgs << llg_sub_problem_pt()->time();
    for(unsigned j=0; j<3; j++) field_avgs << " " << hms[j];
    field_avgs << std::endl;
    field_avgs.close();


    // Get average (and standard deviation) of |m| - 1 and |m|.dm/dn
    double m_error_avg(0), m_error_stddev(0), orthogonality_error_avg(0),
      orthogonality_error_stddev(0);
    llg_sub_problem_pt()->norm_m_error(m_error_avg, m_error_stddev);
    llg_sub_problem_pt()->orthogonality_m_error
      (orthogonality_error_avg, orthogonality_error_stddev);

    // Write them to file
    std::ofstream errors((doc_info.directory()+"/errors").c_str(),std::ios::app);
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
    doc_info.number()++;
  }

  //============================================================
  //
  //============================================================
  template<class FIELD_ELEMENT, class MM_ELEMENT>
  void SemiImplicitHybridMicromagneticsProblem<FIELD_ELEMENT, MM_ELEMENT>::
  average_magnetostatic_field(Vector<double> &average_magnetostatic_field) const
  {
    const unsigned nodal_dim = checked_dynamic_cast<FIELD_ELEMENT*>
      (Phi_problem.mesh_pt()->element_pt(0))->node_pt(0)->ndim();

    // Pick a point in the middle of the element
    const Vector<double> s(nodal_dim, 0.3);
    Vector<double> total_dphidx(nodal_dim,0.0);

    // Loop over all elements calculating the value in the middle of the element
    for(unsigned e=0, ne=Phi_problem.mesh_pt()->nelement(); e < ne; e++)
      {
        FIELD_ELEMENT* ele_pt = checked_dynamic_cast<FIELD_ELEMENT*>
          (Phi_problem.mesh_pt()->element_pt(e));

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
    double nele = double(Phi_problem.mesh_pt()->nelement());
    average_magnetostatic_field.assign(3,0.0);
    for(unsigned j=0; j<nodal_dim; j++)
      {
        average_magnetostatic_field[j] = - total_dphidx[j] / nele;
      }
  }

  // template<class FIELD_ELEMENT, class MM_ELEMENT>
  // void SemiImplicitHybridMicromagneticsProblem<FIELD_ELEMENT, MM_ELEMENT>::
  // doc_magnetostatic_field(DocInfo &doc_info) const
  // {
  //   // Number of plot points
  //   unsigned npts;
  //   npts=2;

  //   std::ofstream soln_file((doc_info.directory() + "/soln"
  //                            + doc_info.number_as_string() + ".dat").c_str());

  //   Phi_problem.doc_solution(soln

  //   soln_file.close();
  // }

} // End of oomph namespace

#endif
