#ifndef OOMPH_IMPLICIT_LLG_EXCHANGE_PROBLEM_H
#define OOMPH_IMPLICIT_LLG_EXCHANGE_PROBLEM_H

#include "generic.h"
#include "./boundary_element_handler.h"
#include "./micromagnetics_element.h"
#include "./micromagnetics_element.cc" //??ds shouldn't need...
#include "./vector_helpers.h"

// Default mesh
#include "meshes/simple_rectangular_tri_mesh.h"


namespace oomph
{

  // ============================================================
  ///
  // ============================================================
  template<class ELEMENT>
  class ImplicitLLGProblem : public Problem
  {
  public:

    // Function pointer for initial magnetisation.
    typedef void (*InitialMFctPt)(const double& t, const Vector<double> &x,
                                  Vector<double> &m);

    // Function pointer for applied field.
    typedef typename ELEMENT::TimeSpaceToDoubleVectFctPt AppliedFieldFctPt;

    /// Default constructor - do nothing except nulling pointers.
    ImplicitLLGProblem() :
      Use_mid_point_method(false), Ele_pt(0), Applied_field_fct_pt(0),
      Renormalise_each_time_step(false)
    {}

    /// Make a problem with the default mesh and timestepper
    ImplicitLLGProblem(const unsigned &nx,
                       const unsigned &ny,
                       const double &lx,
                       const double &ly,
                       AppliedFieldFctPt input_applied_field_pt,
                       const bool adaptive_flag=false,
                       const bool use_mid_point_method=false) :
      Ele_pt(0), Applied_field_fct_pt(input_applied_field_pt),
      Use_mid_point_method(use_mid_point_method)
    {
      // Create timestepper
      TimeStepper* time_stepper_pt = 0;
      if(Use_mid_point_method)
        {
          time_stepper_pt = new BDF<1>(adaptive_flag);
        }
      else
        {
          time_stepper_pt = new BDF<2>(adaptive_flag);
        }
      add_time_stepper_pt(time_stepper_pt);

      // Create mesh
      bulk_mesh_pt() = new SimpleRectangularTriMesh<ELEMENT>(nx,ny,lx,ly,time_stepper_pt);
      bulk_mesh_pt()->setup_boundary_element_info();

      // Set up some simple parameters
      mag_parameters_pt()->set_simple_llg_parameters();

      // Build the problem
      build();
    }

    // /// Alternative constructor with any mesh and timestepper.
    // // (we have to ask for timestepper_pt because it must have already been
    // // set in the mesh but it also needs to be set in the problem).
    // ImplicitLLGProblem(Mesh* input_mesh_pt, TimeStepper* time_stepper_pt,
    //                    AppliedFieldFctPt input_applied_field_pt) : //??ds should be able to passi
    //   Ele_pt(0), Applied_field_fct_pt(input_applied_field_pt)
    // {
    //   // Set mesh, timestepper and applied field
    //   bulk_mesh_pt() = input_mesh_pt;
    //   add_time_stepper_pt(time_stepper_pt);

    //   // Everything else
    //   build();
    // }

    /// Function that does the real work of the constructors.
    void build()
    {
#ifdef PARANOID
      if((bulk_mesh_pt() == 0) || (applied_field_fct_pt() == 0)
         || (time_stepper_pt() == 0))
        {
          std::ostringstream error_msg;
          error_msg << "Must assign mesh, timestepper and applied "
                    << "field pointers must be set "
                    << "before calling build().";
          throw OomphLibError(error_msg.str(),
                              "ImplicitLLGProblem::build()",
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      // Write out parameters data.
      Magnetic_parameters.output(std::cout);

      // Fiddle with newton paramters if needed
      // max_residuals() *= 2;

      // Get a pointer to an element for use later
      Ele_pt = dynamic_cast<ELEMENT*>(bulk_mesh_pt()->element_pt(0));

      // Find out the problem dimensions
      Dim = Ele_pt->nodal_dimension();

      // Boundary conditions - our finite element discretisation requires
      // residual addition of m x dm/dn along the boundary (due to need to
      // reduce the order of the laplacian on m in exchange field).

      // // Create mesh of MicromagFluxElement<ELEMENT> to add boundary
      // // contribution (if any).
      // surface_exchange_mesh_pt() = new Mesh;
      // for(unsigned b=0, nb=bulk_mesh_pt()->nboundary(); b < nb; b++)
      //   {
      //     create_surface_exchange_elements(b);
      //   }

      // Finish off elements
      for(unsigned i=0; i< bulk_mesh_pt()->nelement(); i++)
        {
          ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(bulk_mesh_pt()->element_pt(i));

          // Set pointer to continuous time
          elem_pt->time_pt() = time_pt();

          // Set values for magnetic parameters
          elem_pt->magnetic_parameters_pt() = mag_parameters_pt();

          // Set pointer for an applied field
          elem_pt->applied_field_pt() = applied_field_fct_pt();
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
      add_sub_mesh(bulk_mesh_pt());
      // add_sub_mesh(surface_exchange_mesh_pt());
      build_global_mesh();

      // Do equation numbering
      std::cout << "LLG Number of equations: " << assign_eqn_numbers() << std::endl;
    }

    /// Destructor
    ~ImplicitLLGProblem()
    {
      // mesh is cleaned up by problem base class
      // timestepper is cleaned up by problem base class
    }

    /// Renormalise magnetisation to 1 (needed with BDF2)
    void renormalise_magnetisation()
    {
      for(unsigned nd=0; nd<bulk_mesh_pt()->nnode(); nd++)
        {
          Node* nd_pt = bulk_mesh_pt()->node_pt(nd);

          // Get m vector
          Vector<double> m_values(3,0.0);
          for(unsigned j=0; j<3; j++) m_values[j] = nd_pt->value(m_index(j));

          // Normalise
          VectorOps::normalise(m_values);

          // Write m vector
          for(unsigned j=0; j<3; j++) nd_pt->set_value(m_index(j),m_values[j]);
        }
    }

    // If we are using midpoint then apply the required update (Malidi2005)
    void midpoint_update()
    {
      // Get a pointer
      ELEMENT* bulk_elem_pt =
        checked_dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(0));

      // Get m indicies
      Vector<unsigned> m_indices(3,0);
      for(unsigned j=0; j<3; j++)
        {
          m_indices[j] = bulk_elem_pt->m_index_micromag(j);
        }

      // Apply the update
      for(unsigned i_nd=0; i_nd<mesh_pt()->nnode(); i_nd++)
        {
          Node* nd_pt = mesh_pt()->node_pt(i_nd);
          for(unsigned j=0; j<3; j++)
            {
              // m[n] --> 2*m[n] - m[n-1]
              // because we have finished this timestep and moved forward one
              double new_m = 2*(nd_pt->value(m_indices[j]))
                - (nd_pt->value(1, m_indices[j]));

              nd_pt->set_value(m_indices[j], new_m);
            }
        }

      // Finally push the time forwards another half step. Note that
      // time_pt()->dt() is half a step because out previous step was
      // real_dt/2.
      time() += time_pt()->dt();
    }


    void actions_after_newton_solve()
    {
      // If we're using BDF we need to keep M normalised.
      if(renormalise_each_time_step())
        {
          renormalise_magnetisation();
        }

      // If we're usign midpoint method then do the update needed to
      // convert from BDF1 to a midpoint timestep.
      if(Use_mid_point_method)
        {
          midpoint_update();
        }
    }


    /// Output solution
    void doc_solution(DocInfo &doc_info) const
    {
      std::ofstream some_file;


      // Number of plot points
      unsigned npts = 2;

      // Output solution with specified number of plot points per element
      std::string filename(doc_info.directory() + "/soln" +
                           doc_info.number_as_string() + ".dat");
      some_file.open(filename.c_str());
      bulk_mesh_pt()->output(some_file,npts);
      some_file.close();


      // Get average (and standard deviation) of |m| - 1 and |m|.dm/dn
      double m_error_avg(0), m_error_stddev(0);
      norm_m_error(m_error_avg, m_error_stddev);

      // Write them to file
      std::ofstream errors((doc_info.directory()+"/errors").c_str(), std::ios::app);
      errors << time()
             << " " << m_error_avg
             << " " << m_error_stddev
             << std::endl;
      errors.close();


      // Increment number used as output label
      doc_info.number()++;
    }

    /// Set up an initial M
    void set_initial_condition(const InitialMFctPt initial_m_pt);

    /// Error for adaptive timestepper (rms of nodal error determined by
    /// comparison with explicit timestepper result).
    double global_temporal_error_norm()
    {
      double global_error = 0.0;

      //Find out how many nodes there are in the problem
      unsigned n_node = bulk_mesh_pt()->nnode();

      //Loop over the nodes and calculate the estimated error in the values
      for(unsigned i=0;i<n_node;i++)
        {
          for(unsigned j=0; j<3; j++)
            {
              // Get error in solution: Difference between predicted
              // (explicit timestepper) and actual value for each m.
              double error = bulk_mesh_pt()->node_pt(i)->time_stepper_pt()->
                temporal_error_in_value(bulk_mesh_pt()->node_pt(i),m_index(j));

              //Add the square of the individual error to the global error
              global_error += error*error;
            }
        }

      // Divide by the number of nodes
      global_error /= 3*double(n_node);

      // Return square root...
      return sqrt(global_error);
    }

    // Lots of magnetisation manipulation functions
    // ============================================================

    // Loop over all nodes in bulk mesh and get magnetisations
    void get_nodal_magnetisations(Vector< Vector<double> > &m_list) const
    {
      unsigned nnode = bulk_mesh_pt()->nnode();
      m_list.resize(nnode);
      for(unsigned nd=0; nd<nnode; nd++)
        {
          m_list[nd].assign(3,0.0);
          for(unsigned j=0; j<3; j++)
            {m_list[nd][j] = bulk_mesh_pt()->node_pt(nd)->value(m_index(j));}
        }
    }

    void get_nodal_two_norms(Vector<double> &output) const
    {
      Vector< Vector<double> > m; get_nodal_magnetisations(m);

      output.assign(m.size(), 0.0);
      transform(m.begin(), m.end(), output.begin(), VectorOps::two_norm);
    }

    double mean_nodal_magnetisation_length() const
    {
      Vector<double> ms; get_nodal_two_norms(ms);
      return VectorOps::mean(ms);
    }

    double mean_mx() const
    {
      Vector< Vector<double> > ms; get_nodal_magnetisations(ms);

      Vector<double> mx(ms.size(),0.0);
      for(unsigned i=0; i<ms.size(); i++)
        {mx[i] = ms[i][0];}

      return VectorOps::mean(mx);
    }

    double mean_mz() const
    {
      Vector< Vector<double> > ms; get_nodal_magnetisations(ms);

      Vector<double> mz(ms.size(),0.0);
      for(unsigned i=0; i<ms.size(); i++)
        {mz[i] = ms[i][2];}

      return VectorOps::mean(mz);
    }

    void mean_magnetisation(Vector<double> &m) const
    {
      // Reset m
      m.assign(3,0.0);

      // Sum over all nodes in mesh
      unsigned nnode = bulk_mesh_pt()->nnode();
      for(unsigned nd=0; nd<nnode; nd++)
        {
          for(unsigned j=0; j<3; j++)
            {
              m[j] += bulk_mesh_pt()->node_pt(nd)->value(m_index(j));
            }
        }

      // Divide by number of nodes in mesh
      for(unsigned j=0; j<3; j++)
        {
          m[j] /= nnode;
        }
    }

    /// \short Abs of mean difference of actual m and m given by a function
    /// at the middle of each element.
    double compare_m_with_function(const InitialMFctPt fct_pt) const
    {
      double diff = 0.0;

      // Compare at middle of element
      Vector<double> s(3,0.0);
      for(unsigned j=0; j<dim(); j++) s[j] = 0.5;

      // Sum the difference over all elements
      for(unsigned e=0, ne=bulk_mesh_pt()->nelement(); e < ne; e++)
        {
          ELEMENT* ele_pt = dynamic_cast<ELEMENT* >
            (bulk_mesh_pt()->element_pt(e));

          Vector<double> numerical_m(3,0.0);
          ele_pt->interpolated_m_micromag(s,numerical_m);

          Vector<double> x(dim(),0.0), exact_m(3,0.0);
          ele_pt->interpolated_x(s,x);
          double t = this->time();
          fct_pt(t, x, exact_m);

          for(unsigned j=0; j<3; j++)
            {
              diff += std::abs(numerical_m[j] - exact_m[j]);
            }
        }

      // Divide to get the mean
      diff /= (3 * bulk_mesh_pt()->nelement());

      return diff;
    }

    double mean_norm_m_error() const
    {
      double temp_error = 0;

      unsigned nnode = bulk_mesh_pt()->nnode();
      for(unsigned nd=0; nd<nnode; nd++)
        {
          Node* nd_pt = bulk_mesh_pt()->node_pt(nd);
          Vector<double> m_values(3,0.0);
          for(unsigned j=0; j<3; j++) m_values[j] = nd_pt->value(m_index(j));

          double m_2norm = VectorOps::two_norm(m_values);
          temp_error += m_2norm - 1;
        }

      return temp_error/double(nnode);
    }

    void norm_m_error(double &m_error_avg, double &m_error_stddev) const
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
          temp_stddev += std::pow( m_2norm - 1 - m_error_avg ,2);
        }

      temp_stddev /= nnode;
      m_error_stddev = std::sqrt(temp_stddev);
    }


    /// Elementwise calculation of m . dm/dn on the boundaries. This gives a good measure of
    /// orthogonality of m and dmdn.
    double mean_orthogonality_m_error() const
    {
      unsigned nele = surface_exchange_mesh_pt()->nelement();
      double temp_sum = 0.0;

      // Calculate m . dm/dn for each element
      for(unsigned e=0; e < nele; e++)
        {
          MicromagFluxElement<ELEMENT>* ele_pt
            = checked_dynamic_cast<MicromagFluxElement<ELEMENT>*>
            (surface_exchange_mesh_pt()->element_pt(e));
          Vector<double> s(ele_pt->dim(), 0.5); // middle of element
          temp_sum += ele_pt->interpolated_mdotdmdn_micromag(s);
        }

      // Take average
      return temp_sum / double(nele);
    }

    void orthogonality_m_error(double &orthogonality_error_avg,
                               double &orthogonality_error_stddev) const
    {
      orthogonality_error_avg = mean_orthogonality_m_error();

      unsigned nele = surface_exchange_mesh_pt()->nelement();
      double temp_sum = 0.0;

      // Calculate m . dm/dn for each element
      for(unsigned e=0; e < nele; e++)
        {
          MicromagFluxElement<ELEMENT>* ele_pt
            = checked_dynamic_cast<MicromagFluxElement<ELEMENT>*>
            (surface_exchange_mesh_pt()->element_pt(e));
          Vector<double> s(ele_pt->dim(), 0.5); // middle of element
          temp_sum += pow( ele_pt->interpolated_mdotdmdn_micromag(s)
                           - orthogonality_error_avg, 2);
        }

      // Take stddev
      orthogonality_error_stddev = std::sqrt(temp_sum / double(nele));
    }

    // Access functions
    // ============================================================

    unsigned m_index(const unsigned &j) const
    {return Ele_pt->m_index_micromag(j);}

    unsigned phi_index() const
    {return Ele_pt->phi_index_micromag();}

    unsigned phi_1_index() const
    {return Ele_pt->phi_1_index_micromag();}

    /// \short Get problem dimension (nodal dimension).
    const unsigned dim() const
    {return Dim;}

    /// \short Non-const access function for Applied_field_fct_pt.
    AppliedFieldFctPt& applied_field_fct_pt() {return Applied_field_fct_pt;}

    /// \short Const access function for Applied_field_fct_pt.
    AppliedFieldFctPt applied_field_fct_pt() const {return Applied_field_fct_pt;}

    /// \short Non-const access function for Magnetic_parameters_pt.
    MagneticParameters* mag_parameters_pt()
    {return &Magnetic_parameters;}

    /// \short Const access function for Magnetic_parameters_pt.
    const MagneticParameters* mag_parameters_pt() const
    {return &Magnetic_parameters;}

    /// \short Non-const access function for Renormalise_each_time_step.
    bool& renormalise_each_time_step() {return Renormalise_each_time_step;}

    /// \short Const access function for Renormalise_each_time_step.
    bool renormalise_each_time_step() const {return Renormalise_each_time_step;}

    /// \short Non-const access function for Surface_exchange_mesh_pt.
    Mesh*& surface_exchange_mesh_pt() {return Surface_exchange_mesh_pt;}

    /// \short Const access function for Surface_exchange_mesh_pt.
    Mesh* surface_exchange_mesh_pt() const
    {
#ifdef PARANOID
      if(Surface_exchange_mesh_pt == 0)
        {
          std::ostringstream error_msg;
          error_msg << "Surface exchange mesh pointer not set.";
          throw OomphLibError(error_msg.str(),
                              "",
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      return Surface_exchange_mesh_pt;
    }

    /// \short Non-const access function for Bulk_mesh_pt.
    Mesh*& bulk_mesh_pt() {return Bulk_mesh_pt;}

    /// \short Const access function for Bulk_mesh_pt.
    Mesh* bulk_mesh_pt() const {return Bulk_mesh_pt;}

    /// Use mid-point method (rather than bdf2)? Public to avoid lots of
    /// access function junk, see how it goes...
    bool Use_mid_point_method;

  private:

    /// A pointer to an element (used to get index numbers)
    const ELEMENT* Ele_pt;

    /// The problem dimension
    unsigned Dim;

    /// Pointer to the applied field.
    AppliedFieldFctPt Applied_field_fct_pt;

    /// Magnetic parameters storage. ??ds should maybe go in meshes?
    MagneticParameters Magnetic_parameters;

    /// Normalise magnetisation problem after each step?
    bool Renormalise_each_time_step;

    /// Pointer to bulk mesh (i.e. the magnetic material).
    Mesh* Bulk_mesh_pt;

    /// Pointer to the mesh of face elements for dealing with boundary
    /// conditions caused by reducing order of differentiation in exchange
    /// term.
    Mesh* Surface_exchange_mesh_pt;

    /// Create
    void create_surface_exchange_elements(const unsigned& b);

    /// Inaccessible copy constructor
    ImplicitLLGProblem(const ImplicitLLGProblem & dummy)
    {BrokenCopy::broken_copy("ImplicitLLGProblem");}

    /// Inaccessible assignment operator
    void operator=(const ImplicitLLGProblem &dummy)
    {BrokenCopy::broken_assign("ImplicitLLGProblem");}

  };

  //======================================================================
  /// Set up the initial conditions
  //======================================================================
  template<class ELEMENT>
  void ImplicitLLGProblem<ELEMENT>::
  set_initial_condition(const InitialMFctPt initial_m_pt)
  {
    // Backup time in global Time object
    double backed_up_time=this->time_pt()->time();

    // Past history needs to be established for t=time0-deltat, ...
    // Then provide current values (at t=time0) which will also form
    // the initial guess for the first solve at t=time0+deltat

    // Get M indicies
    Vector<unsigned> m_index_micromag(3,0);
    ELEMENT* elem_pt = dynamic_cast<ELEMENT* >(this->mesh_pt()->element_pt(0));
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
            unsigned dim = mesh_pt()->node_pt(n)->ndim();

            // Get initial value of m from inputs
            Vector<double> m(3,0.0), x(dim,0.0);
            this->mesh_pt()->node_pt(n)->position(t,x);
            initial_m_pt(time,x,m);

            // Set initial condition on m
            for(unsigned i=0; i<3; i++)
              this->mesh_pt()->node_pt(n)->set_value(t,m_index_micromag[i],m[i]);
          }
      }

    // Reset backed up time for global timestepper
    this->time_pt()->time()=backed_up_time;
  }


  //======================================================================
  /// Create
  //======================================================================
  template<class ELEMENT>
  void ImplicitLLGProblem<ELEMENT>::
  create_surface_exchange_elements(const unsigned& b)
  {
    // Loop over the bulk elements adjacent to boundary b
    for(unsigned e=0, ne=bulk_mesh_pt()->nboundary_element(b);e<ne;e++)
      {
        // Get pointer to the bulk element that is adjacent to boundary b
        ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>
          (bulk_mesh_pt()->boundary_element_pt(b,e));

        // What is the index of the face of the bulk element at the boundary
        int face_index = bulk_mesh_pt()->face_index_at_boundary(b,e);

        // Build the corresponding prescribed-flux element
        MicromagFluxElement<ELEMENT>* flux_element_pt =
          new MicromagFluxElement<ELEMENT>(bulk_elem_pt,face_index);

        // Pass a pointer to the flux element to the bulk element
        bulk_elem_pt->add_face_element_pt(flux_element_pt);

        // Add the prescribed-flux element to the mesh
        surface_exchange_mesh_pt()->add_element_pt(flux_element_pt);

      } // End of loop over bulk elements adjacent to boundary b
  }


}

#endif
