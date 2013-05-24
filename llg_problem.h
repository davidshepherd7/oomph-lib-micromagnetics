#ifndef OOMPH_IMPLICIT_LLG_EXCHANGE_PROBLEM_H
#define OOMPH_IMPLICIT_LLG_EXCHANGE_PROBLEM_H

#include "generic.h"
#include "./boundary_element_handler.h"
#include "./micromagnetics_element.h"
#include "./micromagnetics_element.cc" //??ds shouldn't need...

#include "./vector_helpers.h"
#include "./magnetics_helpers.h"
#include "./my_generic_problem.h"

#include "./mallinson_solution.h"

#include <algorithm>


namespace oomph
{

  // ============================================================
  ///
  // ============================================================
  class LLGProblem : public MyProblem
  {
  public:

    /// Default constructor - do nothing except nulling pointers.
    LLGProblem() :
      Compare_with_mallinson(false),
      Swap_solver_large_dt(false),
      Applied_field_fct_pt(0),
      Renormalise_each_time_step(false)
    {}

    /// Function that does the real work of the constructors.
    void build()
    {
#ifdef PARANOID
      if((bulk_mesh_pt() == 0) || (applied_field_fct_pt() == 0)
         || (this->time_stepper_pt() == 0))
        {
          std::ostringstream error_msg;
          error_msg << "Must assign mesh, timestepper and applied "
                    << "field pointers to non-null values "
                    << "before calling build().\n"
                    << "bulk_mesh_pt() = " << bulk_mesh_pt() << "\n"
                    << "applied_field_fct_pt() = " << applied_field_fct_pt() << "\n"
                    << "this->time_stepper_pt() = " << this->time_stepper_pt()
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

      // Finish off elements
      for(unsigned i=0; i< bulk_mesh_pt()->nelement(); i++)
        {
          MicromagEquations* elem_pt = checked_dynamic_cast<MicromagEquations*>
            (bulk_mesh_pt()->element_pt(i));

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
      this->add_sub_mesh(bulk_mesh_pt());
      // add_sub_mesh(surface_exchange_mesh_pt());
      this->build_global_mesh();


      // ??ds For if we want to swap solver for large dt. For now swap if
      // we are using any iterative solver.
      if(dynamic_cast<IterativeLinearSolver*>(linear_solver_pt()) != 0)
        {
          Swap_solver_large_dt = true;
        }
      My_linear_solver_pt = linear_solver_pt();

      // Do equation numbering
      oomph_info << "LLG Number of equations: " << this->assign_eqn_numbers() << std::endl;
      oomph_info << "Number of sub meshes: " << this->nsub_mesh() << std::endl;
    }

    /// Destructor
    virtual ~LLGProblem()
    {
      // mesh is cleaned up by problem base class
      // timestepper is cleaned up by problem base class
      delete Magnetic_parameters_pt; Magnetic_parameters_pt = 0;
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

    void actions_before_newton_step()
    {
      std::cout << std::endl
                << "Newton step " << Nnewton_iter_taken + 1 << std::endl
                << "---------------------------------------" << std::endl;

      if(Swap_solver_large_dt)
        {
          if(this->time_pt()->dt() > 1e-2)
            {
              linear_solver_pt() = Default_linear_solver_pt;
            }
          else
            {
              linear_solver_pt() = My_linear_solver_pt;
            }
        }
    }

    void actions_before_newton_solve()
    {
      // Call lower level actions function
      MyProblem::actions_before_newton_solve();
    }

    void actions_after_newton_solve()
    {

      std::cout << std::endl
                << "Finalising" << std::endl
                << "-------------------------------------------" << std::endl;

      // If we're using BDF we need to keep M normalised.
      if(renormalise_each_time_step())
        {
          std::cout << "Renormalising nodal magnetisations." << std::endl;
          renormalise_magnetisation();
        }

#ifdef PARANOID

      // From the nmag user manual:
      // [Begin quote M Donahue]
      // * if the spin angle is approaching 180 degrees, then the results are completely bogus.
      // * over 90 degrees the results are highly questionable.
      // * Under 30 degrees the results are probably reliable.
      // [end quote]
      // (the spin angle is the angle between two neighbouring magnetisations).

      Vector<double> a = elemental_max_m_angle_variations();
      double max_angle_var = *std::max_element(a.begin(), a.end());
      if(max_angle_var > MathematicalConstants::Pi/4)
        {
          std::string error_msg
            = "Large angle variations of " + to_string(max_angle_var)
            + " > " + to_string(MathematicalConstants::Pi/4)
            + " across a single element,\n";
          error_msg += "this often means that your mesh is not sufficiently refined.";
          OomphLibWarning(error_msg, OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
        }
#endif


    }

    void actions_after_newton_step()
    {
      // Call lower level actions function
      MyProblem::actions_after_newton_step();
    }

    /// Output solution
    void doc_solution_additional(std::ofstream &some_file) const
    {
      // Number of plot points
      unsigned npts = 2;

      // Output solution with specified number of plot points per element
      bulk_mesh_pt()->output(some_file, npts);
    }

    void write_additional_trace_data(std::ofstream& trace_file) const
    {

      // Get average (and standard deviation) of |m| - 1
      double m_error_avg(0), m_error_stddev(0);
      norm_m_error(m_error_avg, m_error_stddev);

      Vector<double> angle_variations = elemental_max_m_angle_variations();
      Vector<double> mean_m = mean_magnetisation();

      trace_file
        << m_error_avg << " " // 21
        << m_error_stddev << " " // 22
        << *std::max_element(angle_variations.begin(), angle_variations.end()) << " " // 23
        << mean_m[0] << " " // 24
        << mean_m[1] << " " // 25
        << mean_m[2] << " " // 26
        << std::endl;
    }

    /// Set up an initial M
    void set_initial_condition(const InitialM::InitialMFctPt initial_m_pt);


    /// \short Return a vector of the maximum angle variation in each
    /// element.
    Vector<double> elemental_max_m_angle_variations() const
    {
      Vector<double> angles;
      for(unsigned e=0, ne=bulk_mesh_pt()->nelement(); e < ne; e++)
        {
          MicromagEquations* ele_pt = dynamic_cast<MicromagEquations*>
            (bulk_mesh_pt()->element_pt(e));
          angles.push_back(ele_pt->max_m_angle_variation());
        }
      return angles;
    }

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

    // Lots of magnetisation manipulation functions
    // ============================================================

    // Loop over all nodes in bulk mesh and get magnetisations
    Vector<Vector<double> > get_nodal_magnetisations() const
    {
      unsigned nnode = bulk_mesh_pt()->nnode();
      Vector< Vector<double> > m_list(nnode, Vector<double>(3, 0.0));

      for(unsigned nd=0; nd<nnode; nd++)
        {
          for(unsigned j=0; j<3; j++)
            {
              m_list[nd][j] = bulk_mesh_pt()->node_pt(nd)->value(m_index(j));
            }
        }

      return m_list;
    }

    void get_nodal_two_norms(Vector<double> &output) const
    {
      Vector< Vector<double> > m = get_nodal_magnetisations();
      output.assign(m.size(), 0.0);
      transform(m.begin(), m.end(), output.begin(), VectorOps::two_norm);
    }

    double mean_nodal_magnetisation_length() const
    {
      Vector<double> ms; get_nodal_two_norms(ms);
      return VectorOps::mean(ms);
    }

    Vector<double> mean_magnetisation() const
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

    /// \short Abs of mean difference of actual m and m given by a function
    /// at the middle of each element.
    double compare_m_with_function(const InitialM::InitialMFctPt fct_pt) const
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
          temp_error += std::abs(m_2norm - 1);
          // std::cout << std::abs(m_2norm - 1) << std::endl;
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
          temp_stddev += std::pow( m_2norm - 1 - m_error_avg, 2);
        }

      temp_stddev /= nnode;
      m_error_stddev = std::sqrt(temp_stddev);
    }


    // /// Elementwise calculation of m . dm/dn on the boundaries. This gives a good measure of
    // /// orthogonality of m and dmdn.
    // double mean_orthogonality_m_error() const
    // {
    //   // unsigned nele = surface_exchange_mesh_pt()->nelement();
    //   // double temp_sum = 0.0;

    //   // // Calculate m . dm/dn for each element
    //   // for(unsigned e=0; e < nele; e++)
    //   //   {
    //   //     MicromagFluxElement<MicromagEquations>* ele_pt
    //   //       = checked_dynamic_cast<MicromagFluxElement<MicromagEquations>*>
    //   //       (surface_exchange_mesh_pt()->element_pt(e));
    //   //     Vector<double> s(ele_pt->dim(), 0.5); // middle of element
    //   //     temp_sum += ele_pt->interpolated_mdotdmdn_micromag(s);
    //   //   }

    //   // // Take average
    //   // return temp_sum / double(nele);

    //   std::ostringstream error_msg;
    //   error_msg << "Not implemented.";
    //   throw OomphLibError(error_msg.str(),
    //                       OOMPH_CURRENT_FUNCTION,
    //                       OOMPH_EXCEPTION_LOCATION);

    // }

    // void orthogonality_m_error(double &orthogonality_error_avg,
    //                            double &orthogonality_error_stddev) const
    // {
    //   orthogonality_error_avg = mean_orthogonality_m_error();

    //   unsigned nele = surface_exchange_mesh_pt()->nelement();
    //   double temp_sum = 0.0;

    //   // Calculate m . dm/dn for each element
    //   for(unsigned e=0; e < nele; e++)
    //     {
    //       MicromagFluxElement<MicromagEquations>* ele_pt
    //         = checked_dynamic_cast<MicromagFluxElement<MicromagEquations>*>
    //         (surface_exchange_mesh_pt()->element_pt(e));
    //       Vector<double> s(ele_pt->dim(), 0.5); // middle of element
    //       temp_sum += pow( ele_pt->interpolated_mdotdmdn_micromag(s)
    //                        - orthogonality_error_avg, 2);
    //     }

    //   // Take stddev
    //   orthogonality_error_stddev = std::sqrt(temp_sum / double(nele));
    // }


    double get_error_norm() const
    {
      if(Compare_with_mallinson)
        {
          using namespace CompareSolutions;

          double time = ele_pt()->node_pt(0)->time_stepper_pt()->time();
          Vector<double> m_now = mean_magnetisation();
          double exact_time = switching_time_wrapper(mag_parameters_pt(), m_now);

          return std::abs(exact_time - time);
        }
      else
        {
          return -1;
        }
    }

    // Access functions
    // ============================================================

    unsigned m_index(const unsigned &j) const
    {return this->ele_pt()->m_index_micromag(j);}

    unsigned phi_index() const {return this->ele_pt()->phi_index_micromag();}

    unsigned phi_1_index() const
    {return this->ele_pt()->phi_1_index_micromag();}

    /// \short Get problem dimension (nodal dimension).
    const unsigned dim() const {return this->Dim;}

    /// \short Non-const access function for Applied_field_fct_pt.
    HApp::HAppFctPt& applied_field_fct_pt() {return Applied_field_fct_pt;}

    /// \short Const access function for Applied_field_fct_pt.
    HApp::HAppFctPt applied_field_fct_pt() const {return Applied_field_fct_pt;}

    /// \short Set function for Magnetic_parameters_pt.
    void set_mag_parameters_pt(MagneticParameters* _magnetic_parameters_pt)
    {
      Magnetic_parameters_pt = _magnetic_parameters_pt;
    }

    /// \short Const access function for Magnetic_parameters_pt.
    const MagneticParameters* mag_parameters_pt() const
    {
#ifdef PARANOID
      if(Magnetic_parameters_pt == 0)
        {
          std::string error_msg = "Magnetic paramters pointer is null.";
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return Magnetic_parameters_pt;
    }

    /// \short Non-const access function for Renormalise_each_time_step.
    bool& renormalise_each_time_step() {return Renormalise_each_time_step;}

    /// \short Const access function for Renormalise_each_time_step.
    bool renormalise_each_time_step() const {return Renormalise_each_time_step;}

    //     /// \short Non-const access function for Surface_exchange_mesh_pt.
    //     Mesh*& surface_exchange_mesh_pt() {return Surface_exchange_mesh_pt;}

    //     /// \short Const access function for Surface_exchange_mesh_pt.
    //     Mesh* surface_exchange_mesh_pt() const
    //     {
    // #ifdef PARANOID
    //       if(Surface_exchange_mesh_pt == 0)
    //         {
    //           std::ostringstream error_msg;
    //           error_msg << "Surface exchange mesh pointer not set.";
    //           throw OomphLibError(error_msg.str(),
    //                               OOMPH_CURRENT_FUNCTION,
    //                               OOMPH_EXCEPTION_LOCATION);
    //         }
    // #endif

    //       return Surface_exchange_mesh_pt;
    //     }

    /// \short Non-const access function for Bulk_mesh_pt.
    void set_bulk_mesh_pt(Mesh* mesh_pt) {Bulk_mesh_pt = mesh_pt;}

    /// \short Const access function for Bulk_mesh_pt.
    Mesh* bulk_mesh_pt() const
    {
#ifdef PARANOID
      if(Bulk_mesh_pt == 0)
        {
          std::string error_msg = "Null bulk mesh pointer!";
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return Bulk_mesh_pt;
    }

    /// \short Get a pointer to a MicromagEquations element.
    MicromagEquations* ele_pt() const
    {
      return checked_dynamic_cast<MicromagEquations*>
        (bulk_mesh_pt()->element_pt(0));
    }

    /// Can we check the solution using Mallinson's exact time + phi
    /// solutions?
    bool Compare_with_mallinson;

    /// \short Should we swap to superlu for large dt solves?
    bool Swap_solver_large_dt;

  private:

    /// \short Storage for initial linear solver: in case we want to swap
    /// to superlu for large dt.
    LinearSolver* My_linear_solver_pt;

    /// Pointer to the applied field.
    HApp::HAppFctPt Applied_field_fct_pt;

    /// Magnetic parameters storage. ??ds should maybe go in meshes?
    MagneticParameters* Magnetic_parameters_pt;

    /// Normalise magnetisation problem after each step?
    bool Renormalise_each_time_step;

    /// Pointer to bulk mesh (i.e. the magnetic material).
    Mesh* Bulk_mesh_pt;

    // /// Pointer to the mesh of face elements for dealing with boundary
    // /// conditions caused by reducing order of differentiation in exchange
    // /// term.
    // Mesh* Surface_exchange_mesh_pt;

    /// Create
    void create_surface_exchange_elements(const unsigned& b);

    /// Inaccessible copy constructor
    LLGProblem(const LLGProblem & dummy)
    {BrokenCopy::broken_copy("LLGProblem");}

    /// Inaccessible assignment operator
    void operator=(const LLGProblem &dummy)
    {BrokenCopy::broken_assign("LLGProblem");}

  };

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
    MicromagEquations* elem_pt = dynamic_cast<MicromagEquations* >(this->mesh_pt()->element_pt(0));
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
  }



  namespace LLGFactories
  {
    /// \short Make a mesh as specified by an input argument. Refined
    /// according to the given refinement level (in some way appropriate
    /// for that mesh type). Assumption: this will be passed into a
    /// problem, which will delete the pointer when it's done.
    Mesh* mesh_factory(const std::string& _mesh_name,
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

      // For some reason we have to call this manually...
      mesh_pt->setup_boundary_element_info();

      // Done: pass out the mesh pointer
      return mesh_pt;
    }

  }

  /// Base command line args processing class for pure llg and semi
  /// implicit (and any other magnetism problems).
  class MMArgs : public MyCliArgs
  {
  public:
    /// Constructor: Initialise pointers to null.
    MMArgs() : initial_m_fct_pt(0), h_app_fct_pt(0),
               magnetic_parameters_pt(0) {}

    virtual void set_flags()
    {
      MyCliArgs::set_flags();

      specify_command_line_flag("-initm", &initial_m_name);
      initial_m_name = "z";

      specify_command_line_flag("-happ", &h_app_name);
      h_app_name = "minus_z";

      specify_command_line_flag("-mag-params", &magnetic_parameters_name);
      magnetic_parameters_name = "simple-llg";

      specify_command_line_flag("-renorm_m", &Renormalise);
      Renormalise = -1;
    }


    virtual void run_factories()
    {
      MyCliArgs::run_factories();

      initial_m_name = to_lower(initial_m_name);
      h_app_name = to_lower(h_app_name);
      magnetic_parameters_name = to_lower(magnetic_parameters_name);

      initial_m_fct_pt = InitialM::initial_m_factory(initial_m_name);
      h_app_fct_pt = HApp::h_app_factory(h_app_name);
      magnetic_parameters_pt =
        Factories::magnetic_parameters_factory(magnetic_parameters_name);
    }

    /// Write out all args (in a parseable format) to a stream.
    virtual void dump_args(std::ostream& out_stream) const
    {
      MyCliArgs::dump_args(out_stream);

      out_stream
        << "initial_m " << initial_m_name << std::endl
        << "h_app " << h_app_name << std::endl
        << "mag_params " << magnetic_parameters_name << std::endl;
    }


    bool renormalise_flag()
    {
      // If flag not set then only do it for bdf timesteppers
      if(Renormalise == -1)
        {
          return (time_stepper_name == "bdf2")
            || (time_stepper_name == "bdf1");
        }
      // Otherwise do what the flag says
      else
        {
          return bool(Renormalise);
        }
    }

    InitialM::InitialMFctPt initial_m_fct_pt;
    HApp::HAppFctPt h_app_fct_pt;
    MagneticParameters* magnetic_parameters_pt;

    // Strings for input to factory functions
    std::string initial_m_name;
    std::string h_app_name;
    std::string magnetic_parameters_name;


    /// Flag to control renormalisation of |m| after each step. -1 =
    /// default for timestepper, 0 = off, 1 = on.
    int Renormalise;
  };


  /// Command line args class for llg problems. Just add the mesh
  /// stuff. ??ds refactor to combine with MMArgs?
  class LLGArgs : public MMArgs
  {
  public:

    /// Constructor: Initialise pointers to null.
    LLGArgs() : mesh_pt(0) {}

    virtual void set_flags()
    {
      MMArgs::set_flags();

      specify_command_line_flag("-mesh", &mesh_name);
      mesh_name = "sq_square";
    }


    virtual void run_factories()
    {
      MMArgs::run_factories();

      mesh_name = to_lower(mesh_name);

      // Do the mesh last of all because it can be slow
      mesh_pt = LLGFactories::mesh_factory(mesh_name, refinement, time_stepper_pt);
    }

    /// Write out all args (in a parseable format) to a stream.
    virtual void dump_args(std::ostream& out_stream) const
    {
      MMArgs::dump_args(out_stream);
      out_stream << "mesh " << mesh_name << std::endl;
    }


    Mesh* mesh_pt;


    // Strings for input to factory functions
    std::string mesh_name;
  };

}

#endif
