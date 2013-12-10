#ifndef OOMPH_SEMI_IMPLICIT_PROBLEM_H
#define OOMPH_SEMI_IMPLICIT_PROBLEM_H

/*
  description of file goes here
*/

#include "boundary_element_handler.h"
#include "generic_poisson_problem.h"
#include "llg_problem.h"

#include "./micromagnetics_boundary_element.h"
#include "./magnetostatic_field_flux_element.h"


namespace oomph
{

  using namespace StringConversion;
  using namespace CommandLineArgs;


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
                           double scaling_factor = 1.0,
                           unsigned nnode1d = 2);


    /// \short Make a mesh of MagnetostaticField elements as specified by an
    /// input argument. Refined according to the given refinement level (in
    /// some way appropriate for that mesh type).
    Mesh* phi_mesh_factory(const std::string& _mesh_name,
                           int refinement_level,
                           TimeStepper* time_stepper_pt,
                           double scaling_factor = 1.0,
                           unsigned nnode1d = 2);


    /// \short Return a factory function which will create the appropriate
    /// "flux mesh" for the bulk element pointer given.
    GenericPoissonProblem::FluxMeshFactoryFctPt
    phi_1_flux_mesh_factory_factory(const FiniteElement* bulk_phi_1_ele_pt);
  }


  /// Base class for problems which solve the llg with magnetostatics as
  /// two separate problems at each step: one Poisson BEM problen for the
  /// field then a normal llg problem for the magnetisation.
  class DecoupledLLGProblem : public LLGProblem
  {
  public:

    // Default constructor
    DecoupledLLGProblem() : Phi_boundary_values_pts()
    {
      Bem_handler_pt = 0;
      Phi_1_problem_pt = 0;
      Phi_problem_pt = 0;
    }

    /// \short Function to do the real work of the constructor.
    void build(Vector<Mesh*>& llg_mesh_pts, Vector<Mesh*>& phi_mesh_pts,
               Vector<Mesh*>& phi_1_mesh_pts, bool pin_phi1=true);

    /// Destructor
    virtual ~DecoupledLLGProblem()
    {
      // Kill boundary value storage vectors
      for(unsigned j=0; j<Phi_boundary_values_pts.size(); j++)
        {
          delete Phi_boundary_values_pts[j];
        }
    }

    /// Interface for whichever function is needed to do the timestepping
    virtual double do_step(const double& dt, const double& tol)=0;

    /// \short Solve for the magnetostatic field.
    void magnetostatics_solve()
    {
      std::cout << std::endl
                << "BEM solve" << std::endl
                << "--------------------------" <<std::endl;

      // solve for phi1
      std::cout << "solving phi1" << std::endl;
      phi_1_problem_pt()->newton_solve();

      // update boundary values of phi
      std::cout << "solving BEM" << std::endl;
      double t_start = TimingHelpers::timer();
      Bem_handler_pt->get_bem_values(Phi_boundary_values_pts);
      double t_end = TimingHelpers::timer();
      std::cout << "BEM time taken: " << t_end - t_start << std::endl;

      // push old phi values back in time (so that we can use them later to
      // get time derivatives of the field). Note that we don't use the
      // problem's shift time values function because we don't want to
      // shift the timestepper (that has been done by the llg problem
      // already) and we don't have any external data to shift.
      phi_problem_pt()->mesh_pt()->shift_time_values();

      // solve for phi
      std::cout << "solving phi" << std::endl;
      phi_problem_pt()->newton_solve();

      std::cout << "mean field is " << average_magnetostatic_field() << std::endl;
    }

    /// Set up an initial M
    void set_initial_condition(const InitialM::InitialMFctPt initial_m_pt);

    /// Initialise timestep: only llg problem has a timestep.
    void initialise_dt(const double &dt)
    {llg_sub_problem_pt()->initialise_dt(dt);}

    /// Output
    void doc_solution_additional(std::ofstream &some_file) const;

    Vector<double> average_magnetostatic_field() const;

    /// Access functions
    // =================================================================

    /// Get pointer to an LLG element for looking up info. All elements in
    /// mesh should be the same type otherwise preconditioning framework
    /// will fail so this should be safe.
    MicromagEquations* llg_element_pt() const
    {return checked_dynamic_cast<MicromagEquations*>
        (llg_mesh_pt()->element_pt(0));}

    /// Get access to magnetic parameters - only relevant in LLG problem so
    /// return the pointer from that problem.
    const MagneticParameters* mag_parameters_pt() const
    {return llg_sub_problem_pt()->mag_parameters_pt();}

    /// \short Const access function for LLG_problem.
    const LLGProblem* llg_sub_problem_pt() const
    {return this;}

    /// \short Non-const acess function for LLG_problem.
    LLGProblem* llg_sub_problem_pt() {return this;}

    const Mesh* bem_mesh_pt() const {return Bem_handler_pt->bem_mesh_pt();}

    const Mesh* llg_mesh_pt() const
    {
      return llg_sub_problem_pt()->mesh_pt();
    }

    const Mesh* phi_1_mesh_pt() const {return phi_1_problem_pt()->mesh_pt();}
    const Mesh* phi_mesh_pt() const {return phi_problem_pt()->mesh_pt();}

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
    // LLGProblem LLG_problem;

    /// Intermediate storage for results of bem (ideally we would have it
    /// call a function to get the boundary values filled in but c++ member
    /// functions pointers are useless...)
    Vector<DoubleVector*> Phi_boundary_values_pts;
  };


  /// Class for solving the implicit llg with magnetostatics handled
  /// explicitly.
  class SemiImplicitHybridMicromagneticsProblem : public DecoupledLLGProblem
  {
    public:
    SemiImplicitHybridMicromagneticsProblem() : DecoupledLLGProblem() {}

    /// Replacement for "newton_solve()" that does a few different solves.
    double semi_implicit_step(const double &dt, const double eps)
    {
      // Solve for the magnetostatic field.
      magnetostatics_solve();

      // solve for m
      if(Use_time_adaptive_newton)
        {
          return llg_sub_problem_pt()->adaptive_unsteady_newton_solve(dt, eps);
        }
      else
        {
          llg_sub_problem_pt()->unsteady_newton_solve(dt);
        }
      return dt;
      //??ds make this smarter?
    }

    /// Function to do a time step: just call semi-implicit step
    double do_step(const double& dt, const double& tol)
      {
        return semi_implicit_step(dt, tol);
      }

  };


  /// Class for an llg problem which can be solved by an explicit time
  /// stepper. This has to be a separate class from
  /// SemiImplicitHybridMicromagneticsProblem because when using midpoint
  /// method we often want to take a single explcit time step as a
  /// predictor, but in that case we don't need to solve the magnetostatics
  /// problem again (because we just did it).
  class ExplicitLLGProblem : public DecoupledLLGProblem
  {
  public:
    ExplicitLLGProblem() : DecoupledLLGProblem() {}

    /// Do an explicit step step (we need extra things because we need to
    /// solve for phi before we can do the normal explicit step).
    void get_inverse_mass_matrix_times_residuals(DoubleVector &Mres)
    {
      // Solve for the magnetostatic field.
      magnetostatics_solve();

      // Now get what we actually wanted
      Problem::get_inverse_mass_matrix_times_residuals(Mres);
    }

    /// Function to do a time step: just call explicit step
    double do_step(const double& dt, const double& tol)
    {
      explicit_timestep(dt);
      return dt;
    }
  };



  /// Command line args class for semi implicit llg problems. Just add the
  /// mesh stuff.
  class SemiImplicitMMArgs : public MMArgs
  {

  public:
    SemiImplicitMMArgs() {}

    virtual void set_flags()
    {
      MMArgs::set_flags();

      specify_command_line_flag("-mesh", &mesh_name);
      mesh_name = "sq_square";

      specify_command_line_flag("-xshift", &xshift);
      xshift = 1.5;

      specify_command_line_flag("-yshift", &yshift);
      yshift = 1.5;

      specify_command_line_flag("-numerical-BEM");
      // automatically defaults to false
    }


    void run_factories()
    {
      MMArgs::run_factories();

      // Store inside class
      use_numerical_integration = Specified_command_line_flag["-numerical-BEM"];

      to_lower(mesh_name);

      // Build the meshes, do this last because they can be SLOW, must be
      // done before factory mesh function selection...

      // If the mesh name is prefixed by "multi_" then build the equivalent
      // (simple, for now?) multiple mesh.
      if(Factories::has_prefix("multi_", mesh_name))
        {
          // Copy the rest of the mesh name into a new string
          std::string single_mesh_name = Factories::rest_of_name("multi_",
                                                                 mesh_name);

          // And use it to build the meshes
          llg_mesh_pts = Factories::simple_multimesh_factory
            (SemiImplicitFactories::llg_mesh_factory,
             single_mesh_name, refinement, time_stepper_pt, xshift);

          phi_mesh_pts = Factories::simple_multimesh_factory
            (SemiImplicitFactories::phi_mesh_factory,
             single_mesh_name, refinement, time_stepper_pt, xshift);

          phi_1_mesh_pts = Factories::simple_multimesh_factory
            (SemiImplicitFactories::phi_mesh_factory,
             single_mesh_name, refinement, time_stepper_pt, xshift);

        }

      // Or with "many" prefix make a load of meshes
      else if(Factories::has_prefix("many_", mesh_name))
        {
          std::string base_name = Factories::rest_of_name("many_", mesh_name);

          llg_mesh_pts = Factories::simple_many_multimesh_factory
            (SemiImplicitFactories::llg_mesh_factory, base_name, refinement,
             time_stepper_pt, xshift, yshift);
          phi_mesh_pts = Factories::simple_many_multimesh_factory
            (SemiImplicitFactories::phi_mesh_factory, base_name, refinement,
             time_stepper_pt, xshift, yshift);
          phi_1_mesh_pts = Factories::simple_many_multimesh_factory
            (SemiImplicitFactories::phi_mesh_factory, base_name, refinement,
             time_stepper_pt, xshift, yshift);
         }

      // Otherwise just build one mesh
      else
        {
          // LLG (magnetism) mesh
          llg_mesh_pts.push_back(SemiImplicitFactories::llg_mesh_factory
                                 (mesh_name, refinement, time_stepper_pt));

          // Make the two phi meshes
          phi_mesh_pts.push_back(SemiImplicitFactories::phi_mesh_factory
                                 (mesh_name, refinement, time_stepper_pt));
          phi_1_mesh_pts.push_back(SemiImplicitFactories::phi_mesh_factory
                                   (mesh_name, refinement, time_stepper_pt));
        }


      // Pick the factory function for creating the phi 1 surface mesh
      phi_1_flux_mesh_factory_fct_pt =
        SemiImplicitFactories::phi_1_flux_mesh_factory_factory
        (phi_1_mesh_pts[0]->finite_element_pt(0));

      // Pick the factory function for creating the BEM elements
      bem_element_factory_fct_pt =
        LLGFactories::bem_element_factory_factory
        (llg_mesh_pts[0]->finite_element_pt(0));


      // Pick the problem itself (explicit or semi implicit)
      if(explicit_flag())
        {
          problem_pt = new ExplicitLLGProblem;
        }
      else
        {
          problem_pt = new SemiImplicitHybridMicromagneticsProblem;
        }
    }

    /// Write out all args (in a parseable format) to a stream.
    virtual void dump_args(std::ostream& out_stream) const
    {
      MMArgs::dump_args(out_stream);
      out_stream << "mesh " << mesh_name << std::endl;
      out_stream << "numerical-BEM " << use_numerical_integration
                 << std::endl;
    }


    Vector<Mesh*> llg_mesh_pts;
    Vector<Mesh*> phi_1_mesh_pts;
    Vector<Mesh*> phi_mesh_pts;

    GenericPoissonProblem::FluxMeshFactoryFctPt phi_1_flux_mesh_factory_fct_pt;
    LLGFactories::BEMElementFactoryFctPt bem_element_factory_fct_pt;

    // Strings for input to factory functions
    std::string mesh_name;

    // Distance to move meshes away from origin (along +- x) in multi-mesh.
    double xshift;
    double yshift;

    bool use_numerical_integration;

    DecoupledLLGProblem* problem_pt;
  };


} // End of oomph namespace

#endif
