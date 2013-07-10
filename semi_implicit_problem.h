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

    /// \short Create a variable order quadrature object based on the
    /// dimension and shape of the element. Only works for
    Integral* variable_order_integrator_factory(const FiniteElement* const el_pt);

    /// \short Make a mesh of Micromag elements as specified by an
    /// input argument. Refined according to the given refinement level (in
    /// some way appropriate for that mesh type).
    Mesh* llg_mesh_factory(const std::string& _mesh_name,
                           int refinement_level,
                           TimeStepper* time_stepper_pt,
                           unsigned nnode1d = 2);


    /// \short Make a mesh of MagnetostaticField elements as specified by an
    /// input argument. Refined according to the given refinement level (in
    /// some way appropriate for that mesh type).
    Mesh* phi_mesh_factory(const std::string& _mesh_name,
                           int refinement_level,
                           TimeStepper* time_stepper_pt,
                           unsigned nnode1d = 2);

    /// \short Function pointer type for function which returns a BEM
    /// element.
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

    /// \short Return a factory function which will create the appropriate
    /// "flux mesh" for the bulk element pointer given.
    GenericPoissonProblem::FluxMeshFactoryFctPt
    phi_1_flux_mesh_factory_factory(const FiniteElement* bulk_phi_1_ele_pt);

    /// \short Return a function which will create the appropriate BEM face
    /// element for the bulk element pointer given (should work for a
    /// pointer to any bulk element type i.e., field or llg).
    BEMElementFactoryFctPt bem_element_factory_factory
    (const FiniteElement* bulk_ele_pt);


  }


  class SemiImplicitHybridMicromagneticsProblem :
    public LLGProblem
  {
  public:

    // Default constructor, who knows what will happen here... ??ds
    SemiImplicitHybridMicromagneticsProblem() :
      Bem_handler_pt(0), Phi_1_problem_pt(), Phi_problem_pt(),
      Phi_boundary_values_pts()
    {}

    /// \short Function to do the real work of the constructor.
    void build(bool pin_phi1 = true);

    /// Destructor
    ~SemiImplicitHybridMicromagneticsProblem()
    {
      // Kill boundary value storage vectors
      for(unsigned j=0; j<Phi_boundary_values_pts.size(); j++)
        {
          delete Phi_boundary_values_pts[j];
        }
    }


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
      Bem_handler_pt->get_bem_values(Phi_boundary_values_pts);

      // solve for phi
      std::cout << "solving phi" << std::endl;
      phi_problem_pt()->newton_solve();
    }


    /// Replacement for "newton_solve()" that does a few different solves.
    double semi_implicit_step(const double &dt, const double eps=0.0)
    {
      // Solve for the magnetostatic field.
      magnetostatics_solve();

      // solve for m
      if(eps != 0.0)
        {
          std::cout << "solving LLG (with time adaptivity)" << std::endl;
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
    void doc_solution_additional(std::ofstream &some_file) const;

    Vector<double> average_magnetostatic_field() const;

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
    // LLGProblem LLG_problem;

    /// Intermediate storage for results of bem (ideally we would have it
    /// call a function to get the boundary values filled in but c++ member
    /// functions pointers are useless...)
    Vector<DoubleVector*> Phi_boundary_values_pts;
  };


  /// Command line args class for semi implicit llg problems. Just add the
  /// mesh stuff.
  class SemiImplicitMMArgs : public MMArgs
  {

  public:
    SemiImplicitMMArgs() : llg_mesh_pt(0), phi_1_mesh_pt(0), phi_mesh_pt(0) {}

    virtual void set_flags()
    {
      MMArgs::set_flags();

      specify_command_line_flag("-mesh", &mesh_name);
      mesh_name = "sq_square";

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

      // LLG (magnetism) mesh
      llg_mesh_pt = SemiImplicitFactories::llg_mesh_factory
        (mesh_name, refinement, time_stepper_pt);

      // Make the two phi meshes
      phi_mesh_pt = SemiImplicitFactories::phi_mesh_factory
        (mesh_name, refinement, time_stepper_pt);
      phi_1_mesh_pt = SemiImplicitFactories::phi_mesh_factory
        (mesh_name, refinement, time_stepper_pt);

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
      MMArgs::dump_args(out_stream);
      out_stream << "mesh " << mesh_name << std::endl;
      out_stream << "numerical-BEM " << use_numerical_integration
                 << std::endl;
    }


    Mesh* llg_mesh_pt;
    Mesh* phi_1_mesh_pt;
    Mesh* phi_mesh_pt;

    GenericPoissonProblem::FluxMeshFactoryFctPt phi_1_flux_mesh_factory_fct_pt;
    SemiImplicitFactories::BEMElementFactoryFctPt bem_element_factory_fct_pt;

    // Strings for input to factory functions
    std::string mesh_name;

    bool use_numerical_integration;
  };


} // End of oomph namespace

#endif
