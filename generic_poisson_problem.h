#ifndef OOMPH_GENERIC_POISSON_PROBLEM_H
#define OOMPH_GENERIC_POISSON_PROBLEM_H


/*
  TODO:
  * solution checking/output/fn pt storage
  * static const for magic number

  * DirichletFctPt -> something better

  */

#include "./template_free_poisson.h"
#include "./template_free_poisson_flux.h"

#include "./my_general_header.h"
#include "./my_generic_problem.h"



namespace oomph
{

  /// Simple struct to store info about Neumann boundaries.
  struct NeumannCondition
  {
    /// Typedef because the full name is far too long
    typedef typename TFPoissonFluxEquations::PoissonPrescribedFluxFctPt
    PoissonFluxFctPt;

    Mesh* mesh_pt;
    unsigned boundary;
    PoissonFluxFctPt prescribed_flux_pt;
  };

  /// Simple struct to store info about Dirichlet boundaries.
  struct DirichletCondition
  {
    /// Use the same function pointer type to store Diriclet conditions as is
    /// used to store exact solutions.
    typedef FiniteElement::SteadyExactSolutionFctPt DirichletFctPt;

    Mesh* mesh_pt;
    unsigned boundary;
    DirichletFctPt dc_fct_pt;
    const DoubleVector* dc_vector_pt;
  };

  class GenericPoissonProblem : public MyProblem
  {

  public:

    /// Use the same function pointer type to store Diriclet conditions as is
    /// used to store exact solutions.
    typedef FiniteElement::SteadyExactSolutionFctPt DirichletFctPt;

    /// Typedef because the full name is far too long
    typedef typename TFPoissonFluxEquations::PoissonPrescribedFluxFctPt
    PoissonFluxFctPt;

    typedef typename TFPoissonEquations::PoissonSourceFctPt SourceFctPt;

    typedef Mesh* (*FluxMeshFactoryFctPt)(Mesh* bulk_mesh_pt,
                                          const Vector<unsigned> &boundaries);


    // Note: typename is needed here and in some places below because the
    // definitions of the things after "typename" are in other files and
    // the compiler needs to know that they are types not functions (I
    // think...).

    /// Constructor
    GenericPoissonProblem()
      : Flux_mesh_factory_fct_pt(0), Poisson_dof_number(0), Source_fct_pt(0),
        Exact_solution_fct_pt(0) {}

    /// Destructor
    virtual ~GenericPoissonProblem() {}

    /// Doc the solution.
    void doc_solution_additional(std::ofstream &soln_file) const;

    /// \short Create flux mesh on boundaries of the mesh, with the flux
    /// given by the output of the function pointer pointer. ??ds currently
    /// does one boundary at a time, change to all boundaries on single
    /// mesh in one go?
    void build_neumann_boundary_mesh(Mesh* mesh_pt,
                                     const unsigned& boundary,
                                     PoissonFluxFctPt const prescribed_flux_pt)
    {
      // Build the mesh
      Mesh* flux_mesh_pt = flux_mesh_factory(mesh_pt,
                                             Vector<unsigned>(1, boundary));

      // Set the prescribed flux function on the elements
      for(unsigned e=0, ne=flux_mesh_pt->nelement(); e < ne; e++)
        {
          TFPoissonFluxEquations* ele_pt
            = checked_dynamic_cast<TFPoissonFluxEquations*>(flux_mesh_pt->element_pt(e));
          ele_pt->flux_fct_pt() = prescribed_flux_pt;
        }

      // Add to global mesh list
      add_sub_mesh(flux_mesh_pt);
    }

    void add_neumann_boundary(Mesh* _mesh_pt,
                              const unsigned& _boundary,
                              PoissonFluxFctPt const _prescribed_flux_pt)
      {
        NeumannCondition n;
        n.mesh_pt = _mesh_pt;
        n.boundary = _boundary;
        n.prescribed_flux_pt = _prescribed_flux_pt;
        Neumann_conditions.push_back(n);
      }

    /// \short Call the stored Flux_mesh_factory_fct_pt (and check if null
    /// in PARANOID).
    Mesh* flux_mesh_factory(Mesh* bulk_mesh_pt,
                            const Vector<unsigned> &boundaries) const
    {
#ifdef PARANOID
      if(Flux_mesh_factory_fct_pt == 0)
        {
          std::string error_msg = "No flux mesh factory function pointer has been set.";
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return Flux_mesh_factory_fct_pt(bulk_mesh_pt, boundaries);
    }

    /// \short Set function for Flux_mesh_factory.
    void set_flux_mesh_factory(FluxMeshFactoryFctPt flux_mesh_factory_fct_pt)
    { Flux_mesh_factory_fct_pt = flux_mesh_factory_fct_pt; }

    /// \short Pin nodes on boundary b of bulk mesh with the value given by a
    /// function pointer.
    void add_dirichlet_boundary(Mesh* mesh_pt,
                                const unsigned &boundary,
                                const DirichletFctPt dc_fct_pt)
    {
      DirichletCondition dcb;
      dcb.mesh_pt = mesh_pt;
      dcb.boundary = boundary;
      dcb.dc_fct_pt = dc_fct_pt;
      dcb.dc_vector_pt = 0;

      Dirichlet_conditions.push_back(dcb);
    }

    /// \short Pin nodes on boundary b of bulk mesh with the value given by
    /// a vector of values.
    void add_dirichlet_boundary_by_vector(Mesh* mesh_pt,
                                          const unsigned& b,
                                          const DoubleVector* dc_vector_pt)
    {
      DirichletCondition dcb;
      dcb.mesh_pt = mesh_pt;
      dcb.boundary = b;
      dcb.dc_fct_pt = 0;
      dcb.dc_vector_pt = dc_vector_pt;

      Dirichlet_conditions.push_back(dcb);
    }


    /// Finish building the problem once everything has been set.
    void build(Vector<Mesh*>& bulk_mesh_pts);

    /// Update dirichlet conditions before Newton solve.
    virtual void actions_before_newton_solve()
    {
      // Call base class version
      MyProblem::actions_before_newton_solve();

      update_dirichlet_conditions();
    }

    /// (Re-)apply the Dirichlet conditions on the boundaries which have
    /// been set as Dirichlet. Note: if values differ where two boundaries
    /// overlap then priority is 1) vector conditions > function
    /// conditions, 2) higher boundary number conditions > lower boundary
    /// number conditions. Might want to change this...
    void update_dirichlet_conditions();

    /// Get the norm of the error (if exact solution has been set and problem
    /// has been solved).
    double get_error_norm() const;

    // Access functions
    // ============================================================

    /// Set function for Source_fct_pt.
    void set_source_fct_pt(SourceFctPt const source_fct_pt)
    {Source_fct_pt = source_fct_pt;}

    /// Const access function for Source_fct_pt.
    SourceFctPt source_fct_pt() const
    {
      if(Source_fct_pt == 0)
        {
          std::ostringstream error_msg;
          error_msg << "Source function pointer for this problem "
                    << "has not been set yet.";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      return Source_fct_pt;
    }

    /// Non-const access function for Poisson_dof_number.
    unsigned& poisson_dof_number() {return Poisson_dof_number;}

    /// Const access function for Poisson_dof_number.
    unsigned poisson_dof_number() const {return Poisson_dof_number;}

    /// \short Non-const access function for Exact_solution_fct_pt.
    FiniteElement::SteadyExactSolutionFctPt& exact_solution_fct_pt()
    {return Exact_solution_fct_pt;}

    /// \short Const access function for Exact_solution_fct_pt.
    FiniteElement::SteadyExactSolutionFctPt exact_solution_fct_pt() const
    {return Exact_solution_fct_pt;}


  private:

    /// \short Function pointer for the function to use to construct a flux
    /// mesh on boundaries of the bulk mesh.
    FluxMeshFactoryFctPt Flux_mesh_factory_fct_pt;

    /// Storage for flux meshes
    Vector<NeumannCondition> Neumann_conditions;

    /// Store the index of the Poisson dof (defaults to zero).
    unsigned Poisson_dof_number;

    /// \short Storage for a pointer to the source function of the Poisson
    /// equations.
    SourceFctPt Source_fct_pt;

    /// \short Storage for a pointer to a function giving the exact solution
    /// (for validation).
    FiniteElement::SteadyExactSolutionFctPt Exact_solution_fct_pt;

    /// \short Storage for a list of Dirichlet boundaries and pointers to
    /// the functions or vectors which define their values.
    Vector<DirichletCondition> Dirichlet_conditions;

    /// Inaccessible copy constructor
    GenericPoissonProblem(const GenericPoissonProblem& dummy)
    {BrokenCopy::broken_copy("GenericPoissonProblem");}

    /// Inaccessible assignment operator
    void operator=(const GenericPoissonProblem& dummy)
    {BrokenCopy::broken_assign("GenericPoissonProblem");}

  };




} // End of oomph namespace

#endif
