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


namespace oomph
{

  class GenericPoissonProblem : public Problem
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
    GenericPoissonProblem() :
      Flux_mesh_factory_fct_pt(0), Bulk_mesh_number(-10), Flux_mesh_number(-10),
      Poisson_dof_number(0), Source_fct_pt(0), Exact_solution_fct_pt(0)
    {}

    /// Destructor
    ~GenericPoissonProblem() {}

    /// Doc the solution.
    virtual void doc_solution(DocInfo& doc_info) const;

    /// \short Create flux mesh on boundaries of bulk mesh, with the flux
    /// given by the output of the function pointer pointer.
    void set_neumann_boundaries(const Vector<unsigned> &boundaries,
                                PoissonFluxFctPt const prescribed_flux_pt)
    {
#ifdef PARANOID
      if(Flux_mesh_number != -10)
        {
          std::string error_msg = "Already have a flux mesh!";
          throw OomphLibError(error_msg, OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

      // Build the mesh
      Mesh* flux_mesh_pt = flux_mesh_factory(bulk_mesh_pt(), boundaries);

      // Set the prescribed flux function on the elements
      for(unsigned e=0, ne=flux_mesh_pt->nelement(); e < ne; e++)
        {
          TFPoissonFluxEquations* ele_pt
            = checked_dynamic_cast<TFPoissonFluxEquations*>(flux_mesh_pt->element_pt(e));
          ele_pt->flux_fct_pt() = prescribed_flux_pt;
        }

      // Store the mesh and its location in the list of meshes.
      Flux_mesh_number = add_sub_mesh(flux_mesh_pt);
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
    void set_dirichlet_boundary(const unsigned &b,
                                const DirichletFctPt condition_fct_pt)
    {
      // Store function pointer to compute values on boundary b
      std::pair<unsigned, DirichletFctPt> dcb;
      dcb.first = b;
      dcb.second = condition_fct_pt;
      Dirichlet_function_conditions.push_back(dcb);
    }


    void set_dirichlet_boundary_by_vector(const unsigned& b,
                                          const DoubleVector* boundary_values_pt)
    {
      // Store pointer to vector where boundary values will be put
      std::pair<unsigned, const DoubleVector*> dcb;
      dcb.first = b;
      dcb.second = boundary_values_pt;
      Dirichlet_vector_conditions.push_back(dcb);
    }


    /// Finish building the problem once everything has been set.
    void build();

    /// Update dirichlet conditions before Newton solve.
    virtual void actions_before_newton_solve()
    {
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
    void set_source_fct_pt
    (SourceFctPt const source_fct_pt)
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

    /// Add a bulk mesh and set up numbering
    void set_bulk_mesh(Mesh* b_mesh)
    {
#ifdef PARANOID
      if(Bulk_mesh_number != -10)
        {
          std::ostringstream error_msg;
          error_msg << "Already have a bulk mesh for this problem."
                    << "It's probably possible to change it but might be messy...";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      Bulk_mesh_number = add_sub_mesh(b_mesh);
    }

    /// Const access function for Bulk_mesh_pt.
    Mesh* bulk_mesh_pt() const
    {
#ifdef PARANOID
      if(Bulk_mesh_number == -10)
        {
          std::ostringstream error_msg;
          error_msg << "No bulk mesh set for this problem.";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return mesh_pt(Bulk_mesh_number);
    }

    /// Const access function for Neumann_condition_mesh_pt.
    Mesh* flux_mesh_pt() const
    {
#ifdef PARANOID
      if(Flux_mesh_number == -10)
        {
          std::ostringstream error_msg;
          error_msg << "No flux mesh set for this problem.";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      return mesh_pt(Flux_mesh_number);
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

    /// \short Store the number of the bulk mesh in the global mesh_pt
    /// array. Defaults to magic number -10 so that we can tell if it is
    /// unset.
    int Bulk_mesh_number;

    /// \short Store the number of the flux (Neumann condition enforcing)
    /// mesh in the global mesh_pt array. Defaults to magic number -10 so
    /// that we can tell if it is unset.
    int Flux_mesh_number;

    /// Store the index of the Poisson dof (defaults to zero).
    unsigned Poisson_dof_number;

    /// \short Storage for a pointer to the source function of the Poisson
    /// equations.
    SourceFctPt Source_fct_pt;

    /// \short Storage for a pointer to a function giving the exact solution
    /// (for validation).
    FiniteElement::SteadyExactSolutionFctPt Exact_solution_fct_pt;

    /// \short Storage for a list of Dirichlet boundaries and pointers to
    /// the functions which define their values.
    Vector<std::pair<unsigned, DirichletFctPt> > Dirichlet_function_conditions;

    /// \short Storage for a list of Dirichlet boundaries and pointers to
    /// the vectors which define their values.
    Vector<std::pair<unsigned, const DoubleVector*> > Dirichlet_vector_conditions;

    /// Inaccessible copy constructor
    GenericPoissonProblem(const GenericPoissonProblem& dummy)
    {BrokenCopy::broken_copy("GenericPoissonProblem");}

    /// Inaccessible assignment operator
    void operator=(const GenericPoissonProblem& dummy)
    {BrokenCopy::broken_assign("GenericPoissonProblem");}

  };




} // End of oomph namespace

#endif
