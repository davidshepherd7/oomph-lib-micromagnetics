#ifndef OOMPH_GENERIC_POISSON_PROBLEM_H
#define OOMPH_GENERIC_POISSON_PROBLEM_H


/*
  TODO:
  * solution checking/output/fn pt storage
  * static const for magic number

  * DirichletFctPt -> something better

  */

#include "generic.h"

#include "poisson.h"

using namespace oomph;

namespace oomph
{

 template<class ELEMENT>
  class GenericPoissonProblem : public Problem
 {

 public:

  /// Use the same function pointer type to store Diriclet conditions as is
  /// used to store exact solutions.
  typedef FiniteElement::SteadyExactSolutionFctPt DirichletFctPt;

  /// Typedef because the full name is far too long
  typedef typename PoissonFluxElement<ELEMENT>::PoissonPrescribedFluxFctPt
   PoissonFluxFctPt;

  // Note: typename is needed here and in some places below because the
  // definitions of the things after "typename" are in other files and the
  // compiler needs to know that they are types not classes or
  // functions. (I think...)

  /// Constructor: Pass pointer to source function
 GenericPoissonProblem() :
  Bulk_mesh_number(-10), Flux_mesh_number(-10),
   Poisson_dof_number(0), Source_fct_pt(0),
   Exact_solution_fct_pt(0)
    {}

  /// Destructor (empty)
  ~GenericPoissonProblem() {}

  /// Doc the solution.
  virtual void doc_solution(DocInfo& doc_info) const;

  /// \short Create flux elements on boundary b of bulk mesh with the flux
  /// given by the value of the pointer (store them in the surface mesh).
  void set_neumann_boundary(const unsigned &b,
                            PoissonFluxFctPt const prescribed_flux_pt);

  /// Pin nodes on boundary b of bulk mesh.
  void set_dirichlet_boundary(const unsigned &b,
                              const DirichletFctPt condition_fct_pt)
  {
   // Pin all nodes on boundary b
   unsigned n_node = bulk_mesh_pt()->nboundary_node(b);
   for(unsigned nd=0;nd<n_node;nd++)
    {
     // We assume that Poisson elements only have one dof, pin it:
     bulk_mesh_pt()->boundary_node_pt(b,nd)->pin(Poisson_dof_number);
    }

   // Store function pointer to compute values on boundary b
   std::pair<unsigned, DirichletFctPt> dcb;
   dcb.first = b;
   dcb.second = condition_fct_pt;
   Dirichlet_conditions.push_back(dcb);
  }

  /// Finish building the problem once everything has been set.
  void build();

  /// Update dirichlet conditions before Newton solve.
  virtual void actions_before_newton_solve()
  {
   update_dirichlet_conditions();
  }


  void update_dirichlet_conditions()
  {
   for(unsigned i=0; i<Dirichlet_conditions.size(); i++)
    {
     unsigned b = Dirichlet_conditions[i].first;
     DirichletFctPt fct_pt = Dirichlet_conditions[i].second;

     for(unsigned nd=0; nd< bulk_mesh_pt()->nboundary_node(b); nd++)
      {
       Node* node_pt = bulk_mesh_pt()->boundary_node_pt(b,nd);

       // Get the value at this point
       Vector<double> x(node_pt->ndim());
       node_pt->position(x);
       Vector<double> value(1);
       fct_pt(x,value);

       // Assign it
       node_pt->set_value(Poisson_dof_number, value[0]);
      }
    }
  }

  /// Get the norm of the error (if exact solution has been set and problem
  /// has been solved).
  double get_error_norm() const;

  // Access functions
  // ============================================================

  /// Set function for Source_fct_pt.
  void set_source_fct_pt
   (typename ELEMENT::PoissonSourceFctPt const source_fct_pt)
  {Source_fct_pt = source_fct_pt;}

  /// Const access function for Source_fct_pt.
  typename ELEMENT::PoissonSourceFctPt source_fct_pt() const
   {
    if(Source_fct_pt == 0)
     {
      std::ostringstream error_msg;
      error_msg << "Source function pointer for this problem "
                << "has not been set yet.";
      throw OomphLibError(error_msg.str(),
                          "",
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
                         "GenericPoissonProblem::set_bulk_mesh",
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
                         "GenericPoissonProblem::bulk_mesh_pt",
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
                         "GenericPoissonProblem::flux_mesh_pt",
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
  typename ELEMENT::PoissonSourceFctPt Source_fct_pt;

  /// \short Storage for a pointer to a function giving the exact solution
  /// (for validation).
  FiniteElement::SteadyExactSolutionFctPt Exact_solution_fct_pt;

  /// \short Storage a list of Dirichlet boundaries and pointers to the
  /// functions which define their values.
  Vector<std::pair<unsigned, DirichletFctPt> > Dirichlet_conditions;

  /// Inaccessible copy constructor
  GenericPoissonProblem(const GenericPoissonProblem& dummy)
   {BrokenCopy::broken_copy("GenericPoissonProblem");}

  /// Inaccessible assignment operator
  void operator=(const GenericPoissonProblem& dummy)
   {BrokenCopy::broken_assign("GenericPoissonProblem");}

 };

 // =================================================================
 /// Finish off building the problem (once everything has been set).
 // =================================================================
 template<class ELEMENT>
  void GenericPoissonProblem<ELEMENT>::
  build()
  {
#ifdef PARANOID
   for(unsigned i=0; i<Dirichlet_conditions.size(); i++)
    {
     if(Dirichlet_conditions[i].second == 0)
      {
       std::ostringstream error_msg;
       error_msg << "Dirichlet function pointer number "
                 << i << " is unassigned.";
       throw OomphLibError(error_msg.str(),
                           "",
                           OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

   // Loop over the bulk elements to set up element-specific things
   unsigned n_element = bulk_mesh_pt()->nelement();
   for(unsigned i=0;i<n_element;i++)
    {
     // Upcast from GeneralisedElement to the present element
     ELEMENT *el_pt = dynamic_cast<ELEMENT*>(bulk_mesh_pt()->element_pt(i));

     // Set the source function pointer
     el_pt->source_fct_pt() = Source_fct_pt;
    }

   // Set the values on Dirichlet boundaries in case they are needed before
   // any Newton solves take place (and so we pass self tests).
   update_dirichlet_conditions();

   // Combine the meshes into a global mesh
   build_global_mesh();

   // Set up equation numbering scheme
   std::cout <<"Poisson number of equations: " << assign_eqn_numbers() << std::endl;
  }

 // =====================================================================
 /// Create Poisson flux elements on the b-th boundary of the bulk mesh,
 /// add the elements to the flux mesh and set the prescribed flux pointer.
 //=======================================================================
 template<class ELEMENT>
  void GenericPoissonProblem<ELEMENT>::set_neumann_boundary
  (const unsigned &b, PoissonFluxFctPt const prescribed_flux_pt)
  {
   // If we don't have a flux mesh yet then make one
   if(Flux_mesh_number == -10)
    {
     Flux_mesh_number = add_sub_mesh(new Mesh);
    }

   // Loop over the bulk elements adjacent to boundary b
   unsigned n_element = bulk_mesh_pt()->nboundary_element(b);
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk element that is adjacent to boundary b
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>
      (bulk_mesh_pt()->boundary_element_pt(b,e));

     // What is the index of the face of the bulk element e on bondary b
     int face_index = bulk_mesh_pt()->face_index_at_boundary(b,e);

     // Build the corresponding prescribed-flux element
     PoissonFluxElement<ELEMENT>* flux_element_pt = new
      PoissonFluxElement<ELEMENT>(bulk_elem_pt,face_index);

     // Add the prescribed-flux element to the surface mesh
     flux_mesh_pt()->add_element_pt(flux_element_pt);

     // Set the prescribed flux on this element
     flux_element_pt->flux_fct_pt() = prescribed_flux_pt;

    } //end of loop over bulk elements adjacent to boundary b

  } // end of create_flux_elements


 // =================================================================
 /// Doc the solution.
 // =================================================================
 template<class ELEMENT>
  void GenericPoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
  const
  {

   std::ofstream some_file;
   char filename[100];

   // Number of plot points
   unsigned npts;
   npts=5;

   // Output solution
   sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
           doc_info.number());
   some_file.open(filename);
   bulk_mesh_pt()->output(some_file,npts);
   some_file.close();

   // If we have an exact solution then use it:
   if(exact_solution_fct_pt() != 0)
    {

     // Output exact solution
     sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
             doc_info.number());
     some_file.open(filename);
     bulk_mesh_pt()->output_fct(some_file, npts, exact_solution_fct_pt());
     some_file.close();


     // Doc error and return of the square of the L2 error
     double error,norm;
     sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
             doc_info.number());
     some_file.open(filename);
     bulk_mesh_pt()->compute_error(some_file, exact_solution_fct_pt(),
                                   error, norm);
     some_file.close();

     // Doc L2 error and norm of solution
     std::cout << "\nNorm of error   : " << sqrt(error) << std::endl;
     std::cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;

    }
   else
    {
     std::cout << "No exact soution pointer given." << std::endl;
    }
  }

 // =================================================================
 /// Get the norm of the error (requires exact solution fct_pt). Used for
 /// testing purposes.
 // =================================================================
 template<class ELEMENT>
  double GenericPoissonProblem<ELEMENT>::get_error_norm() const
  {
   double error, norm;
   std::ofstream dummy_file;

   if(exact_solution_fct_pt() != 0)
    {
     bulk_mesh_pt()->compute_error(dummy_file, exact_solution_fct_pt(),
                                   error, norm);
    }
   else
    {
     throw OomphLibError("No exact solution set.",
                         "GenericPoissonProblem::get_error_norm",
                         OOMPH_EXCEPTION_LOCATION);
    }
   return sqrt(error);
  }


} // End of oomph namespace

#endif
