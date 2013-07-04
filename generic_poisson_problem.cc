#include "generic_poisson_problem.h"

namespace oomph
{

// =================================================================
/// Finish off building the problem (once everything has been set).
// =================================================================
void GenericPoissonProblem::
build()
{
#ifdef PARANOID
  for(unsigned i=0; i<Dirichlet_function_conditions.size(); i++)
    {
      if(Dirichlet_function_conditions[i].second == 0)
        {
          std::ostringstream error_msg;
          error_msg << "Dirichlet function pointer number "
                    << i << " is unassigned.";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }

  for(unsigned i=0; i<Dirichlet_vector_conditions.size(); i++)
    {
      if(Dirichlet_vector_conditions[i].second == 0)
        {
          std::ostringstream error_msg;
          error_msg << "Dirichlet vector pointer number "
                    << i << " is unassigned.";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }
#endif

  // Pin Dirichlet boundaries set by functions
  for(unsigned i=0, n=Dirichlet_function_conditions.size(); i<n; i++)
    {
      unsigned b = Dirichlet_function_conditions[i].first;
      for(unsigned nd=0, nnd=bulk_mesh_pt()->nboundary_node(b); nd<nnd; nd++)
        {
          Node* nd_pt = bulk_mesh_pt()->boundary_node_pt(b,nd);
          nd_pt->pin(Poisson_dof_number);
        }
    }

  // Pin Dirichlet boundaries set by vectors
  for(unsigned i=0, n=Dirichlet_vector_conditions.size(); i<n; i++)
    {
      unsigned b = Dirichlet_vector_conditions[i].first;
      for(unsigned nd=0, nnd=bulk_mesh_pt()->nboundary_node(b); nd<nnd; nd++)
        {
          Node* nd_pt = bulk_mesh_pt()->boundary_node_pt(b,nd);
          nd_pt->pin(Poisson_dof_number);
        }
    }


  // Loop over the bulk elements to set up element-specific things
  unsigned n_element = bulk_mesh_pt()->nelement();
  for(unsigned i=0;i<n_element;i++)
    {
      // Upcast from GeneralisedElement to the present element
      TFPoissonEquations *el_pt =
        checked_dynamic_cast<TFPoissonEquations*>(bulk_mesh_pt()->element_pt(i));

      // Set the source function pointer
      el_pt->source_fct_pt() = Source_fct_pt;
    }

  // Set the values on Dirichlet boundaries in case they are needed before
  // any Newton solves take place (and so we pass self tests).
  update_dirichlet_conditions();

  // Combine the meshes into a global mesh
  build_global_mesh();

  // Always use CG with amg w/ poisson settings (it's very good).
  linear_solver_pt() = Factories::linear_solver_factory("cg");
  checked_dynamic_cast<IterativeLinearSolver*>(linear_solver_pt())
    ->preconditioner_pt() = Factories::preconditioner_factory("poisson-amg");

  // Set up equation numbering scheme
  std::cout << "Poisson number of equations: " << assign_eqn_numbers() << std::endl;
}

// =====================================================================
/// Create Poisson flux elements on the b-th boundary of the bulk mesh,
/// add the elements to the flux mesh and set the prescribed flux pointer.
//=======================================================================



// =================================================================
/// Doc the solution.
// =================================================================
void GenericPoissonProblem::doc_solution(DocInfo& doc_info)
  const
{

  using namespace StringConversion;

  // Number of plot points
  unsigned npts;
  npts=5;

  // Output solution
  std::ofstream soln_file((doc_info.directory() + "/soln"
                           + to_string(doc_info.number()) + ".dat").c_str());
  bulk_mesh_pt()->output(soln_file,npts);
  soln_file.close();

  // If we have an exact solution then use it:
  if(exact_solution_fct_pt() != 0)
    {

      // Output exact solution
      std::ofstream exact_file((doc_info.directory() + "/exact_soln"
                                + to_string(doc_info.number()) + ".dat").c_str());
      bulk_mesh_pt()->output_fct(exact_file, npts, exact_solution_fct_pt());
      exact_file.close();


      // Doc error and return of the square of the L2 error
      double error,norm;
      std::ofstream error_file((doc_info.directory() + "/error"
                                + to_string(doc_info.number()) + ".dat").c_str());
      bulk_mesh_pt()->compute_error(error_file, exact_solution_fct_pt(),
                                    error, norm);
      error_file.close();

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
double GenericPoissonProblem::get_error_norm() const
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
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
  return sqrt(error);
}

  /// (Re-)apply the Dirichlet conditions on the boundaries which have
  /// been set as Dirichlet. Note: if values differ where two boundaries
  /// overlap then priority is 1) vector conditions > function
  /// conditions, 2) higher boundary number conditions > lower boundary
  /// number conditions. Might want to change this...
  void GenericPoissonProblem::update_dirichlet_conditions()
  {
    // Conditions assigned by function pointer
    for(unsigned i=0; i<Dirichlet_function_conditions.size(); i++)
      {
        unsigned b = Dirichlet_function_conditions[i].first;
        DirichletFctPt fct_pt = Dirichlet_function_conditions[i].second;

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

    // Conditions assigned by vector
    for(unsigned i=0; i<Dirichlet_vector_conditions.size(); i++)
      {
        const unsigned b = Dirichlet_vector_conditions[i].first;
        const DoubleVector* const vec_pt = Dirichlet_vector_conditions[i].second;

        for(unsigned nd=0; nd< bulk_mesh_pt()->nboundary_node(b); nd++)
          {
            Node* node_pt = bulk_mesh_pt()->boundary_node_pt(b,nd);
            node_pt->set_value(Poisson_dof_number, (*vec_pt)[nd]);
          }
      }
  }

}
