#include "generic_poisson_problem.h"

namespace oomph
{

// =================================================================
/// Finish off building the problem (once everything has been set).
// =================================================================
void GenericPoissonProblem::
build(Vector<Mesh*>& bulk_mesh_pts)
{
  // Call the underlying build to deal with adding meshes
  MyProblem::build(bulk_mesh_pts);

  // Check that all the dirichlet conditions have exactly one condition
  // assigned.
#ifdef PARANOID
  for(unsigned i=0; i<Dirichlet_conditions.size(); i++)
    {
      DirichletCondition* d = &Dirichlet_conditions[i];
      if((d->dc_fct_pt == 0) && (d->dc_vector_pt == 0))
        {
          std::ostringstream error_msg;
          error_msg << "Dirichlet function pointer number "
                    << i << " has no condition assigned.";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      else if((d->dc_fct_pt != 0) && (d->dc_vector_pt != 0))
        {
          std::ostringstream error_msg;
          error_msg << "Dirichlet function pointer number "
                    << i << " has two conditions assigned.";
          throw OomphLibError(error_msg.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
    }
#endif


  // Pin Dirichlet boundaries
  for(unsigned i=0, n=Dirichlet_conditions.size(); i<n; i++)
    {
      unsigned b = Dirichlet_conditions[i].boundary;
      Mesh* mesh_pt = Dirichlet_conditions[i].mesh_pt;

      for(unsigned nd=0, nnd=mesh_pt->nboundary_node(b); nd<nnd; nd++)
        {
          Node* nd_pt = mesh_pt->boundary_node_pt(b,nd);
          nd_pt->pin(Poisson_dof_number);
        }
    }


  // If there are no Dirichlet conditions on a given mesh then we need to
  // pin some other value to zero to avoid a singular Jacobian.
  // Loop over all meshes in problem
  for(unsigned msh=0, nmsh=bulk_mesh_pts.size(); msh<nmsh; msh++)
    {

      // Look for this mesh in the Dirichlet conditions
      bool flag = false;
      for(unsigned i=0, n=Dirichlet_conditions.size(); i<n; i++)
        {
          if(Dirichlet_conditions[i].mesh_pt == bulk_mesh_pts[msh])
            {
              flag = true;
              break;
            }
        }

      // If we didn't find any dirichlet conditions then pin something.
      if(!flag)
        {
          Node* pinned_node_pt = bulk_mesh_pts[msh]->get_some_non_boundary_node();
          pinned_node_pt->pin(poisson_dof_number());
          pinned_node_pt->set_value(poisson_dof_number(), 0.0);
        }
    }



  // Loop over the bulk elements to set up element-specific things
  for(unsigned msh=0, nmsh=bulk_mesh_pts.size(); msh<nmsh; msh++)
    {
      Mesh* mesh_pt = bulk_mesh_pts[msh];
      for(unsigned i=0; i<mesh_pt->nelement(); i++)
        {
          // Upcast from GeneralisedElement to the present element
          TFPoissonEquations *el_pt =
            checked_dynamic_cast<TFPoissonEquations*>(mesh_pt->element_pt(i));

          // Set the source function pointer
          el_pt->source_fct_pt() = Source_fct_pt;
        }
    }


  // Build flux meshes
  for(unsigned i=0; i<Neumann_conditions.size(); i++)
    {
      build_neumann_boundary_mesh(Neumann_conditions[i].mesh_pt,
                                  Neumann_conditions[i].boundary,
                                  Neumann_conditions[i].prescribed_flux_pt);
    }

  // Combine the meshes into a global mesh
  build_global_mesh();
  mesh_pt()->setup_boundary_element_info();

  // Set the values on Dirichlet boundaries in case they are needed before
  // any Newton solves take place (and so we pass self tests).
  update_dirichlet_conditions();

  // Always use CG with amg w/ poisson settings (it's very good).
  linear_solver_pt() = Factories::linear_solver_factory("cg", "cr", 1e-8,
                                                        200, true);
  checked_dynamic_cast<IterativeLinearSolver*>(linear_solver_pt())
    ->preconditioner_pt() = Factories::preconditioner_factory("poisson-amg");

  // Poisson problems are linear so we can reuse Jacobians
  enable_jacobian_reuse();

  // And we can say so to skip a residual evaluation
  Problem_is_nonlinear = false;

  // Set up equation numbering scheme
  std::cout << "Poisson number of equations: " << assign_eqn_numbers() << std::endl;
}



// =================================================================
/// Doc the solution.
// =================================================================
void GenericPoissonProblem::doc_solution_additional(std::ofstream &soln_file)
  const
{

  using namespace StringConversion;

  // Number of plot points
  unsigned npts = 2;

  // Output solution with specified number of plot points per element
  mesh_pt()->output(soln_file, npts);

  // If we have an exact solution then use it:
  if(exact_solution_fct_pt() != 0)
    {
      // Output exact solution
      std::ofstream exact_file((Doc_info.directory() + "/exact_soln"
                                + to_string(Doc_info.number()) + ".dat").c_str());
      exact_file.precision(Output_precision);
      exact_file.close();


      // Doc error and return of the square of the L2 error
      double error,norm;
      std::ofstream error_file((Doc_info.directory() + "/error"
                                + to_string(Doc_info.number()) + ".dat").c_str());
      error_file.precision(Output_precision);

      // Loop over all bulk meshes in problem
      for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
        {
          // Skip non-bulk meshes
          FiniteElement* ele_pt = mesh_pt(msh)->finite_element_pt(0);
          if(ele_pt->nodal_dimension() != ele_pt->dim()) continue;

          mesh_pt(msh)->compute_error(error_file, exact_solution_fct_pt(),
                                      error, norm);
        }

      error_file.close();

      // Doc L2 error and norm of solution
      std::cout << "\nNorm of error   : " << sqrt(error) << std::endl;
      std::cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;

    }
  else
    {
      std::cout << "No exact solution pointer given." << std::endl;
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
      // Loop over all bulk meshes in problem
      for(unsigned msh=0, nmsh=nsub_mesh(); msh<nmsh; msh++)
        {

          // Skip non-bulk meshes
          FiniteElement* ele_pt = mesh_pt(msh)->finite_element_pt(0);
          if(ele_pt->nodal_dimension() != ele_pt->dim()) continue;

          mesh_pt(msh)->compute_error(dummy_file, exact_solution_fct_pt(),
                                      error, norm);
        }
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
    for(unsigned i=0; i<Dirichlet_conditions.size(); i++)
      {
        // Extract values
        Mesh* mesh_pt = Dirichlet_conditions[i].mesh_pt;
        unsigned b = Dirichlet_conditions[i].boundary;
        DirichletFctPt fct_pt = Dirichlet_conditions[i].dc_fct_pt;
        const DoubleVector* vec_pt = Dirichlet_conditions[i].dc_vector_pt;

        // Assign values from function
        if(fct_pt != 0)
          {
            for(unsigned nd=0; nd<mesh_pt->nboundary_node(b); nd++)
              {
                Node* node_pt = mesh_pt->boundary_node_pt(b,nd);

                // Get the value at this point
                Vector<double> x(node_pt->ndim());
                node_pt->position(x);
                Vector<double> value(1);
                fct_pt(x,value);

                // Assign it
                node_pt->set_value(Poisson_dof_number, value[0]);
              }
          }
        // Or copy from vector
        else if(vec_pt != 0)
          {

            for(unsigned nd=0; nd<mesh_pt->nboundary_node(b); nd++)
              {
                Node* node_pt = mesh_pt->boundary_node_pt(b,nd);
                node_pt->set_value(Poisson_dof_number, (*vec_pt)[nd]);
              }
          }
        else
          {
            std::string err = "No dirichlet condition assigned!";
            throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                                OOMPH_CURRENT_FUNCTION);
          }
      }
  }

}
